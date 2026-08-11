import numpy as np
#import pandas as pd
import heapq as hq
from scipy.stats import erlang
#from scipy.integrate import ode
#from lmfit import minimize, Parameters
import sys
from util import *

class Mpox_endemic_model:

	def __init__(self, degrees, counts_bc, vacc1_t, vacc2_t, parameters, demographic_influx = 0):
		#pass
		values, counts = np.unique(degrees, return_counts = True)

		# network specific
		self.P_k = counts
		self.k = values
		self.p_k = self.P_k/np.sum(self.P_k)
		self.corr = (self.k  * self.p_k)/np.mean(degrees)
		self.n = len(degrees)
		self.degrees = degrees

		# simulation specific
		self.time_queue = []
		self.t = 0
		self.rates = np.zeros(5)
		self.successful_finished = True
		self.epsilon = 10**-6
		self.max_number_updates = 10**7
		self.time_dependent_diagnosis = parameters["time_dependent_diagnosis"]
		self.infection_id = 1
		self.vaccination_id = -1
		self.imports = 0
		self.demographic_influx = demographic_influx
		self.import_check = False
		self.import_start = 0

		# fixed parameters
		self.return_rate = parameters["return_rate"]
		self.brf = parameters["brf"]
		#self.initial_behaviour_reduction = parameters["initial_behaviour_reduction"]
		self.I0 = parameters["I0"]
		#self.waning = parameters["waning"]
		self.scale_E = parameters["scale_E"]
		self.scale_I = parameters["scale_I"]
		self.dpsi = parameters["dpsi"]
		self.ve = parameters["ve"]
		self.ve_I = parameters["ve_I"]
		self.share0 = parameters["share0"]
		self.diag_2nd_factor = 1
		self.infectiousness_2nd_factor = 1
		self.long_mpox_prob = parameters["long_mpox_prob"]

		# outbreak specific
		self.vacc1_t = np.round((1-self.share0) * vacc1_t)
		self.vacc2_t = np.round((1-self.share0) * vacc2_t)
		#self.imports_t = imports_t
		self.SBR = counts_bc
		self.initial_behaviour_reduction = np.sum(counts_bc)/(self.n - self.P_k[0])

		# immune waning
		ts = np.arange(0, 52 * 4)
		self.ecdf1 = 1 - efficacy(ab_concentration(ts, parameters["x01"], parameters["f1"], parameters["delta_s"], parameters["delta_l"]), parameters["ic50"])
		self.ecdf2 = 1 - efficacy(ab_concentration(ts, parameters["x02"], parameters["f2"], parameters["delta_s"], parameters["delta_l"]), parameters["ic50"])
		self.ecdf_inf = 1 - efficacy(ab_concentration(ts, parameters["x0inf"], parameters["f2"], parameters["delta_s"], parameters["delta_l"]), parameters["ic50"])

		# bookkeeping
		self.importID = 18
		self.vaccinationIDs = [13,14,15]
		self.first_vaccinationIDs = [13,14]
		self.second_vaccID = 15
		self.infectionIDs = [1,2,3,4]
		self.waning_vaccinationID = 17
		self.waning_infectionID = 16
		self.recoverIDs = [7,9,11]
		self.diagnosisIDs = [8, 12]


	def initialize_reaction_dict(self, diag):
		tmp = np.zeros_like(self.k)
		self.SSA_reaction_dict = {
		#		name					eventID		from 		to			recalculate propensities 			triggered queue event 	rate 	
			1:	["infection1S", 		1,			self.S, 	self.E1, 	[True, False, False, False, False],	6,						tmp],
			2:	["infection1SBR", 		2,			self.SBR, 	self.E1, 	[False, True, False, False, True],	6,						tmp],
			3:	["reinfectionI", 		3,			self.S2I, 	self.E2, 	[False, False, True, False, False],	10,						tmp],
			4:	["reinfectionV", 		4,			self.S2V, 	self.E2, 	[False, False, False, True, False],	10,						tmp],
			5:	["behaviour_return", 	5,			self.SBR, 	self.S, 	[True, True, False, False, True],	None,					tmp],
			101:["infection1S", 		101,		self.S, 	self.I1, 	[True, False, False, False, False],	6,						tmp],
			102:["infection1SBR", 		102,		self.SBR, 	self.I1, 	[False, True, False, False, True],	6,						tmp],
			103:["reinfectionI", 		103,		self.S2I, 	self.I2, 	[False, False, True, False, False],	10,						tmp],
			104:["reinfectionV", 		104,		self.S2V, 	self.I2, 	[False, False, False, True, False],	10,						tmp],
			}
		self.scheduled_reaction_dict = {
		#		name					eventID		from 		to		 		recalculate propensities				triggered queue event	rate 							a	
			6: ["incubation1", 			6,			self.E1, 	self.I1,		[True, True, True, True, False],		[7, 8],					self.scale_E,					1	],
			7: ["recovery1", 			7,			self.I1, 	self.RI, 		[True, True, True, True, False],		None,					self.scale_I, 					5 	],
			8: ["diagnosis1", 			8,			self.I1, 	self.D, 		[True, True, True, True, False],		None,					1/diag,		 					1 	],
			9: ["recoveryD", 			9,			self.D, 	self.RI, 		[False, False, False, False, False],	None,					self.scale_I, 					5 	],
			10:["incubation2", 			10,			self.E2, 	self.I2, 		[True, True, True, True, False],		[11, 12],				self.scale_E, 					1 	],
			11:["recovery2", 			11,			self.I2, 	self.RI, 		[True, True, True, True, False],		None,					self.scale_I/self.dpsi,			5 	],
			12:["diagnosis2", 			12,			self.I2, 	self.D, 		[True, True, True, True, False],		None,					1/(self.diag_2nd_factor*diag),	1 	],
			13:["vaccination1S", 		13,			self.S, 	self.RV, 		[True, False, False, False, False],		None, 					None,							None],
			14:["vaccination1SBR",		14,			self.SBR,	self.RV, 		[False, True, False, False, True],		None, 					None, 							None],
			15:["vaccination2", 		15,			self.RV, 	self.RV, 		[True, False, False, False, False],		None, 					None,							None],
			16:["waning_infection", 	16,			self.RI, 	self.S2I, 		[False, False, True, False, False],		None,					None,							None],
			17:["waning_vaccination", 	17,			self.RV, 	self.S2V, 		[False, False, False, True, False],		None,					None,							None],
			18:["stimulate_import", 	18,			None, 		None, 			[True, True, True, True, False],		None,					None,							None],
			}

	def update_SSA_propensities(self, beta, diag, bool_vec):
		# if diagnosis is time dependent, perform update
		if self.time_dependent_diagnosis:
			tmp_diag = 1/(max(np.log10(np.sum(self.D) + self.epsilon), 1) * diag)
			self.scheduled_reaction_dict[2][6] = tmp_diag
			self.scheduled_reaction_dict[6][6] = tmp_diag

		if bool_vec[0] or bool_vec[1] or bool_vec[2] or bool_vec[3]:
			infection_pressure = np.sum(self.corr * (self.I1 + self.infectiousness_2nd_factor*self.I2)/self.P_k) + self.imports
		
		if bool_vec[0]:
			# infection1S
			self.SSA_reaction_dict[1][6] = self.S * self.k * beta * infection_pressure
			self.rates[0] = np.sum(self.SSA_reaction_dict[1][6])
		if bool_vec[1]:
			# infection1SBR
			self.SSA_reaction_dict[2][6] = self.SBR * self.brf * self.k * beta * infection_pressure
			self.rates[1] = np.sum(self.SSA_reaction_dict[2][6])
		if bool_vec[2]:
			# reinfectionI
			self.SSA_reaction_dict[3][6] = self.S2I * (1 - self.ve_I) * self.k * beta * infection_pressure
			self.rates[2] = np.sum(self.SSA_reaction_dict[3][6])
		if bool_vec[3]:
			# reinfectionV
			self.SSA_reaction_dict[4][6] = self.S2V * (1 - self.ve) * self.k * beta * infection_pressure
			self.rates[3] = np.sum(self.SSA_reaction_dict[4][6])
		if bool_vec[4]:	
			# behaviour return
			if self.t >= 10:
				self.SSA_reaction_dict[5][6] = self.SBR * self.return_rate
			else:
				self.SSA_reaction_dict[5][6] = np.zeros_like(self.SBR)
			self.rates[4] = np.sum(self.SSA_reaction_dict[5][6])

		# update rates
		#self.rates = np.array([np.sum(values[3]) for values in self.SSA_reaction_dict.values()])


	def initialize(self, diag):
		# distribute initially infected/imports uniformly over all groups
		# no infections in k = 0 group
		self.p_k0 = self.P_k/(self.n - self.P_k[0])
		self.p_k0[0] = 0
		self.E1 = np.random.multinomial(np.ceil(self.I0/2), self.p_k0).astype(int)
		self.I1 = np.random.multinomial(np.floor(self.I0/2), self.p_k0).astype(int)

		# all other agents are susceptible at start
		#S = self.P_k - self.E1 - self.I1

		# a percentage of initially S is behaviour reduced
		#self.SBR = np.round(self.initial_behaviour_reduction * S).astype(int)

		self.S = np.maximum(self.P_k - self.E1 - self.I1 - self.SBR, 0).astype(int)

		# all other compartments are emtpy at the start
		self.D = np.zeros_like(self.S, dtype = int)
		self.RI = np.zeros_like(self.S, dtype = int)
		self.RV = np.zeros_like(self.S, dtype = int)
		self.S2I = np.zeros_like(self.S, dtype = int)
		self.S2V = np.zeros_like(self.S, dtype = int)
		self.E2 = np.zeros_like(self.S, dtype = int)
		self.I2 = np.zeros_like(self.S, dtype = int)

		# result dataframe with all events
		self.result_df = []

		# incubation period, diagnosis delay and recovery period is sampled explicitly and modelled with queue
		for i in range(len(self.E1)):
			n = self.E1[i]
			while n > 0:
				inf_id = self.infection_id
				# sample incubation period 
				hq.heappush(self.time_queue, (erlang.rvs(a = 3, scale = self.scale_E), 6, i, inf_id))
				# add infection event to result
				self.result_df.append([self.t, 1, i, inf_id])
				self.infection_id += 1
				n -= 1

		for i in range(len(self.I1)):
			n = self.I1[i]
			while n > 0:
				inf_id = self.infection_id
				# sample recovery period 
				recovery_period = erlang.rvs(a = 5, scale = self.scale_I)
				diagnosis_delay = np.random.exponential(1/diag)
				# add infection and incubation event to queue
				self.result_df.append([self.t, 1, i, inf_id])
				self.result_df.append([self.t, 6, i, inf_id])

				if recovery_period < diagnosis_delay:
					# no diagnosis
					hq.heappush(self.time_queue, (recovery_period, self.recoverIDs[0], i, inf_id))
				else:
					# diagnosis and then recovery				
					hq.heappush(self.time_queue, (diagnosis_delay, self.diagnosisIDs[0], i, inf_id))
					hq.heappush(self.time_queue, (recovery_period, self.recoverIDs[1], i, inf_id))
				self.infection_id += 1
				n -= 1

		self.S[0] = 0
		self.SBR[0] = 0

		###############################################################################################################################################################################
		# vaccination logic: degree class is randomly assigned and we hope for the best
		# remove zero contacts

		# sample first dose
		# second dose will be sampled when vaccination occurs
		for i,count in enumerate(self.vacc1_t):
			while count > 0:
				tmp_time = i + np.random.rand()
				# decide between normal behaviour and reduced behaviour
				if np.random.rand() < self.initial_behaviour_reduction:
					tmp_k = np.random.choice(len(self.k), p=self.SBR/(np.sum(self.SBR)))
					eventID = self.vaccinationIDs[1]
				else:
					tmp_k = np.random.choice(len(self.k), p=self.S/(np.sum(self.S)))
					eventID = self.vaccinationIDs[0]
				hq.heappush(self.time_queue, (tmp_time, eventID, tmp_k, self.vaccination_id))
				count -= 1
				self.vaccination_id -= 1

	def initialize_imports(self, imports_t):
		###############################################################################################################################################################################
		# schedule time to increase force of infection to stimulate infection due to import
		for i,count in enumerate(imports_t):
			while count > 0:
				tmp_time = self.t + i + np.random.rand() * 0.05 # make sure that events do not happen exactly at the same time
				hq.heappush(self.time_queue, (tmp_time, self.importID, None, 0))
				count -= 1


	def simulate(self, beta, diag, T):
		# simulate system until time T
		marker = 1
		while self.t < T:
			if self.t > marker:
			#	print(self.t)
				sys.stdout.write(f"\rCurrent simulation time: {self.t} weeks")
				sys.stdout.flush()
				marker = np.ceil(self.t)
				# perform demographic change if possible
				if self.demographic_influx > 0:
					tmp_S = np.random.multinomial(self.demographic_influx, self.p_k0).astype(int)
					self.S += tmp_S
					self.P_k += tmp_S
					self.p_k = self.P_k/np.sum(self.P_k)
					self.corr = (self.k  * self.p_k)/np.mean(self.degrees)
					print("Sum:",np.sum(self.S + self.SBR+ self.E1 + self.E2 + self.I1 + self.I2 + self.RI + self.RV + self.S2I + self.S2V + self.D))
					self.update_SSA_propensities(beta, diag, [True, True, True, True, True, True])

			# break after max number of updates (RAM footprint becomes too large otherwise)
			if len(self.result_df) >= self.max_number_updates:
				self.successful_finished = False
				break

			if self.import_check:
				# stop if 50 follow up infections occured or if no infected is left
				if self.infection_id - self.import_start > 50 or np.sum(self.E1 + self.I1 + self.E2 + self.I2) == 0:
					break

			# sample next SSA event
			r0 = np.sum(self.rates)
			# make sure there is event to sample
			if r0 < self.epsilon:
				t_ssa = T + 1
			else:
				t_ssa = self.t + np.random.exponential(1/r0)

			# next event on the queue 
			if len(self.time_queue) > 0:
				t_queue = self.time_queue[0][0]
			else:
				t_queue = T + 1

			if min(t_ssa, t_queue) >= T:
				break

			if t_ssa < t_queue:
				# update time
				self.t = t_ssa
				# execute SSA event
				# select event type
				u = r0 * np.random.rand()
				index = np.searchsorted(np.cumsum(self.rates), u, side="right") + 1 # because dictionary starts at 1
				event = self.SSA_reaction_dict[index]
				eventID = event[1]
				# select degree within event
				if index > 1:
					u -= np.cumsum(self.rates)[index - 2]
				index_k = np.searchsorted(np.cumsum(self.SSA_reaction_dict[index][6]), u, side="right")


				if eventID in self.infectionIDs:
					# check if it is an import
					if self.imports > 0:
						self.imports -= 1
						inf_id = self.infection_id
						self.infection_id += 1
						eventID = event[1] + 100
						# select peron with more than 5 expected partners (consistent with 2022 project)
						#while (index_k < 5) and event[2][index_k] > 0:
						#	index_k = np.random.choice(np.arange(5,60))

						# remove person from compartment
						event[2][index_k] -= 1
						# directly add it to infectious compartment
						if event[1] < 3:
							self.I1[index_k] += 1
							diag_id = self.diagnosisIDs[0]
							rec_id = self.recoverIDs[1]
						else:
							self.I2[index_k] += 1
							diag_id = self.diagnosisIDs[1]
							rec_id = self.recoverIDs[1]

						# add time to recovery and time to diagnosis (100% diag rate)
						recovery_period = erlang.rvs(a = 5, scale = self.scale_I)
						diagnosis_delay = np.random.rand() * recovery_period
						hq.heappush(self.time_queue, (self.t + diagnosis_delay, diag_id, index_k, inf_id))
						hq.heappush(self.time_queue, (self.t + recovery_period, rec_id, index_k, inf_id))

						self.update_SSA_propensities(beta, diag, [True, True, True, True, False])
						self.result_df.append([self.t, eventID, index_k, inf_id])
						continue


				# execute event
				event[2][index_k] -= 1
				event[3][index_k] += 1
				# check if event has to be added to queue 
				if event[5] is not None:
					inf_id = self.infection_id
					queue_event = self.scheduled_reaction_dict[event[5]]
					# add event to queue
					if queue_event[7] == 1:
						hq.heappush(self.time_queue, (self.t + np.random.exponential(queue_event[6]), event[5], index_k, inf_id))
					else:
						hq.heappush(self.time_queue, (self.t + erlang.rvs(a = queue_event[7], scale = queue_event[6]), event[5], index_k, inf_id))
					self.infection_id += 1
				else:
					inf_id = 0


			##################################################################################################################################################################################################################				
			else:
				# upadate time
				self.t = t_queue
				event = self.scheduled_reaction_dict[self.time_queue[0][1]]
				eventID = event[1]
				index_k = self.time_queue[0][2]
				inf_id = self.time_queue[0][3]

				# remove event from queue 
				hq.heappop(self.time_queue)
				# execute event
				# vaccination exception
				if eventID in self.vaccinationIDs:
					# check if vaccination is possible or if agents already got infected
					if event[2][index_k] > 0:
						event[2][index_k] -= 1
						event[3][index_k] += 1

						# check if second vaccination is possible in 2-4 weeks
						if eventID in self.first_vaccinationIDs:
							t_index = int(np.floor(self.t)) + 2
							t_inf_max = min(t_index + 2, len(self.vacc2_t)-1)
							tmp_bool_waning = True
							while t_index <= t_inf_max:
								if (self.vacc2_t[t_index] > 0):
									# shedule 2nd vacc event
									hq.heappush(self.time_queue, (t_index + np.random.rand(), self.second_vaccID, index_k, inf_id))
									self.vacc2_t[t_index] -= 1
									tmp_bool_waning = False
									break
								t_index += 1

							if tmp_bool_waning:
								# no 2nd vacc, schedule waning
								hq.heappush(self.time_queue, (self.t + np.searchsorted(self.ecdf1, np.random.rand(), side="right"), self.waning_vaccinationID, index_k, inf_id))

						# for 2nd vaccine: schedule waning
						else:
							hq.heappush(self.time_queue, (self.t + np.searchsorted(self.ecdf2, np.random.rand(), side="right"), self.waning_vaccinationID, index_k, inf_id))

					else:
						# add failed vacc event to timeline
						self.result_df.append([self.t, eventID+100, index_k, inf_id])
						continue



				elif eventID == self.importID: # stimulate import
					self.imports += 1

				else:
					event[2][index_k] -= 1
					event[3][index_k] += 1


				# recovery exception
				if eventID in self.recoverIDs:
					# schedule waning
					hq.heappush(self.time_queue, (self.t + np.searchsorted(self.ecdf_inf, np.random.rand(), side="right"), self.waning_infectionID, index_k, inf_id))


				# check if event has to be added to queue 
				if event[5] is not None:
					tmp_queue = []
					for index_event in event[5]:
						queue_event = self.scheduled_reaction_dict[index_event]
						if queue_event[7] == 1:
							hq.heappush(tmp_queue, (self.t + np.random.exponential(queue_event[6]), index_event, index_k, inf_id))
						else:
							# check for long mpox 
							long_mpox = False
							if self.long_mpox_prob > 0:
								if np.random.rand() < self.long_mpox_prob:
									long_mpox = True
							if long_mpox:
								hq.heappush(tmp_queue, (self.t + 28 + (180-28) * np.random.rand(), index_event, index_k, inf_id))	# uniformly sample between (28--180) days
							else:
								hq.heappush(tmp_queue, (self.t + erlang.rvs(a = queue_event[7], scale = queue_event[6]), index_event, index_k, inf_id))
					# select only first event for actual event queue
					hq.heappush(self.time_queue, tmp_queue[0])
					if len(tmp_queue) == 2:
						# special case: if second event in tmp_queue is recovery, add it as well
						# this is now recovery from the diagnosed compartment (first tmp_queue event was diagnosis)
						if tmp_queue[1][1] in [7,11]:
							hq.heappush(self.time_queue, (tmp_queue[1][0], 9, tmp_queue[1][2], tmp_queue[1][3]))


			# add update to dataframe
			self.result_df.append([self.t, eventID, index_k, inf_id])

			# update propensities
			self.update_SSA_propensities(beta, diag, event[4])

	def execute_reaction(self, eventID, degree):
		if eventID <= 5 or eventID > 100:
			self.SSA_reaction_dict[eventID][2][degree] -= 1
			self.SSA_reaction_dict[eventID][3][degree] += 1
		else:
			self.scheduled_reaction_dict[eventID][2][degree] -= 1
			self.scheduled_reaction_dict[eventID][3][degree] += 1




