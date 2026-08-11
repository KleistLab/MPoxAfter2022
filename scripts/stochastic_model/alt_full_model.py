import numpy as np
#import pandas as pd
import heapq as hq
from scipy.stats import erlang
#from scipy.integrate import ode
#from lmfit import minimize, Parameters

class Mpox_endemic_model:

	def __init__(self, degrees, scale_E, scale_I, vacc_t):
		#pass
		values, counts = np.unique(degrees, return_counts = True)

		self.P_k = counts
		self.p_k = counts/np.sum(counts)
		self.corr = (values * counts/np.sum(counts))/np.mean(degrees)
		self.k = values
		self.vacc_t = vacc_t
		#self.imports_t = imports_t
		#self.rate_E = 1/scale_E
		#self.rate_I = 1/scale_I
		self.scale_E = scale_E
		self.scale_I = scale_I
		self.n = len(degrees)
		#self.T = T
		self.time_queue = []
		self.t = 0
		self.rates = np.zeros(7)

		# S, SBR, E1, I1, D, RI, RV, S2I, S2V, E2, I2, CBR, CI, CI2, CIV, CD, CD2, CV compartments 
		# each compartment for each degree
		self.N = 26
		self.M = self.corr.shape[0]
		self.return_rate = 0.03 # estimated from first mpox paper, hopefully not important
		#self.behaviour_reduction_factor = 10
		self.initial_behaviour_reduction = 0.22 # estimate from first mpox paper
		self.I0 = 53 # source, Bartel import paper
		self.waning = 0.23 # source, nat comms predicting...

	def initialize_reaction_dict(self, diag, dpsi):
		tmp = np.zeros_like(self.k)
		self.SSA_reaction_dict = {
		#		name					eventID		from 		to			recalculate propensities 							triggered queue event 	rate
			0: ["infection1S", 			1,			self.S, 	self.E1, 	[True, False, False, False, False, False, False],	0,						tmp],
			1: ["infection1SBR", 		2,			self.SBR, 	self.E1, 	[False, True, False, False, False, False, True],	0,						tmp],
			2: ["waning_infection", 	3,			self.RI, 	self.S2I, 	[False, False, True, False, True, False, False],	None,					tmp],
			3: ["waning_vaccination", 	4,			self.RV, 	self.S2V, 	[False, False, False, True, False, True, False],	None,					tmp],
			4: ["reinfectionI", 		5,			self.S2I, 	self.E2, 	[False, False, False, False, True, False, False],	4,						tmp],
			5: ["reinfectionV", 		6,			self.S2V, 	self.E2, 	[False, False, False, False, False, True, False],	4,						tmp],
			6: ["behaviour_return", 	7,			self.SBR, 	self.S, 	[True, True, False, False, False, False, True],		None,					tmp]
			}
		self.scheduled_reaction_dict = {
		#		name					eventID		from 		to		 	recalculate propensities							triggered queue event	rate 				a	
			0: ["incubation1", 			8,			self.E1, 	self.I1,	[True, True, False, False, True, True, False],		[1, 2],					self.scale_E,		3 	],
			1: ["recovery1", 			9,			self.I1, 	self.RI, 	[True, True, True, False, True, True, False],		None,					self.scale_I, 		3 	],
			2: ["diagnosis1", 			10,			self.I1, 	self.D, 	[True, True, False, False, True, True, False],		[3],					1/diag,		 		1 	],
			3: ["recoveryD", 			11,			self.D, 	self.RI, 	[False, False, True, False, False, False, False],	None,					self.scale_I, 		1 	],
			4: ["incubation2", 			12,			self.E2, 	self.I2, 	[True, True, False, False, True, True, False],		[5, 6],					self.scale_E, 		3 	],
			5: ["recovery2", 			13,			self.I2, 	self.RI, 	[True, True, True, False, True, True, False],		None,					self.scale_I/dpsi,	3 	],
			6: ["diagnosis2", 			14,			self.I2, 	self.D, 	[True, True, False, False, True, True, False],		[3],					1/diag,		 		1 	],
			7: ["vaccinationS", 		15,			self.S, 	self.RV, 	[True, False, False, True, False, False, False],	None, 					None,				None],
			8: ["vaccinationSBR",		16,			self.SBR,	self.RV, 	[False, True, False, True, False, False, True],		None, 					None, 				None]
			
			}

	def update_SSA_propensities(self, brf, beta, ve, ve_I, bool_vec):
		infection_pressure = np.sum(self.corr * (self.I1 + self.I2)/self.P_k)

		if bool_vec[0] or bool_vec[1] or bool_vec[4] or bool_vec[5]:
			infection_pressure = np.sum(self.corr * (self.I1 + self.I2)/self.P_k)
		
		if bool_vec[0]:
			# infection1S
			self.SSA_reaction_dict[0][6] = self.S * self.k * beta * infection_pressure
			self.rates[0] = np.sum(self.SSA_reaction_dict[0][6])
		if bool_vec[1]:
			# infection1SBR
			self.SSA_reaction_dict[1][6] = self.SBR * brf * self.k * beta * infection_pressure
			self.rates[1] = np.sum(self.SSA_reaction_dict[1][6])
		if bool_vec[2]:
			# waning infection
			self.SSA_reaction_dict[2][6] = self.RI * self.waning
			self.rates[2] = np.sum(self.SSA_reaction_dict[2][6])
		if bool_vec[3]:
			# waning ve
			self.SSA_reaction_dict[3][6] = self.RV * self.waning
			self.rates[3] = np.sum(self.SSA_reaction_dict[3][6])
		if bool_vec[4]:
			# reinfectionI
			self.SSA_reaction_dict[4][6] = self.S2I * (1 - ve_I) * self.k * beta * infection_pressure
			self.rates[4] = np.sum(self.SSA_reaction_dict[4][6])
		if bool_vec[5]:
			# reinfectionV
			self.SSA_reaction_dict[5][6] = self.S2V * (1 - ve) * self.k * beta * infection_pressure
			self.rates[5] = np.sum(self.SSA_reaction_dict[5][6])
		if bool_vec[6]:	
			# behaviour return
			self.SSA_reaction_dict[6][6] = self.SBR * self.return_rate
			self.rates[6] = np.sum(self.SSA_reaction_dict[6][6])

		# update rates
		#self.rates = np.array([np.sum(values[3]) for values in self.SSA_reaction_dict.values()])


	def initialize(self, diag):
		#E1, I1 = np.zeros_like(self.ks, dtype = int), np.zeros_like(self.ks, dtype = int)
		# distribute initially infected/imports uniformly over all groups
		# half of the initially infected/imports should be E, the other hald I
		#self.E1 = np.random.multinomial(np.ceil(self.I0/2), np.ones_like(self.k)/len(self.k)).astype(int)
		#self.I1 = np.random.multinomial(np.floor(self.I0/2), np.ones_like(self.k)/len(self.k)).astype(int)

		# no infections in k = 0 group
		p_k0 = self.P_k/(self.n - self.P_k[0])
		p_k0[0] = 0
		self.E1 = np.random.multinomial(np.ceil(self.I0/2), p_k0).astype(int)
		self.I1 = np.random.multinomial(np.floor(self.I0/2), p_k0).astype(int)

		# all other agents are susceptible at start
		S = self.P_k - self.E1 - self.I1

		# a percentage of initially S is behaviour reduced
		self.SBR = np.round(self.initial_behaviour_reduction * S).astype(int)

		self.S = np.maximum(S - self.SBR, 0).astype(int)

		# all other compartments are emtpy at the start
		self.D = np.zeros_like(S, dtype = int)
		self.RI = np.zeros_like(S, dtype = int)
		self.RV = np.zeros_like(S, dtype = int)
		self.S2I = np.zeros_like(S, dtype = int)
		self.S2V = np.zeros_like(S, dtype = int)
		self.E2 = np.zeros_like(S, dtype = int)
		self.I2 = np.zeros_like(S, dtype = int)

		# incubation period, diagnosis delay and recovery period is sampled explicitly and modelled with queue
		for i in range(len(self.E1)):
			n = self.E1[i]
			while n > 0:
				# sample incubation period 
				hq.heappush(self.time_queue, (erlang.rvs(a = 3, scale = self.scale_E), 0, i))
				n -= 1

		for i in range(len(self.I1)):
			n = self.I1[i]
			while n > 0:
				# sample incubation period 
				incubation_period = erlang.rvs(a = 3, scale = self.scale_I)
				diagnosis_delay = np.random.exponential(1/diag)

				if incubation_period < diagnosis_delay:
					# no diagnosis
					hq.heappush(self.time_queue, (incubation_period, 1, i))
				else:
					# diagnosis and then recovery				
					hq.heappush(self.time_queue, (diagnosis_delay, 2, i))
					#hq.heappush(self.time_queue, (incubation_period, 3, i))
				n -= 1


		# result dataframe with all events
		#self.result_df = pd.DataFrame(columns=['time', 'event', 'degree'])
		self.result_df = []
		
		# vaccination logic: degree class is randomly assigned upon vaccination
		#tmp_vacc = vacc_t + np.random.rand(len(vacc_t))
		share0 = self.P_k[0]/self.n
		# remove zero contacts
		self.S[0] = 0
		self.SBR[0] = 0

		for i,count in enumerate(self.vacc_t):
			while count > 0:
				if np.random.rand() > share0:
					tmp_time = i + np.random.rand()
					if np.random.rand() < self.initial_behaviour_reduction:
						tmp_k = np.random.choice(len(self.k), p=self.SBR/(np.sum(self.SBR)))
						hq.heappush(self.time_queue, (tmp_time, 8, tmp_k))
					else:
						tmp_k = np.random.choice(len(self.k), p=self.S/(np.sum(self.S)))
						hq.heappush(self.time_queue, (tmp_time, 7, tmp_k))
				count -= 1

	def simulate(self, brf, beta, diag, ve, ve_I, dpsi, T):
		# simulate system until time T
		#self.initialize(brf, beta, diag, waning_I, ve, ve_I, dpsi)
		print(np.sum(self.SBR))
		print(np.sum(self.E1 + self.I1 + self.S + self.D + self.RI + self.RV + self.RI + self.S2I + self.S2V + self.E2 + self.I2 + self.SBR))
		#print("Number scheduled events:")
		#print(len(self.time_queue))
		marker = 1
		while self.t < T:
			#break
			#if self.t > marker:
			#	print(self.t)
			#	marker = np.ceil(self.t)
			# sample next SSA event
			r0 = np.sum(self.rates)
			if abs(r0) < 10**-6:
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
				index = np.searchsorted(np.cumsum(self.rates), u, side="right")
				# select degree within event
				if index > 0:
					u -= np.cumsum(self.rates)[index - 1]
				index_k = np.searchsorted(np.cumsum(self.SSA_reaction_dict[index][6]), u, side="right")
				# execute event
				event = self.SSA_reaction_dict[index]
				event[2][index_k] -= 1
				event[3][index_k] += 1
				# check if event has to be added to queue 
				if event[5] is not None:
					queue_event = self.scheduled_reaction_dict[event[5]]
					# add event to queue 
					hq.heappush(self.time_queue, (self.t + erlang.rvs(a = queue_event[7], scale = queue_event[6]), event[5], index_k))

				
			else:
				# upadate time
				self.t = t_queue
				event = self.scheduled_reaction_dict[self.time_queue[0][1]]
				index_k = self.time_queue[0][2]
				# execute event
				if self.time_queue[0][1] > 6: # vaccination
					if event[2][index_k] > 0:
						event[2][index_k] -= 1
						event[3][index_k] += 1
					else:
						# remove event from queue and move on
						hq.heappop(self.time_queue)
						continue

				else:
					event[2][index_k] -= 1
					event[3][index_k] += 1

				# remove event from queue 
				hq.heappop(self.time_queue)

				# check if event has to be added to queue 
				if event[5] is not None:
					tmp_queue = []
					for index_event in event[5]:
						queue_event = self.scheduled_reaction_dict[index_event]
						hq.heappush(tmp_queue, (self.t + erlang.rvs(a = queue_event[7], scale = queue_event[6]), index_event, index_k))
					# select only first event for actual event queue
					hq.heappush(self.time_queue, tmp_queue[0])


			# add update to dataframe
			#self.result_df.loc[len(self.result_df)] = [self.t, event[0], index_k]
			#self.result_df.append([self.t, event[0], index_k])
			self.result_df.append([self.t, event[1], index_k])

			# update propensities
			self.update_SSA_propensities(brf, beta, ve, ve_I, event[4])




