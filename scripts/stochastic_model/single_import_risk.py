import numpy as np
import pandas as pd
import os
import sys
import time
import copy

import heapq as hq

import tracemalloc
tracemalloc.start()

import argparse

from full_model import Mpox_endemic_model


parser = argparse.ArgumentParser()

#parser.add_argument('index', type=int)
parser.add_argument('T_end', type=int)
parser.add_argument('model', type=str)
args = parser.parse_args()

# load fixed parameters
from parameters import parameter_sets
parameters2022 = parameter_sets["model2022"]
parameters = parameter_sets[args.model]
parameters2022["I0"] = 0

repo_root = os.path.dirname(os.path.abspath(__file__))
data_directory = os.path.join(repo_root, '../../data')
file_path = os.path.join(data_directory, 'mpox_2022_paper', 'population.csv')
population = pd.read_csv(file_path)
degrees_real = population["SEX_PARTNERS_MALE_ANAL_UNCODED"]/12
n = len(degrees_real)

population_reduced = population[population["MPX_VERHALTEN_REDUCED"] == 1]
_, counts_bc = np.unique(population_reduced["SEX_PARTNERS_MALE_ANAL_UNCODED"], return_counts=True)


# load data
file_path = os.path.join(data_directory, 'berlin_mpox_all.csv')
berlin_mpox_all = pd.read_csv(file_path)


t_label = pd.to_datetime(berlin_mpox_all["start_date_week"], format="%Y-%m-%d")
data = berlin_mpox_all["reported_cases"].values

# filter until end date
berlin_mpox_all = berlin_mpox_all[berlin_mpox_all['start_date_week'] <= parameters2022["end_date"]]

T_end = args.T_end
T2 = len(data)

# result matrix 
res_mat = np.zeros((60, 100 * 10))

# load parameters
data_directory = os.path.join(repo_root, '../../inputs/'+args.model)
result_directory = os.path.join(repo_root, '../../results/'+args.model+'/final_run/')
file_path = os.path.join(data_directory, 'final_run.csv')
df = pd.read_csv(file_path)
result = []

for index in range(0,100):

	start = time.time()

	# load result file 
	file_name = "result_"+str(index) + ".csv"
	file_path = os.path.join(result_directory, file_name)
	result = pd.read_csv(file_path)

	# initialize model
	beta = df.loc[index, "beta"]
	diag = df.loc[index, "diag"]
	seed = df.loc[index, "seed"]
	diag2 = diag
	beta2 = beta
	try:
		dpsi = df.loc[i, "dpsi"]
	except:
		dpsi = 1
	try:
		infectiousness_2nd_factor = df.loc[i, "infectiousness_2nd_factor"]
	except:
		infectiousness_2nd_factor = 1

	model = Mpox_endemic_model(degrees_real, counts_bc.copy(), np.zeros(len(data)), np.zeros(len(data)), parameters2022)
	model.initialize(diag = diag)
	#model.initialize_imports(imports_t)
	model.initialize_reaction_dict(diag = diag)

	# replay model until T_end
	i = 0
	max_inf_id = 0
	max_vacc_id = 0
	print("Start replay simulation.")
	while result.loc[i, "time"] < T_end:
		eventID = result.loc[i, "event"]
		if eventID in [18, 113, 114]:
			i += 1
			continue
		degree = int(result.loc[i, "degree"])
		# execute event
		model.execute_reaction(eventID, degree)
		# save max inf ids that have been dealt with so far
		max_inf_id = max(max_inf_id, result.loc[i, "inf_id"])
		max_vacc_id = min(max_vacc_id, result.loc[i, "inf_id"])
		i += 1

	model.t = result.loc[i-1, "time"]
	print("Simulation replayed.")
	print("Remove currently infected and keep track of waning times.")
	# schedule immune waning reactions and perform recovery for currently infected/exposed nodes (so up until max_inf_id)
	# also keep vaccinations scheduled
	for k in range(i,len(result)):
		inf_id = result.loc[k, "inf_id"]
		if inf_id > max_inf_id:
			continue
		eventID = result.loc[k, "event"]
		if eventID in [16, 17]:
			if eventID == 17 and result.loc[k, "inf_id"] < max_vacc_id:
				continue
			# schedule waning
			hq.heappush(model.time_queue, (result.loc[k, "time"], eventID, int(result.loc[k, "degree"]), result.loc[k, "inf_id"]))
		elif eventID in [6,7,8,9,10,11,12]:
			# remove all currently infected
			degree = int(result.loc[k, "degree"])
			model.execute_reaction(eventID, degree)
		elif eventID in [13,14]:
			# schedule vaccinations
			hq.heappush(model.time_queue, (result.loc[k, "time"], eventID, int(result.loc[k, "degree"]), result.loc[k, "inf_id"]))

	# check if there are infected left
	for degree in range(0,61):
		model.RI[degree] += model.E1[degree] + model.I1[degree] + model.E2[degree] + model.I2[degree]
		model.E1[degree] = 0
		model.I1[degree] = 0
		model.E2[degree] = 0
		model.I2[degree] = 0

	# save the state of the model 
	model.dpsi = dpsi
	model.infectiousness_2nd_factor = infectiousness_2nd_factor
	model.infection_id = max_inf_id + 1
	model.import_check = True
	model.import_start = model.infection_id


	for degree in range(1,61):
		for rep in range(10):
			model_tmp = copy.deepcopy(model)
			# introduce one additional infection in degree class 
			if model_tmp.S[degree] > 0:
				model_tmp.S[degree] -= 1
				model_tmp.E1[degree] += 1
				triggered_event = 6
			elif model_tmp.SBR[degree] > 0:
				model_tmp.SBR[degree] -= 1
				model_tmp.E1[degree] += 1
				triggered_event = 6
			elif model_tmp.S2I[degree] > 0:
				model_tmp.S2I[degree] -= 1
				model_tmp.E2[degree] += 1
				triggered_event = 10
			elif model_tmp.S2V[degree] > 0:
				model_tmp.S2V[degree] -= 1
				model_tmp.E2[degree] += 1
				triggered_event = 10

			# add event to queue
			hq.heappush(model_tmp.time_queue, (model_tmp.t + np.random.exponential(1), triggered_event, degree, model_tmp.infection_id))
			model_tmp.infection_id += 1

			# update propensities
			model_tmp.initialize_reaction_dict(diag = diag2)
			model_tmp.update_SSA_propensities(
				beta = beta2, 
				diag = diag2,
				bool_vec = [True, True, True, True, True, True])

			model_tmp.simulate(
				beta = beta2, 
				diag = diag2,  
				T = T_end + 50)
			end = time.time()

			res_mat[degree-1, index*10 + rep] = model_tmp.infection_id - model.infection_id

	print("\nResult "+ str(index) + " done. Time: "+ str(end - start))

print(res_mat)


df = pd.DataFrame(res_mat)
file_path = os.path.join(result_directory, 'single_import_'+str(T_end)+'.csv')
df.to_csv(file_path, index=False)






