import numpy as np
import pandas as pd
import os
import sys
import time
import copy

import tracemalloc
tracemalloc.start()


import argparse

from full_model import Mpox_endemic_model


parser = argparse.ArgumentParser()

parser.add_argument('index_start', type=int)
parser.add_argument('index_end', type=int)
parser.add_argument('model', type=str)
parser.add_argument('run', type=str)
parser.add_argument('--final_run', action='store_true')

args = parser.parse_args()

# load fixed parameters
# load fixed parameters
from parameters import parameter_sets
parameters2022 = parameter_sets["model2022"]
parameters = parameter_sets[args.model]

repo_root = os.path.dirname(os.path.abspath(__file__))
data_directory = os.path.join(repo_root, '../../data')
file_path = os.path.join(data_directory, 'mpox_2022_paper', 'population.csv')
population = pd.read_csv(file_path)
degrees_real = population["SEX_PARTNERS_MALE_ANAL_UNCODED"]/12
n = len(degrees_real)

population_reduced = population[population["MPX_VERHALTEN_REDUCED"] == 1]
_, counts_bc = np.unique(population_reduced["SEX_PARTNERS_MALE_ANAL_UNCODED"], return_counts=True)

if parameters["bc_off"]:
	# turn off behaviour change but keep track of the people that should have changed behaviour
	parameters2022["brf"] = 1.0

# load data
file_path = os.path.join(data_directory, 'berlin_mpox_all.csv')
berlin_mpox_all = pd.read_csv(file_path)

vacc1_t = berlin_mpox_all["first_dosis"].values
vacc2_t = berlin_mpox_all["second_dosis"].values

try:
	if parameters["vacc_off"]:
		vacc1_t = np.zeros_like(vacc1_t)
		vacc2_t = np.zeros_like(vacc2_t)
except:
	pass


t_label = pd.to_datetime(berlin_mpox_all["start_date_week"], format="%Y-%m-%d")
data = berlin_mpox_all["reported_cases"].values
imports_t2 = berlin_mpox_all[berlin_mpox_all['start_date_week'] > parameters2022["end_date"]]["imports_reported"].values

# filter until end date
berlin_mpox_all = berlin_mpox_all[berlin_mpox_all['start_date_week'] <= parameters2022["end_date"]]
imports_t = berlin_mpox_all["imports_reported"].values

if parameters["imports"] == False:
	imports_t2 = np.zeros_like(imports_t2)


T1 = max(t_label[t_label<=parameters2022["end_date"]].index) # fist simulate until end of 2022
T2 = len(data) # then simulate for the rest, with different seed and maybe different parameters
# load data


# load data
data_directory = os.path.join(repo_root, '../../inputs/'+args.model)
file_path = os.path.join(data_directory, args.run+'.csv')
df = pd.read_csv(file_path)
result = []
result_directory = os.path.join(repo_root, '../../results/'+args.model+"/"+args.run+'/')


if (df.loc[args.index_start, "beta"] != df.loc[args.index_start+1, "beta"]) or args.final_run == True:
	print("Different parametrization for 2022. Need to simulate 2022 in each iteration.")
	simulate22 = True
else:
	print("Same parametrization for 2022. Only one simulation for 2022.")
	simulate22 = False
	beta = df.loc[args.index_start, "beta"]
	diag = df.loc[args.index_start, "diag"]
	seed = df.loc[args.index_start, "seed"]
	# simulate 2022 model
	np.random.seed(seed)
	model = Mpox_endemic_model(degrees_real, counts_bc.copy(), vacc1_t, vacc2_t, parameters2022)
	model.initialize(diag = diag)
	model.initialize_imports(imports_t)
	model.initialize_reaction_dict(diag = diag)

	print("Start simulating 2022.")

	model.update_SSA_propensities(
		beta = beta, 
		diag = diag,
		bool_vec = [True, True, True, True, True, True])

	model.simulate(
		beta = beta, 
		diag = diag,  
		T = T1)

	rng_state = np.random.get_state()



for i in range(args.index_start, args.index_end):
	# set seed to index
	print(i)
	sys.stdout.flush()
	beta = df.loc[i, "beta"]
	diag = df.loc[i, "diag"]
	seed = df.loc[i, "seed"]

	diag2 = diag
	try:
		dpsi = df.loc[i, "dpsi"]
	except:
		dpsi = 1
	try:
		diag_factor = df.loc[i, "diag_factor"]
		diag2 = diag * diag_factor
	except:
		pass
	try:
		diag2 = df.loc[i, "diag2"]
	except:
		pass
	try:
		beta2 = df.loc[i, "beta2"]
	except:
		beta2 = beta

	try:
		diag_2nd_factor = df.loc[i, "diag_2nd_factor"]
	except:
		diag_2nd_factor = 1

	try:
		infectiousness_2nd_factor = df.loc[i, "infectiousness_2nd_factor"]
	except:
		infectiousness_2nd_factor = 1

	try:
		long_mpox_prob = df.loc[i, "long_mpox_prob"]
	except:
		long_mpox_prob = 0	
		

	if beta != beta2:
		print(seed, beta, diag, dpsi, beta2, diag2)

	else:
		print(seed, beta, diag, dpsi, diag_2nd_factor, infectiousness_2nd_factor)

	np.random.seed(seed)

	if simulate22:
		model_tmp = Mpox_endemic_model(degrees_real, counts_bc.copy(), vacc1_t, vacc2_t, parameters2022)
		model_tmp.initialize(diag = diag)
		model_tmp.initialize_imports(imports_t)
		model_tmp.initialize_reaction_dict(diag = diag)

		model_tmp.update_SSA_propensities(
			beta = beta, 
			diag = diag,
			bool_vec = [True, True, True, True, True, True])

		model_tmp.simulate(
			beta = beta, 
			diag = diag,  
			T = T1)
	else:
		model_tmp = copy.deepcopy(model)


	try:
		np.random.seed(df.loc[i, "seed2"])
	except:
		if simulate22 == False:
			np.random.set_state(rng_state)

	# update parameters and continue simulation
	model_tmp.dpsi = dpsi
	model_tmp.ve_I = parameters["ve_I"]
	model_tmp.initialize_imports(imports_t2)
	model_tmp.diag_2nd_factor = diag_2nd_factor
	model_tmp.infectiousness_2nd_factor = infectiousness_2nd_factor
	model_tmp.long_mpox_prob = long_mpox_prob


	model_tmp.initialize_reaction_dict(diag = diag2)

	start = time.time()
	model_tmp.update_SSA_propensities(
		beta = beta2, 
		diag = diag2,
		bool_vec = [True, True, True, True, True, True])

	model_tmp.simulate(
		beta = beta2, 
		diag = diag2,  
		T = T2)
	end = time.time()

	print("\nTime elapsed:", end - start)

	# calculate distance to data
	tmp_df = pd.DataFrame(model_tmp.result_df, columns = ['time', 'event', 'degree', 'inf_id'])
	
	#file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')

	if args.final_run:
		# save everything
		file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')
		tmp_df.to_csv(file_path, index = False)
	else:
		if model_tmp.successful_finished:
			diagnosis = tmp_df[(tmp_df["event"] == 8)|(tmp_df["event"] == 12)].reset_index()
			event_counts, _ = np.histogram(diagnosis["time"], bins = np.arange(0, T2 + 1))
			distance2 = np.linalg.norm(np.cumsum(data) - np.cumsum(event_counts))
			distance = np.linalg.norm(data - event_counts)
		else:
			distance = 10**6
			distance2 = 10**6
		# add distance to dataframe!
		result.append([i, distance, distance2])
		print("Distance:", distance)


if args.final_run == False:
	#result_directory = os.path.join(repo_root, '../../results/ABC_SMC/'+args.run+'/')
	df = pd.DataFrame(result, columns=["index", "distance", "distance_cum"])
	file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')
	df.to_csv(file_path, index=False)

current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
tracemalloc.stop()





