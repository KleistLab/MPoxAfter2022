import numpy as np
import pandas as pd
import os
import sys
import time

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
from parameters import parameter_sets
parameters = parameter_sets[args.model]

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


vacc1_t = berlin_mpox_all["first_dosis"].values
vacc2_t = berlin_mpox_all["second_dosis"].values


# filter until end date
berlin_mpox_all = berlin_mpox_all[berlin_mpox_all['start_date_week'] <= parameters["end_date"]]
data = berlin_mpox_all["reported_cases"].values

if parameters["imports"]:
	imports_t = berlin_mpox_all["imports_reported"].values
else:
	imports_t = np.zeros_like(vacc_t)



T = len(data) # cw1 is happening between [0,1)
# load data
data_directory = os.path.join(repo_root, '../../inputs/'+args.model)
file_path = os.path.join(data_directory, args.run+'.csv')
df = pd.read_csv(file_path)

result = []
result_directory = os.path.join(repo_root, '../../results/'+args.model+"/"+args.run+'/')
for i in range(args.index_start, args.index_end):
	# set seed to index
	print(i)
	sys.stdout.flush()
	beta = df.loc[i, "beta"]
	diag = df.loc[i, "diag"]
	seed = df.loc[i, "seed"]
	print(seed, beta, diag)
	np.random.seed(seed)

	model = Mpox_endemic_model(degrees_real, counts_bc.copy(), vacc1_t, vacc2_t, parameters)
	model.initialize(diag = diag)
	model.initialize_imports(imports_t)
	model.initialize_reaction_dict(diag = diag)

	model.update_SSA_propensities(
		beta = beta, 
		diag = diag,
		bool_vec = [True, True, True, True, True, True])

	start = time.time()

	model.simulate(
		beta = beta, 
		diag = diag,  
		T = T)
	end = time.time()
	print("Time elapsed:")
	print(end - start)

	# calculate distance to data
	tmp_df = pd.DataFrame(model.result_df, columns = ['time', 'event', 'degree', 'inf_id'])
	
	#file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')

	if args.final_run:
		# save everything
		#tmp_df = model.result_df
		file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')
		#tmp_df = pd.DataFrame(model.result_df, columns = ['time', 'event', 'degree'])
		tmp_df.to_csv(file_path, index = False)
	else:
		if model.successful_finished:
			diagnosis = tmp_df[(tmp_df["event"] == 8)|(tmp_df["event"] == 12)].reset_index()
			event_counts, _ = np.histogram(diagnosis["time"], bins = np.arange(0, T + 1))
			#distance = np.linalg.norm(data - event_counts)
			distance2 = np.linalg.norm(np.cumsum(data) - np.cumsum(event_counts))
			distance = np.linalg.norm(data - event_counts)
		else:
			distance = 10**6
			distance2 = 10**6
		# add distance to dataframe!
		result.append([i, distance, distance2])
		print("Distance:")
		print(distance)

	#tmp_df.to_csv(file_path, index = False)

if args.final_run == False:
	#result_directory = os.path.join(repo_root, '../../results/ABC_SMC/'+args.run+'/')
	df = pd.DataFrame(result, columns=["index", "distance", "distance_cum"])
	file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')
	df.to_csv(file_path, index=False)

current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
tracemalloc.stop()





