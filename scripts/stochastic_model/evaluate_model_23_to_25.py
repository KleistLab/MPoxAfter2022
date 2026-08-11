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

parser.add_argument('index', type=int)
parser.add_argument('reps', type=int)
parser.add_argument('model', type=str)

args = parser.parse_args()

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

# load data
file_path = os.path.join(data_directory, 'berlin_mpox_all.csv')
berlin_mpox_all = pd.read_csv(file_path)


vacc1_t = berlin_mpox_all["first_dosis"].values
vacc2_t = berlin_mpox_all["second_dosis"].values

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
data_directory = os.path.join(repo_root, '../../inputs/model2022')
file_path = os.path.join(data_directory, 'final_run.csv')
df = pd.read_csv(file_path)

result = []
result_directory = os.path.join(repo_root, '../../results/'+args.model+'/results23_25/')


beta = df.loc[args.index, "beta"]
diag = df.loc[args.index, "diag"]
seed1 = df.loc[args.index, "seed"]

print(seed1, beta, diag)
np.random.seed(seed1)

model = Mpox_endemic_model(degrees_real, counts_bc.copy(), vacc1_t, vacc2_t, parameters2022)
model.initialize(diag = diag)
model.initialize_imports(imports_t)
model.initialize_reaction_dict(diag = diag)

model.update_SSA_propensities(
	beta = beta, 
	diag = diag,
	bool_vec = [True, True, True, True, True, True])


model.simulate(
	beta = beta, 
	diag = diag,  
	T = T1)


for i in range(args.reps):
	# set seed to index
	print(i)
	sys.stdout.flush()
	start = time.time()
	seed2 = i
	#model_tmp = model
	model_tmp = copy.deepcopy(model)
	# update parameters and simulate with fresh seed
	model_tmp.return_rate = parameters["return_rate"]
	model_tmp.dpsi = parameters["dpsi"]
	model_tmp.ve_I = parameters["ve_I"]
	model_tmp.initialize_imports(imports_t2)
	model_tmp.initialize_reaction_dict(diag = parameters["diag_factor"] * diag)

	model_tmp.update_SSA_propensities(
		beta = beta, 
		diag = parameters["diag_factor"] * diag,
		bool_vec = [True, True, True, True, True, True])

	print("Start simulation")
	np.random.seed(seed2)
	model_tmp.simulate(
		beta = beta, 
		diag = parameters["diag_factor"] * diag,  
		T = T2)
	end = time.time()
	print("Time elapsed:")
	print(end - start)

	# calculate distance to data
	tmp_df = pd.DataFrame(model_tmp.result_df, columns = ['time', 'event', 'degree', 'inf_id'])
	
	#file_path = os.path.join(result_directory, 'result_'+str(i)+'.csv')
	
	# save everything
	#tmp_df = model.result_df
	file_path = os.path.join(result_directory, 'result_'+str(args.index)+"_"+str(i)+'.csv')
	#tmp_df = pd.DataFrame(model.result_df, columns = ['time', 'event', 'degree'])
	tmp_df.to_csv(file_path, index = False)
	


current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
tracemalloc.stop()





