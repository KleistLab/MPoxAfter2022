import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from lmfit import Parameters
import sys
import time

import tracemalloc
tracemalloc.start()

import cProfile

# parameters
model_name = "2ndinf"
sims = 100

repo_root = os.path.dirname(os.path.abspath(__file__))
data_directory = os.path.join(repo_root, '../../inputs')
result_directory = os.path.join(repo_root, '../../results/model_die_out_global/')
file_path = os.path.join(data_directory, 'model_die_out_global', 'die_out_local.csv')
die_out_local = pd.read_csv(file_path)
die_out_prob = die_out_local[model_name].values


from model_global import Global_die_out_model

number_nodes = 44
reporting_countries_0 = 7
c = [0.134752,-0.013194,-1.489550]#[0.085894,-0.023508,0] #[0.055283,0,0]	# param from poisson regression

np.random.rand(1)

result = []

for i in range(sims):

	model = Global_die_out_model(number_nodes, die_out_prob, reporting_countries_0, c)
	model.run_model()
	result.append(model.rc_t)


df = pd.DataFrame(result)
file_path = os.path.join(result_directory, 'result_'+model_name+'.csv')
df.to_csv(file_path, index=False)

