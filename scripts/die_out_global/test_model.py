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


from model_global import Global_die_out_model

number_nodes = 10
die_out_prob = np.ones(100) * 0.1
reporting_countries_0 = 10
c = 0.055283	# param from poisson regression
I0 = number_nodes

model = Global_die_out_model(number_nodes, die_out_prob, reporting_countries_0, c, I0)

model.run_model()

print(model.result)
plt.plot(model.rc_t)
plt.show()

'''
from full_model import Mpox_endemic_model

# load fixed parameters
from parameters import parameter_sets
parameters = parameter_sets["model2022"]#["model1"]

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

vacc1_t = berlin_mpox_all["first_dosis"].values
vacc2_t = berlin_mpox_all["second_dosis"].values
data = berlin_mpox_all["reported_cases"].values
imports_t2 = berlin_mpox_all[berlin_mpox_all['start_date_week'] > parameters["end_date"]]["imports_reported"].values

# filter until end date
berlin_mpox_all_tmp = berlin_mpox_all[berlin_mpox_all['start_date_week'] <= parameters["end_date"]]

imports_t = berlin_mpox_all_tmp["imports_reported"].values

if parameters["imports"] == False:
	imports_t2 = np.zeros_like(imports_t2)



#T = len(data)
T1 = max(t_label[t_label<=parameters["end_date"]].index) # fist simulate until end of 2022
T2 = len(data) # then simulate for the rest, with different seed and maybe different parameters

beta = 1.7817066581863767
diag = 0.1042623850605849


#ve = 1.0
#ve_I = ve


# set seed globally
np.random.seed(9378)
model = Mpox_endemic_model(degrees_real, counts_bc, vacc1_t, vacc2_t, parameters)
#model = Mpox_endemic_model(degrees_real, vacc_t, imports_t, parameters)
model.initialize(diag = diag)
model.initialize_imports(imports_t)
model.initialize_reaction_dict(diag = diag)

print(np.sum(model.S+model.SBR+model.E1+model.I1))



model.update_SSA_propensities(
	beta = beta, 
	diag = diag,
	bool_vec = [True, True, True, True, True, True])

start = time.time()

model.simulate(
	beta = beta, 
	diag = diag,  
	T = T1)


model.dpsi = 1.204082
diag_factor = 0.363265

model.initialize_imports(imports_t2)
model.initialize_reaction_dict(diag = diag_factor * diag)

model.demographic_influx = 10

model.update_SSA_propensities(
	beta = beta, 
	diag = diag,
	bool_vec = [True, True, True, True, True, True])

model.simulate(
	beta = beta, 
	diag = diag_factor * diag,  
	T = T2)
end = time.time()

end = time.time()
print("Time elapsed")
print(end - start)


print(model.S)
print(model.SBR)
print(model.E1)
print(model.I1)
print(model.D)
print(model.RI)
print(model.RV)
print(model.S2I)
print(model.S2V)
print(model.E2)
print(model.I2)

result_directory = os.path.join(repo_root, '../../results/stochastic_test')
file_path = os.path.join(result_directory, 'test_T1.csv')
tmp_df = pd.DataFrame(model.result_df, columns = ['time', 'event', 'degree', 'inf_id'])
tmp_df.to_csv(file_path, index = False)
#model.result_df.to_csv(file_path, index = False)

diagnosis = tmp_df[(tmp_df["event"] == 8)|(tmp_df["event"] == 12)].reset_index()
event_counts, _ = np.histogram(diagnosis["time"], bins = np.arange(0, T2 + 1))
#distance = np.linalg.norm(np.cumsum(data) - np.cumsum(event_counts))
distance = np.linalg.norm(data - event_counts)
print("Distance cases:", distance)

distance = np.linalg.norm(np.cumsum(data) - np.cumsum(event_counts))
print("Distance cum cases:", distance)

print(np.sum(model.E1 + model.I1 + model.S + model.D + model.RI + model.RV + model.S2I + model.S2V + model.E2 + model.I2 + model.SBR))

current, peak = tracemalloc.get_traced_memory()
print(f"Current memory usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
tracemalloc.stop()

diagnosis = tmp_df[tmp_df["event"].isin([8,12])].reset_index(drop = True)
diag_counts, _ = np.histogram(diagnosis["time"], bins = np.arange(0, len(data) +1))

fig, ax = plt.subplots(1, 2, figsize=(9, 4))

ax[0].plot(t_label[_[:-1]], diag_counts, color = "coral", label = "Stoch")
ax[1].plot(t_label[_[:-1]], np.cumsum(diag_counts), color = "coral", label = "Stoch")

ax[0].scatter(berlin_mpox_all[berlin_mpox_all["reported_cases"]>0]["start_date_week"], berlin_mpox_all[berlin_mpox_all["reported_cases"]>0]["reported_cases"], color = "forestgreen", alpha = 0.3, label = "Berlin")
ax[0].set_ylabel("Reported cases")
ax[1].plot(t_label, np.cumsum(data), color = "forestgreen", label = "Data")
ax[1].set_ylabel("Cumulative cases")

plt.show()
'''


