from tqdm import tqdm
import time
import argparse
import numpy as np
import pandas as pd
import os
import pickle

def get_result_dict(t_max):
    result_dict = {}

    result_dict["diagnosis"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["infection"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["infectionS"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["infectionSBR"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["behaviour_return"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["incubation"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["recovery"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["firstinfection"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["reinfection"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["reinfectionV"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["waningI"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["waningV"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["vaccinations1"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["vaccinations2"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["vaccinations_failed"] = pd.DataFrame(columns = np.arange(0, t_max))
    
    result_dict["diagnosis_adjusted"] = pd.DataFrame(columns = np.arange(0, t_max))
    
    result_dict["infection_risk"] = pd.DataFrame(columns = np.arange(0, 60))
    result_dict["reinfection_risk"] = pd.DataFrame(columns = np.arange(0, 60))
    return(result_dict)

def add_sim_to_dict(result_df, event_name, eventIDs, result_dict):
    tmp_df = result_df[result_df["event"].isin(eventIDs)].reset_index(drop = True)
    event_counts, _ = np.histogram(tmp_df["time"], bins = np.arange(0, len(result_dict[event_name].columns)+1))
    new_row_df = pd.DataFrame([event_counts], columns=result_dict[event_name].columns)
    result_dict[event_name] = pd.concat([result_dict[event_name], new_row_df], ignore_index=True)
    return(result_dict)

def add_sim_to_dict_degreebased(result_df, event_name, eventIDs, result_dict):
    tmp_df = result_df[result_df["event"].isin(eventIDs)].reset_index(drop = True)
    _, counts = np.unique(tmp_df["degree"], return_counts=True)
    new_row_df = pd.DataFrame([counts], columns=result_dict[event_name].columns)
    result_dict[event_name] = pd.concat([result_dict[event_name], new_row_df], ignore_index=True)
    return(result_dict)

event_to_id = {}

event_to_id["diagnosis"] = [8,12]
event_to_id["infection"] = [1,2,3,4,101,102,103,104]
event_to_id["infectionS"] = [1,101]
event_to_id["infectionSBR"] = [2,102]
event_to_id["behaviour_return"] = [5]
event_to_id["incubation"] = [6,10]
event_to_id["recovery"] = [7, 11, 9]
event_to_id["firstinfection"] = [1, 2, 101, 102]
event_to_id["reinfection"] = [3, 103]
event_to_id["reinfectionV"] = [4, 104]
event_to_id["vaccinations1"] = [13, 14]
event_to_id["vaccinations2"] = [15]
event_to_id["vaccinations_failed"] = [113,114]

# year dates
year_start_dates = ["01-01-2022","01-01-2023","01-01-2024","01-01-2025","01-01-2026"]

# load data
repo_root = os.path.dirname(os.path.abspath(__file__))
data_directory = os.path.join(repo_root, '../data')
file_path = os.path.join(data_directory, 'berlin_mpox_all.csv')
berlin_mpox_all = pd.read_csv(file_path)
t_label = pd.to_datetime(berlin_mpox_all["start_date_week"], format="%Y-%m-%d")
data = berlin_mpox_all["reported_cases"].values

parser = argparse.ArgumentParser()
parser.add_argument('model', type=str)
args = parser.parse_args()

data_directory = os.path.join(repo_root, '../results')

file_path = os.path.join(data_directory, args.model, 'results23_25')
print(file_path)


csv_files = [f for f in os.listdir(file_path) if f.endswith('.csv')]
print("Files found:",len(csv_files))


result_dict_full = get_result_dict(len(data))

tmp = []

# Wrapping item_list with tqdm to create a progress bar
for file in tqdm(csv_files, desc="Processing files"):
    file_path_tmp = os.path.join(file_path, file)
    result_df = pd.read_csv(file_path_tmp)

    for event in event_to_id.keys():
        result_dict_full = add_sim_to_dict(result_df, event, event_to_id[event], result_dict_full)

    result_dict_full = add_sim_to_dict_degreebased(result_df, "infection_risk", event_to_id["firstinfection"], result_dict_full)
    result_dict_full = add_sim_to_dict_degreebased(result_df, "reinfection_risk", event_to_id["infection"], result_dict_full)

    # ugly code
    tmp.append(np.linalg.norm(data - result_dict_full["diagnosis"].iloc[-1].values))

    #for i in range(len(dates)-1):
    #inf_tmp = [row.sum() for _,row in result_dict_full["infection"].iloc[:, np.where((t_label >= dates[i]) & (t_label < dates[i+1]))[0]].iterrows()]
    #diag_tmp = [row.sum() for _,row in result_dict_full["diagnosis_adjusted"].iloc[:, np.where((t_label >= dates[i]) & (t_label < dates[i+1]))[0]].iterrows()]


result_dict_full["distances"] = np.array(tmp)

file_path = os.path.join(data_directory, args.model, 'result_dict_full.pkl')

with open(file_path, 'wb') as f:
    pickle.dump(result_dict_full, f)


