import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from pathlib import Path
import os

current_dir = Path().resolve()

#data_directory = current_dir.parent / 'data'
file_path = current_dir.parent / 'data/mpox_2022_paper/population.csv'
population = pd.read_csv(file_path)
degrees_real = population["SEX_PARTNERS_MALE_ANAL_UNCODED"]/12
n = len(degrees_real)

file_path = current_dir.parent / 'data/berlin_mpox_all.csv'
berlin_mpox_all = pd.read_csv(file_path)

vacc1_t = berlin_mpox_all["first_dosis"].values
vacc2_t = berlin_mpox_all["second_dosis"].values
imports_t = berlin_mpox_all["imports_reported"].values
data = berlin_mpox_all["reported_cases"].values#/len(degrees_real)
t = np.arange(len(data))
t_label = pd.to_datetime(berlin_mpox_all["start_date_week"], format="%Y-%m-%d")

# load accepted parameters
file_path = current_dir.parent / 'inputs/model2022/final_run.csv'
accepted_params = pd.read_csv(file_path)

# load behaviour change date
population_reduced = population[(population["MPX_VERHALTEN_REDUCED"] == 1)&(population["SEX_PARTNERS_MALE_ANAL_UNCODED"]>0)]
degree_bc, counts_bc = np.unique(population_reduced["SEX_PARTNERS_MALE_ANAL_UNCODED"], return_counts=True)

# calculate network statistics 
values, counts = np.unique(degrees_real, return_counts = True)

# network specific
P_k = counts
k = values
p_k = P_k/np.sum(P_k)
corr = (k  * p_k)/np.mean(degrees_real)

P_k0 = counts[1:]
p_k0 = P_k0/np.sum(P_k0)
K = k[k>0] * 12

#N = len(degrees_real)
N = len(degrees_real[degrees_real>0])
K_tot = np.sum(degrees_real * 12)
k0 = K_tot/N
K_tot_S = K_tot - np.sum(degree_bc * counts_bc)

n0 = len(degrees_real[degrees_real>0])

def get_result_dict(t_max):
    result_dict = {}

    result_dict["diagnosis"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["diagnosis1"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["diagnosis2"] = pd.DataFrame(columns = np.arange(0, t_max))
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
    result_dict["imports"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["imports_re"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["imports_rev"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["immune_gain"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["immune_loss"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["br_leaving"] = pd.DataFrame(columns = np.arange(0, t_max))

    result_dict["kr"] = pd.DataFrame(columns = np.arange(0, t_max))
    result_dict["krr"] = pd.DataFrame(columns = np.arange(0, t_max))
    
    result_dict["diagnosis_adjusted"] = pd.DataFrame(columns = np.arange(0, t_max))

    result_dict["die_out_prob"] = pd.DataFrame(columns = np.arange(0, t_max))

    result_dict["infection_risk"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["reinfection_risk"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["vaccination1_risk"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["vaccination2_risk"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["br_risk"] = pd.DataFrame(columns = np.arange(1, 61))

    result_dict["infection_risk22"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["reinfection_risk22"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["vaccination1_risk22"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["vaccination2_risk22"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["br_risk22"] = pd.DataFrame(columns = np.arange(1, 61))

    result_dict["infection_risk23"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["reinfection_risk23"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["vaccination1_risk23"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["vaccination2_risk23"] = pd.DataFrame(columns = np.arange(1, 61))
    result_dict["br_risk23"] = pd.DataFrame(columns = np.arange(1, 61))

    return(result_dict)

event_to_id = {}

event_to_id["diagnosis"] = [8,12]
event_to_id["diagnosis1"] = [8]
event_to_id["diagnosis2"] = [12]
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
event_to_id["imports"] = [101,102,103,104]
event_to_id["imports_re"] = [103]
event_to_id["imports_rev"] = [104]
event_to_id["immune_gain"] = [1,2,3,4,101,102,103,104,13,14]
event_to_id["immune_loss"] = [16,17]
event_to_id["waningI"] = [16]
event_to_id["waningV"] = [17]
event_to_id["br_leaving"] = [2,5,14]


def add_sim_to_dict(result_df, event_name, eventIDs, result_dict):
    tmp_df = result_df[result_df["event"].isin(eventIDs)].reset_index(drop = True)
    event_counts, _ = np.histogram(tmp_df["time"], bins = np.arange(0, len(result_dict[event_name].columns)+1))
    new_row_df = pd.DataFrame([event_counts], columns=result_dict[event_name].columns)
    result_dict[event_name] = pd.concat([result_dict[event_name], new_row_df], ignore_index=True)
    return(result_dict)

def add_sim_to_dict_degreebased(result_df, event_name, eventIDs, result_dict, date = 200):
    # filter for events before some time 
    result_df = result_df[result_df["time"] <= date]
    tmp_df = result_df[result_df["event"].isin(eventIDs)].reset_index(drop = True)

    degrees = result_dict[event_name].columns
    degree_counts = np.zeros(len(degrees), dtype=int)
    unique_degrees, counts = np.unique(tmp_df["degree"], return_counts=True)
    for u_degree, count in zip(unique_degrees, counts):
        degree_idx = degrees.get_loc(u_degree)
        degree_counts[degree_idx] = count
    new_row_df = pd.DataFrame([degree_counts], columns=degrees)
    result_dict[event_name] = pd.concat([result_dict[event_name], new_row_df], ignore_index=True)
    return(result_dict)

def add_sim_geometry(result_df, result_dict, bc, pos_off):
    #N = len(degrees_real)
    #K_tot = np.sum(degrees_real * 12)
    degrees, S0 = np.unique(degrees_real[degrees_real>0]*12, return_counts = True)
    if bc:
        S0 -= counts_bc
        K_tot_tmp = K_tot_S
    else:
        K_tot_tmp = K_tot

    # kr for frailty 
    if pos_off:
        pos_events = []
    else:
        pos_events = [16, 17, 5]   # -> +1 16,17 are immune waning, 5 is behaviour return
    neg_events = [1,2,3,4,101, 102, 103, 104, 13, 14]   # -> -1

    result_df["sign"] = np.select(
        [result_df["event"].isin(pos_events), result_df["event"].isin(neg_events)],
        [1, -1],
        default=0
    ).astype(int)

    result_df["degree_effective"] = result_df["degree"] * result_df["sign"]
    result_df["cum_change"] = result_df["sign"].cumsum()
    result_df["cum_deg_change"] = result_df["degree_effective"].cumsum()
    result_df["nS"] = np.sum(S0) + result_df["cum_change"]
    result_df["KS"] = K_tot_tmp + result_df["cum_deg_change"]
    result_df["kr"] = result_df["KS"] / result_df["nS"]

    frailty_int = pd.merge_asof(
        pd.DataFrame({"time": np.arange(0, len(data), dtype=float)}),
        result_df.sort_values("time")[["time", "kr"]],
        on="time",
        direction="backward"
    )["kr"].ffill().to_numpy()

    new_row_df = pd.DataFrame([frailty_int], columns=result_dict["kr"].columns)
    result_dict["kr"] = pd.concat([result_dict["kr"], new_row_df], ignore_index=True)

    # krr for interference
    dS = (
    result_df.pivot_table(index=result_df.index, columns="degree", values="sign",
                   aggfunc="sum", fill_value=0)
    )
    #dS = dS.reindex(result_df.index, fill_value=0) 
    dS = dS.reindex(columns=degrees, fill_value=0)
    S = dS.cumsum().add(S0, axis=1)
    
    # probability degree k node is S
    p_Sk = S.div(S0, axis=1)   
    # corr * p_Sk    
    theta_k = p_Sk.mul(corr[corr>0], axis=1)
    # sum over corr * p_Sk
    theta = theta_k.sum(axis=1)
    # sum_k P(k) * theta * r_k
    result_df["krr"] = result_df["kr"] * theta

    interference_int = pd.merge_asof(
        pd.DataFrame({"time": np.arange(0, len(data), dtype=float)}),
        result_df.sort_values("time")[["time", "krr"]],
        on="time",
        direction="backward"
    )["krr"].ffill().to_numpy()

    new_row_df = pd.DataFrame([interference_int], columns=result_dict["krr"].columns)
    result_dict["krr"] = pd.concat([result_dict["krr"], new_row_df], ignore_index=True)

    return(result_dict)


def get_die_out_prob(result_df, result_dict):
    pos_events = [1, 2, 3, 4, 101, 102, 103, 104]   # -> +1 all kind of infections
    neg_events = [7,9,11]   # -> -1 # all events that lead to end of infectiousness

    result_df["sign"] = np.select(
        [result_df["event"].isin(pos_events), result_df["event"].isin(neg_events)],
        [1, -1],
        default=0
    ).astype(int)
    result_df["I_active"] = result_df["sign"].cumsum()

    result_df['time_int'] = result_df['time'].astype(int)
    def determine_indicator(group):
        return 0 if (group['I_active'] == 0).any() else 1

    indicator_df = result_df.groupby('time_int', group_keys=False).apply(lambda group: pd.Series({
        'indicator': determine_indicator(group)
    }), include_groups=False).reset_index()

    new_row_df = pd.DataFrame([(1- indicator_df["indicator"].values)], columns=result_dict["die_out_prob"].columns)
    result_dict["die_out_prob"] = pd.concat([result_dict["die_out_prob"], new_row_df], ignore_index=True)

    return(result_dict)




def final_run_dict(model, bc = True, pos_off = False):

    file_tmp = 'results/'+model+'/final_run/'
    result_files_dict =  current_dir.parent / file_tmp
    csv_files = [f for f in os.listdir(result_files_dict) if f.endswith('.csv')]
    print(len(csv_files))

    model_dict = get_result_dict(len(data))

    for file in csv_files:
        if "result" in file:
            file_path = os.path.join(result_files_dict, file)
            result_df = pd.read_csv(file_path)
        else:
            continue
        
        for event in event_to_id.keys():
            model_dict = add_sim_to_dict(result_df, event, event_to_id[event], model_dict)

        model_dict = add_sim_to_dict_degreebased(result_df, "infection_risk", event_to_id["firstinfection"], model_dict)
        model_dict = add_sim_to_dict_degreebased(result_df, "reinfection_risk", event_to_id["infection"], model_dict)
        model_dict = add_sim_to_dict_degreebased(result_df, "vaccination1_risk", event_to_id["vaccinations1"], model_dict)
        model_dict = add_sim_to_dict_degreebased(result_df, "vaccination2_risk", event_to_id["vaccinations2"], model_dict)
        model_dict = add_sim_to_dict_degreebased(result_df, "br_risk", event_to_id["br_leaving"], model_dict)

        model_dict = add_sim_to_dict_degreebased(result_df, "infection_risk22", event_to_id["firstinfection"], model_dict, 26)
        model_dict = add_sim_to_dict_degreebased(result_df, "reinfection_risk22", event_to_id["infection"], model_dict, 26)
        model_dict = add_sim_to_dict_degreebased(result_df, "vaccination1_risk22", event_to_id["vaccinations1"], model_dict, 26)
        model_dict = add_sim_to_dict_degreebased(result_df, "vaccination2_risk22", event_to_id["vaccinations2"], model_dict, 26)
        model_dict = add_sim_to_dict_degreebased(result_df, "br_risk22", event_to_id["br_leaving"], model_dict, 26)

        model_dict = add_sim_to_dict_degreebased(result_df, "infection_risk23", event_to_id["firstinfection"], model_dict, 61)
        model_dict = add_sim_to_dict_degreebased(result_df, "reinfection_risk23", event_to_id["infection"], model_dict, 61)
        model_dict = add_sim_to_dict_degreebased(result_df, "vaccination1_risk23", event_to_id["vaccinations1"], model_dict, 61)
        model_dict = add_sim_to_dict_degreebased(result_df, "vaccination2_risk23", event_to_id["vaccinations2"], model_dict, 61)
        model_dict = add_sim_to_dict_degreebased(result_df, "br_risk23", event_to_id["br_leaving"], model_dict, 61)

        #model_dict = add_sim_geometry(result_df, model_dict, bc, pos_off)

        model_dict = get_die_out_prob(result_df, model_dict)
    return(model_dict)


def plot_bar(data, pos, color, ax, width = 0.05):
    b = data.mean()
    y_errormin = [b - data.quantile(0.025)]
    y_errormax = [data.quantile(0.975) - b]
    c = [y_errormin, y_errormax]
    ax.bar(pos, b, width = width, color = color)
    ax.errorbar(pos, b, yerr=c, fmt='.', color="r",  capsize=5, capthick=1, ms = 0)
    print("Mean:", np.round(b,2), "95\% CI:", np.round(data.quantile(0.025),2), np.round(data.quantile(0.975),2))



def ci_by_bootstrap(df, n_iterations = 1000, alpha = 0.05, seed = 1):
    np.random.rand(seed)
    n_size = len(df)  # Size of each sample
    bootstrap_means = []
    for _ in range(n_iterations):
        sample = df.sample(n=n_size, replace=True)
        bootstrap_means.append(sample.mean())

    bootstrap_means = pd.DataFrame(bootstrap_means)

    lower_bound = bootstrap_means.quantile(alpha/2)
    upper_bound = bootstrap_means.quantile(1-alpha/2)

    return(lower_bound, upper_bound)

def jan_only_formatter(x, pos=None):
    d = mdates.num2date(x)
    return d.strftime('%m/%Y') if d.month == 7 else ''


def scatter_hist(x, y):
    fig, axs = plt.subplot_mosaic([['histx', '.'],
                               ['scatter', 'histy']],
                              figsize=(6, 6),
                              width_ratios=(4, 1), height_ratios=(1, 4),
                              layout='constrained')
    # no labels
    ax_histx = axs['histx']
    ax_histy = axs['histy']
    ax = axs['scatter']
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # the scatter plot:
    ax.scatter(x, y)

    # now determine nice limits by hand:
    #binwidth = 0.01
    #xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    #lim = (int(xymax/binwidth) + 1) * binwidth

    #bins = np.arange(-lim, lim + binwidth, binwidth)
    ax_histx.hist(x, bins=100)
    ax_histy.hist(y, bins=100, orientation='horizontal')
    return(axs)
