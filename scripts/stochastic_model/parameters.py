# model 0 has no behaviour change (at all, not even in 2022)
model0_parameters = {
	'diag_factor': 1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': True
}


# model 1 is continuation of 2022 parameters
model1_parameters = {
	'diag_factor': 1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}


# model 2 has no imports after calibration phase
model2_parameters = {
	'diag_factor': 1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': False,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}


# model 3 has no additional behaviour return after calibration phase
model3_parameters = {
	'diag_factor': 1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# stop behaviour return
	'return_rate': 0.0,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}



# model 4 no reinfection!
model4_parameters = {
	'diag_factor': 1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 1,

	'bc_off': False
}



########################################################################################################
# these models have (at least) one additional free parameter

# model 5 has shorter duration of 2nd infection
model_dpsi_parameters = {
	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}

# model 6 has fewer reported cases after major outbreak
model_diag_parameters = {
	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}

# model 7 has shorter duration of 2nd infection and overall reduced diagnosis 
model7_parameters = {
	'diag_factor': 0.35,
	'dpsi': 1.5,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}

# model 8 has higher diag rate
model8_parameters = {
	'diag_factor':1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': True,

	# behaviour return
	'return_rate': 0.03,

	'long_mpox_prob': 0,

	've_I': 0,

	'bc_off': False
}

# model final size is there to calculate the maximum outbreak size on a network of that shape
# behaviour change is turned off, so are vaccinations
model_fs_parameters = {
	'diag_factor':1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': False,

	# behaviour return
	'return_rate': 0.03,

	've_I': 1, # no reinfection

	'long_mpox_prob': 0,

	'bc_off': True,

	'vacc_off': True
}


# model stable point helps understanding the stable/endemic point of the unrestricted (no vacc no bc) system
model_sp_parameters = {
	'diag_factor':1,
	'dpsi': 1,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2025-12-31',

	# use imports
	'imports': False,

	# behaviour return
	'return_rate': 0.03,

	've_I': 0, # no reinfection

	'long_mpox_prob': 0,

	'bc_off': True,

	'vacc_off': True
}


# model 2022 
model2022_parameters = {
	'return_rate': 0.03, # estimated from first mpox paper, hopefully not important
	'initial_behaviour_reduction': 0.3448645539169543,#0.22, # estimate from first mpox paper
	'I0': 2, # source: 2022 paper
	#'waning': 0.23, # source, nat comms predicting...
	'brf': 0.031275514799677356, # from our mpox 2022 paper
	'share0': 21150/59394,

	# fitted
	'scale_E': 1,# 0.46525995626931915, # from first mpox paper!
	'scale_I': 2/5,#0.375866851595007, # from first mpox paper!

	# default parameters
	've': 0,
	've_I': 0,
	've_factor': 1,
	'dpsi': 1,

	# immune parameters
	'delta_l': 0.0028,
	'delta_s': 0.23,
	'x01': 80,
	'x02': 673,
	'f1': 0.89,
	'f2': 0.94,
	'x0inf': 1500,
	#'x0inf': 2000,
	'finf': 0.94,
	'ic50': 23.980360677686917,

	# adaptive diagnosis rate
	'time_dependent_diagnosis': False,

	# end date
	'end_date': '2022-10-31',

	# use imports
	'imports': True,

	'long_mpox_prob': 0,

	'bc_off': False
}


parameter_sets = {
	"model0": model0_parameters,
	"model1": model1_parameters,
	"model2": model2_parameters,
	"model3": model3_parameters,
	"model4": model4_parameters,
	"model_dpsi": model_dpsi_parameters,
	"model_diag": model_diag_parameters,
	#"model6": model6_parameters,
	"model7": model7_parameters,
	"model8": model8_parameters,
	"model_finalsize": model_fs_parameters,
	"model_stable_point": model_sp_parameters,
	"model_stable_point_dpsi": model_sp_parameters,
	"model_stable_point_2ndinf": model_sp_parameters,
	"model_long": model_dpsi_parameters,
	"model_testi": model_dpsi_parameters,
	#"model_diag": model1_parameters,
	"model_refit": model1_parameters,
	"model_2ndinf": model1_parameters,
	"model_2nddiag": model1_parameters,
	"model_base": model1_parameters,
	"model2022": model2022_parameters,
	"model22alt": model2022_parameters,
	'fitting': model1_parameters
}



