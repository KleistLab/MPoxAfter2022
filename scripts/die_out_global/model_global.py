import numpy as np
#import pandas as pd
import sys
#from util import *

class Global_die_out_model:

	def __init__(self, number_nodes, die_out_prob, reporting_countries_0, c):
		self.n = number_nodes
		self.rc_t = np.zeros_like(die_out_prob)
		self.rc_t[0] = reporting_countries_0
		self.dop_t = die_out_prob
		self.c = c[0] 		# coefficient of reporting countries 
		self.c_t = c[1] 	# coefficient of time 
		self.intercept = c[2]

		self.result = np.zeros((self.n, len(self.rc_t)))
		self.result[:reporting_countries_0,0] = 1

		self.weights = np.exp(-0.15 * np.arange(self.n)) # np.ones(self.n)
		self.grace_period = np.zeros(self.n)	# grace period turned off right now


	def run_model(self):
		# iterate through every day
		for i in range(1, len(self.dop_t)):
			# active cells: calculate die out probability
			filt = (self.result[:,i-1]==1) & (self.grace_period == 0)
			n_active = np.sum(filt)
			u = np.random.rand(n_active)
			tmp = (u >= self.dop_t[i]).astype(int) # switch condition, such that 1 is survival and 0 is die out
			self.result[filt, i] = tmp

			# inactive cells: calculate reintroduction probability
			n_passive = np.sum(self.result[:,i-1]==0)
			u = np.random.rand(n_passive)

			lam_import = self.weights[self.result[:,i-1]==0] * np.exp(self.intercept + self.c * self.rc_t[i-1] + self.c_t * i) 	# poisson regression log(lam) = c1*rc_t + c2*t + c3

			tmp = (u < 1 - np.exp(-np.maximum(lam_import,0))).astype(int)
			self.result[self.result[:,i-1]==0, i] = tmp

			# update reporting countries
			self.rc_t[i] = np.sum(self.result[:, i])
			#self.grace_period = np.zeros(self.n)
			#self.grace_period[self.result[:,i-1]==0] = tmp





