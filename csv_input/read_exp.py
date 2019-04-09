import pandas as pd
import numpy as np
import math as mp
import sys


def exp_data_fitting_3(species_list,time,folder):
	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'
	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+species_list[k]+'.csv',header=None)
	curve_fitting = {}
	exp_data = user_data
	for k_new in species_list:
		
		time_step = []
		times = []
		values = []	
		#Peak point
		peak_loc = user_data[k_new].iloc[user_data[k_new][1].idxmax()]
		peak2 = user_data[k_new].loc[user_data[k_new][0] == peak_loc[0]].index

		time_step.append(peak2[0])
		times.append(round(peak_loc[0],6))
		values.append(peak_loc[1])

		#times.append(peak_loc[0])
		#values.append(peak_loc[1])
		data = {}
		data['time_step'] = time_step
		data['times'] = times
		data['values'] = values

		curve_fitting[k_new] = data
	#sys.exit()
	return curve_fitting

def exp_data_fitting_2(species_list,time,folder):
	user_data = {}
	species_list = species_list[:len(species_list)]
	print(species_list)
	time_steps = int(time)
	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+species_list[k]+'.csv',header=None)
	curve_fitting = {}
	exp_data = user_data
	for k_new in species_list:
		time_step = []
		times = []
		values = []
		for j in range(0,time_steps):
			time_step.append(j)
			times.append(round(user_data[k_new].iloc[j,0],6))
			values.append(user_data[k_new].iloc[j,1])

		data = {}
		data['time_step'] = time_step
		data['times'] = times
		data['values'] = values
		curve_fitting[k_new] = data
	return curve_fitting

def exp_data_fitting(species_list,time,folder):
	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'
	print(species_list)
	
	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+species_list[k]+'.csv',header=None)
	curve_fitting = {}
	exp_data = user_data
	for k_new in species_list:
		
		time_step = []
		times = []
		values = []	
		#First data point
		time_step.append(1)
		times.append(round(user_data[k_new].iloc[1,0],6))
		values.append(user_data[k_new].iloc[1,1])
	
		#Peak point
		peak_loc = user_data[k_new].iloc[user_data[k_new][1].idxmax()]
		peak2 = user_data[k_new].loc[user_data[k_new][0] == peak_loc[0]].index

		test3 = int(round((peak2[0]+1)/2,0))
		mid_loc = user_data[k_new].iloc[test3,:]
		time_step.append(test3)
		times.append(round(mid_loc[0],6))
		values.append(mid_loc[1])

		time_step.append(peak2[0])
		times.append(round(peak_loc[0],6))
		values.append(peak_loc[1])

		#%50 point
		#value_test = 0.85*peak_loc[1]
		#sort_fif = user_data[k_new].iloc[(user_data[k_new][1]-value_test).abs().argsort()[:]]

		#for k in range(0,sort_fif.index[0]):
		#	if sort_fif.iloc[k,0] > peak_loc[0]:
		#		fif_point = sort_fif.iloc[k]
		#		break
		#	else:
		#		pass

		#test2 = user_data[k_new].loc[user_data[k_new][0] == fif_point[0]].index
		#time_step.append(test2[0])
		#times.append(round(fif_point[0],6))
		#values.append(fif_point[1])
		#%20 point
		#value_test = 0.2*peak_loc[1]
		#sort_twen = user_data[k_new].iloc[(user_data[k_new][1]-value_test).abs().argsort()[:]]
		
		#for k in range(0,sort_twen.index[0]):
		#	if sort_twen.iloc[k,0] > peak_loc[0]:
		#		twen_point = sort_twen.iloc[k]
		#		break
		#	else:
		#		pass

		#test2 = user_data[k_new].loc[user_data[k_new][0] == twen_point[0]].index
		#time_step.append(test2[0])
		#times.append(round(twen_point[0],6))
		#values.append(twen_point[1])
		

		#times.append(peak_loc[0])
		#values.append(peak_loc[1])
		data = {}
		data['time_step'] = time_step
		data['times'] = times
		data['values'] = values

		curve_fitting[k_new] = data
	#sys.exit()
	return curve_fitting

#input_list = ['CO','Inert']

#print(type(output_fitting['Inert']['time_step']))