from tap_sim import tap_simulation_function
from math import *
import os
import math
import numpy as np
import time
import pandas as pd
import sys
import matplotlib.pyplot as plt
from python_RRM import R_RRM_func, petal_plots_RRM
from experimental_data_odes import MKM_graphs, petal_plots_exp
from RRM_MKM import jacobian_visual

def base_run():
	
	user_data = pd.read_csv('./input_file.csv',header=None)
	
	rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
	rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
	rows_4, cols_4 = np.where(user_data == 'Reaction_Information')
	
	reactor_info = user_data.iloc[(2+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = user_data.iloc[2+rows_2[0]:rows_3[0]-1,:]
	data_storage = user_data.iloc[2+rows_3[0]:rows_4[0]-1,:]
	reaction_info = user_data.iloc[2+rows_4[0]:,:]
	
	reactor_kinetics_input = {}
	kinetic_parameters = {}
	
	for k in range(0,len(reactor_info.index)):
		try:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]

	for k in range(0,len(data_storage.index)):
		try:
			reactor_kinetics_input[data_storage.iloc[k,0]] = float(data_storage.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[data_storage.iloc[k,0]] = data_storage.iloc[k,1]

	for k in range(0,len(feed_surf_info.index)):
		try:
			reactor_kinetics_input[feed_surf_info.iloc[k,0]] = float(feed_surf_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[feed_surf_info.iloc[k,0]] = feed_surf_info.iloc[k,1]

	reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['length_reac']/2 -  0.5*(reactor_kinetics_input['factor'])*reactor_kinetics_input['length_reac']
	reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['factor'])*reactor_kinetics_input['length_reac'] 
	reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['length_reac']/2 -  0.5*(reactor_kinetics_input['factor'])*reactor_kinetics_input['length_reac']

	reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()

	for j in range(0,len(reaction_info.index)):
		kinetic_parameters['Ke'+str(j)] = float(reaction_info.iloc[j,1])
		if str(reaction_info.iloc[j,2]) != 'nan':
			kinetic_parameters['Kd'+str(j)] = float(reaction_info.iloc[j,2])
		else:
			pass
	kin_in = kinetic_parameters.copy()
	time_of_sim, graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)

	if reactor_kinetics_input['mkm_analysis'].lower() == 'true':
		MKM_graphs(kin_in,reactor_kinetics_input['reactions_test'],reactor_kinetics_input['output_file_name']+'_folder/thin_data',reactor_kinetics_input['Display_figure'].lower())
		if reactor_kinetics_input['petal_plots'].lower() == 'true':
			print('true_mkm')
	if reactor_kinetics_input['RRM_analysis'].lower() == 'true':
		R_RRM_func(legend_label[0:len(legend_label)-1],os.getcwd(),reactor_kinetics_input['output_file_name']+'_folder')
		if reactor_kinetics_input['petal_plots'].lower() == 'true':
			print('true_RRM')
	if reactor_kinetics_input['mkm_analysis'].lower() == 'true' and reactor_kinetics_input['RRM_analysis'].lower() == 'true':
		jacobian_visual(kin_in,reactor_kinetics_input['reactions_test'],reactor_kinetics_input['output_file_name']+'_folder/thin_data',reactor_kinetics_input['Display_figure'].lower(),legend_label[0:len(legend_label)-1],os.getcwd(),reactor_kinetics_input['output_file_name']+'_folder',legend_label[:-1],reactor_kinetics_input['number_of_pulses'])

	#sys.exit()

base_run()	