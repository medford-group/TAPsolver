
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
from fenics_adjoint import *
import sys
import pandas as pd
import numpy as np
import math as mp
#import dijitso
import time
import os
import scipy
import copy
import warnings

from structures import *
from file_io import *
from mechanism_construction import *
from initial_conditions import *

# Return reduced functional

warnings.simplefilter(action='ignore', category=FutureWarning)
set_log_level(30)
tol = 1E-14

def forward_problem(pulse_time, reactor_data: reactor, mechanism_data: mechanism, initial_conditions_data: initial_conditions):



	######################## DEFINING INITIAL CONSTANTS #############################

	time_steps = pulse_time*1000
	dk = Constant(pulse_time/time_steps)
	cat_location = 1 - reac_input['Catalyst Location']

	# Standard Constants
	kbt = 1.38064852e-23
	hb = 6.62607004e-34
	Rgas = 8.314
	reference_pulse_size = 1
	

	######################## DEFINING OUTPUT FOLDER INFORMATION #############################

	def generateFolder(path_name):
		try:  
			os.mkdir(path_name)
		except OSError:  
			pass
		else: 
			pass

	path = './'+reactor_data.output_name+'_folder/'
	generateFolder(path)
			
	path_3 = './'+reactor_data.output_name+'_folder/flux_data/'
	generateFolder(path_3)

	user_data = pd.read_csv(input_file,header=None)
	user_data.to_csv(path+input_file,header=None,index=False)
	original_input_structure = user_data.copy()

	######################## READING KINETIC PARAMETER INFORMATION ###########################

	# Declare and define the constants of interest
	mechanism_copy = copy.deepcopy(mechanism_data)

	if mechanism_copy.rate_array == None:
		mechanism_copy = mechanism_constructor(mechanism_copy)
	mechanism_processes = mechanism_copy.elementary_processes.keys()
	
	r_Ga_in = {}
	r_dG_in = {}
	r_const = {}
	r_Ao = {}
	r_Ea = {}
	r_links = {}

	for j in mechanism_processes:
		if mechanism_copy.elementary_processes[j].forward.Ga['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'g':
			r_Ga_in["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.Ga['value'])
		if mechanism_copy.elementary_processes[j].forward.dG['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'g':
			r_dG_in["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.dG['value'])
		if mechanism_copy.elementary_processes[j].forward.k['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'k':
			r_const["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.k['value'])
		if mechanism_copy.elementary_processes[j].forward.Ao['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'a':
			r_Ao["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.Ao['value'])
		if mechanism_copy.elementary_processes[j].forward.Ea['value'] != None and mechanism_copy.elementary_processes[j].forward.default == 'a':
			r_Ea["kf"+str(j)] = Constant(mechanism_copy.elementary_processes[j].forward.Ea['value'])

		if mechanism_copy.elementary_processes[j].backward.k['value'] != None and mechanism_copy.elementary_processes[j].backward.default == 'k':
			r_const["kb"+str(j)] = Constant(mechanism_copy.elementary_processes[j].backward.k['value'])
		if mechanism_copy.elementary_processes[j].backward.Ao['value'] != None and mechanism_copy.elementary_processes[j].backward.default == 'a':
			r_Ao["kb"+str(j)] = Constant(mechanism_copy.elementary_processes[j].backward.Ao['value'])
		if mechanism_copy.elementary_processes[j].backward.Ea['value'] != None and mechanism_copy.elementary_processes[j].backward.default == 'a':
			r_Ea["kb"+str(j)] = Constant(mechanism_copy.elementary_processes[j].backward.Ea['value'])

	#print(r_Ga_in)
	#print(r_dG_in)
	#print(r_const)
	#print(r_Ao)
	#print(r_Ea)
	#print(r_links)


	for jnum,j in enumerate(r_Ga_in):
		#print(((kbt*reac_input['Reactor Temperature']/hb)*exp(-Ga_in["Ga"+str(jnum)])))
		#print(((kbt*reac_input['Reactor Temperature']/hb)*exp(-(Ga_in["Ga"+str(jnum)] - dG_in["dG"+str(jnum)]))))

		r_Ga_in[j] = Constant(r_Ga_in[j])
	
	for j in r_dG_in:
		r_dG_in[j] = Constant(r_dG_in[j])
	
	for j in r_const:
		r_const[j] = Constant(r_const[j])
	
	for j in r_Ao:
		r_Ao[j] = Constant(r_Ao[j])
	
	for j in r_Ea:
		r_Ea[j] = Constant(r_Ea[j])
	for j in r_links:
		r_links[j] = Constant(r_links[j])


	# Initialize the grid system, time step size, pulse size
	r_param, dx_r, dx2_r, frac_length, cat_location = establishMesh(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['Mesh Size'])
	