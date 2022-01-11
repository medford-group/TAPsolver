# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
from fenics_adjoint import *
import sys
import pandas as pd
import numpy as np
import math as mp
import time
import os
import scipy
import copy
import warnings
import matplotlib.pyplot as plt

#from structures import *
from .define_adspecies import define_adspecies
from .define_gas import define_gas
from .display_gasses import display_gasses
from .display_surface import display_surface
from .experimental_data import experimental_data
from .mechanism import mechanism
from .reactor import reactor
from .reactor_species import reactor_species
#from .read_old_input import read_old_input
from .TAPobject import TAPobject

#from file_io import *
from .new_experiments import new_experiments
from .read_CSV_mechanism import read_CSV_mechanism
from .read_CSV_reactor import read_CSV_reactor
from .read_CSV_reactor_species import read_CSV_reactor_species
from .read_experimental_data_object import read_experimental_data_object
from .read_mechanism_object import read_mechanism_object
from .read_reactor_object import read_reactor_object
from .read_reactor_species_object import read_reactor_species_object 
from .read_TAPobject import read_TAPobject 
from .read_transient_sensitivity import read_transient_sensitivity 
from .save_object import save_object
#from .vary_input_file import vary_input_file

#from .mechanism_construction import *
#from .construct_batch_equation import make_batch_equation
from .construct_f_equation import construct_f_equation
from .construct_f_equation_multiple_experiments import construct_f_equation_multiple_experiments
from .construct_rate_equations import rateEqs
from .display_elementary_processes import display_elementary_processes
from .elementary_process import elementary_process
from .elementary_process_details import elementary_process_details
from .mechanism_constructor import mechanism_constructor
from .mechanism_reactants import mechanism_reactants

#from reactor_species import *
#from reference_parameters import *
from .reference_parameters import load_standard_parameters

#from simulation_notes import *
from .timing_details import *
from .error_details import *
from .generate_folders import *

#from inverse_problem import *
from .define_fitting_species import curveFitting
#from .point_objective import point_objective
from .std_objective import stdEstablishment
from .total_objective import curveFitting

import numpy as np
import pandas as pd
import sys
from .forward_problem import forward_problem
import copy
import re

def transient_sensitivity(pulse_time, pulse_number, TAPobject_data_original: TAPobject):
	# method = 'tl' or 'fd'
	transient_data = {}
	TAPobject_data = copy.deepcopy(TAPobject_data_original)
	qoi = copy.deepcopy(TAPobject_data.parameters_of_interest)
	if TAPobject_data.tangent_linear_sensitivity == True:	
		transient_data = {}
		for k in TAPobject_data.reactor_species.gasses:
			transient_data[k] = {}
			for j in range(0,pulse_number):
				transient_data[k][j] = {}
				for n in qoi:
					transient_data[k][j][n] = []
		for j in qoi:
			TAPobject_data.parameters_of_interest = [j]
			new_data = forward_problem(pulse_time,pulse_number,TAPobject_data)
			for k in TAPobject_data.reactor_species.gasses:
				transient_data[k][0][n] = new_data[k]
		print(transient_data['CO'][0]['TAPobject_data.mechanism.elementary_processes[0].forward.k'])
		save_object(transient_data,'./'+TAPobject_data.output_name+'/TAP_test.json')
		deriviative_information = read_transient_sensitivity('./'+TAPobject_data.output_name+'/TAP_test.json')	
		
	if TAPobject_data.finite_difference_trans_sensitivty == True:
		original_name = copy.deepcopy(TAPobject_data.output_name)
		forward_problem(pulse_time,pulse_number,TAPobject_data)
		old_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')
		TAPobject_data.output_name = 'df_flux'
		for k in TAPobject_data.reactor_species.gasses:
			transient_data[k] = {}
			for j in range(0,pulse_number):
				transient_data[k][j] = {}
				for n in qoi:
					transient_data[k][j][n] = []
		for j in qoi:

			old_value = copy.deepcopy(eval(j))
			
			if 'delay' not in j:
				new_value = old_value+0.0001*old_value
			else:
				new_value = old_value+0.001

			if 'elementary_processes' in j:
				elementary_step = float(re.split(r"[\[\]]",j)[1])
				if 'forward' in j:
					if 'k' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.k = new_value
					elif 'Ao' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ao = new_value
					elif 'Ea' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ea = new_value
					elif 'Ga' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ga = new_value
					elif 'dG' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.k = new_value
				
				if 'backward' in j:
					if 'k' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.k = new_value
					elif 'Ao' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ao = new_value
					elif 'Ea' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ea = new_value
					elif 'Ga' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ga = new_value
					elif 'dG' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.k = new_value
			elif 'temperature' in j:
				TAPobject_data.reactor_species.temperature = new_value
			elif 'delay' in j:
				gas_name = re.split(r'[\"\"]',re.split(r"[\[\]]",j)[1])[1]
				#print(TAPobject_data.reactor_species.gasses['CO'].delay)
				TAPobject_data.reactor_species.gasses[str(gas_name)].delay = new_value
			elif 'intensity' in j:
				gas_name = re.split(r'[\"\"]',re.split(r"[\[\]]",j)[1])[1]
				TAPobject_data.reactor_species.gasses[gas_name].intensity = new_value

			change_parameter = new_value

			step_size = new_value - old_value
			forward_problem(pulse_time,pulse_number,TAPobject_data)

			new_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')

			for x in TAPobject_data.reactor_species.gasses:
				for z in new_data[x]:
					temp_change = [(a_i - b_i)/change_parameter for a_i, b_i in zip(new_data[x][z], old_data[x][z])]
					transient_data[x][int(z)][j] = temp_change

			if 'elementary_processes' in j:
				elementary_step = float(re.split(r"[\[\]]",j)[1])
				if 'forward' in j:
					if 'k' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.k = old_value
					elif 'Ao' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ao = old_value
					elif 'Ea' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ea = old_value
					elif 'Ga' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ga = old_value
					elif 'dG' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].forward.k = old_value
				
				if 'backward' in j:
					if 'k' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.k = old_value
					elif 'Ao' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ao = old_value
					elif 'Ea' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ea = old_value
					elif 'Ga' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ga = old_value
					elif 'dG' in j:
						TAPobject_data.mechanism.elementary_processes[elementary_step].backward.k = old_value
			elif 'temperature' in j:
				TAPobject_data.reactor_species.temperature = old_value
			elif 'delay' in j:
				gas_name = re.split(r'[\"\"]',re.split(r"[\[\]]",j)[1])[1]
				TAPobject_data.reactor_species.gasses[gas_name].delay = old_value
			elif 'intensity' in j:
				gas_name = re.split(r'[\"\"]',re.split(r"[\[\]]",j)[1])[1]
				TAPobject_data.reactor_species.gasses[gas_name].intensity = old_value
		print(transient_data['CO'][0]['TAPobject_data.mechanism.elementary_processes[0].forward.k'])
		print(transient_data['O2'][0]['TAPobject_data.mechanism.elementary_processes[0].forward.k'])
		print(transient_data['O2'][0]['TAPobject_data.mechanism.elementary_processes[0].forward.k'])
		
		TAPobject_data.output_name = original_name
		save_object(transient_data,'./'+original_name+'/TAP_test.json')
		deriviative_information = read_transient_sensitivity('./'+original_name+'/TAP_test.json')

	if True == False:
	
		q_matricies = {}
	
		for k in TAPobject_data.reactor_species.gasses:
			parameter_options = list(deriviative_information[k][0].keys())
			temp_sensitivity_storage = {}
			for j in parameter_options:
				temp_sensitivity_storage[j] = transient_data[k][0][j]
			for j in range(1,pulse_number):
				for z in parameter_options:
					temp_sensitivity_storage[z].extend(transient_data[k][j][z])
			parameter_options_2 = list(temp_sensitivity_storage.keys())
			for jnum,j in enumerate(parameter_options_2):
				if jnum == 0:
					q_matricies[k] = np.array([temp_sensitivity_storage[j]])
				else:
					q_matricies[k] = np.append(q_matricies[k], np.array([temp_sensitivity_storage[j]]), 0)
		
		
			#for j in deriviative_information[k][0]:
	
		for jnum,j in enumerate(list(q_matricies.keys())):

			for znum,z in enumerate(list(q_matricies.keys())):
				if jnum == 0 and znum == 0:
					final_q_array = q_matricies[j].dot(np.transpose(q_matricies[z]))
				else:
					final_q_array += q_matricies[j].dot(np.transpose(q_matricies[z]))
 
		pd.DataFrame(final_q_array).to_csv('./'+original_name+'/covariance_matrix.csv')
		with open('./'+original_name+'/optimal_criteria.txt', 'w') as f:
			f.write('Determinant')
			f.write('\n')
			f.write(str(np.linalg.det(np.linalg.inv(final_q_array))))
			f.write('\n')
			f.write('Trace')
			f.write('\n')
			f.write(str(np.trace(np.linalg.inv(final_q_array))))
		f.close()