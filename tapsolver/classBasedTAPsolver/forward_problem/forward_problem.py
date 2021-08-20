
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
from fenics_adjoint import *
import sys
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import os
import scipy

from structures import *
from file_io import *
from mechanism_construction import *
from initial_conditions import *

def tap_simulate(reactor_data: reactor, mechanism_data: mechanism, initial_conditions_data: initial_conditions):

	kbt = 1.38064852e-23
	hb = 6.62607004e-34
	Rgas = 8.314

	warnings.simplefilter(action='ignore', category=FutureWarning)
	tape2 = Tape()
	tape2.clear_tape()
	set_working_tape(tape2)

	######################## DEFINING OUTPUT FOLDER INFORMATION #############################

	path = './'+reactor_data.output_name+'_folder/'
	generateFolder(path)
			
	path_3 = './'+reactor_data.output_name+'_folder/flux_data/'
	generateFolder(path_3)

	######################## READING KINETIC PARAMETER INFORMATION ###########################

	# Declare and define the constants of interest
	r_Ga_in = Ga_in.copy()
	r_dG_in = dG_in.copy()
	r_const = constants_input.copy()
	r_Ao = Ao_in.copy()
	r_Ea = Ea_in.copy()
	r_links = reactor_kinetics_input['linked parameters'].copy()
	linkForward = reactor_kinetics_input['link forward'].copy()
	linkBackard = reactor_kinetics_input['link backward'].copy()
	r_fit = fitting_input.copy()

	for jnum,j in enumerate(Ga_in):
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

	doe_form_pulse = True

	if reac_input['Fit Parameters'].lower() == 'true' or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in r_fit:
			if j.find('Ga') > -1: # 'Ga' and 'dG'
				controls.append(Control(r_Ga_in[j]))
	
			elif j.find('dG') > -1: # 'Ga' and 'dG'			
				controls.append(Control(r_dG_in[j]))
				
			elif j.find('Ao') > -1:
				controls.append(Control(r_Ao[j]))
	
			elif j.find('Ea') > -1:
				controls.append(Control(r_Ea[j]))

			elif j.find('{') > -1:
				controls.append(Control(r_links[j[1:-1]]))
	
			else:
				controls.append(Control(r_const[j]))
			legend_2.append(j)



