from fenics import *
from fenics_adjoint import *
from func_sim import *
from vari_form import *
from reac_odes import *
import mpmath
import matplotlib.pyplot as plt
import matplotlib
#matplotlib.use('Agg')
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import imageio
import csv
import ast
import shutil
import sys
import os
import scipy
import pip
import pkg_resources
import ufl
import re
from ufl import sqrt,exp,ln
from shutil import copyfile
#from hippylib import *
import warnings
#from ipyopt import *
#from ipopt import *
#from pyadjoint import ipopt

#from gryphon import *

def general_run(timeFunc,uncertainty_quantificaiton=None,optimization=None,fitting_gif=None,input_file = './input_file.csv',sensitivityType=None,noise=False):
	
	fit_temperature = False

	warnings.simplefilter(action='ignore', category=FutureWarning)
	
	#############################################################
	################ TAPsolver FEniCS Function ##################
	#############################################################
	
	def tap_simulation_function(reactor_kinetics_input,constants_input,Ao_in,Ea_in,Ga_in,dG_in,fitting_input,arrForward,arrBackward,gForward):
	
		# Define or clear working Tape of FEniCS / Dolfin-Adjoint
		tape2 = Tape()
		tape2.clear_tape()
		set_working_tape(tape2)

		#############################################################
		#### GENERATE OUTPUT FOLDERS (depending on input values) ####
		#############################################################
		
		kVals = constants_input.copy()
		reac_input = reactor_kinetics_input
	
		reac_input['Pulse Duration'] = timeFunc
		reac_input['Time Steps'] = reac_input['Pulse Duration']*1000
	
		sens_type = sensitivityType
		if (sens_type != 'trans' and sens_type != 'total') and reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true':
			print('Error: sensitivity analysis must be "trans" or "total"')
			sys.exit()

		#reac_input['Optimization Method'] = 'nonlinear'
		reac_input['Optimization Method'] = 'L-BFGS-B' #
		reac_input['Objective Points'] = 'all'
		reac_input['Store Outlet Flux'] = 'TRUE'
		reac_input['Store Graph'] = 'FALSE'
		reac_input['Display Graph'] = 'FALSE'
		reac_input['Scale Output'] = 'FALSE'
		reactor_kinetics_input['Reactor Type'] = 'tap'

		if 'thermo equations' in reac_input.keys():
			if len(reac_input['thermo equations']) > 0:
				thermoConstraints = True

			else:
				thermoConstraints = False

		else:
			thermoConstraints = False

		if thermoConstraints == True:
			thermoReactions = {}
			thermoStoich = {}
			for znum,z in enumerate(reac_input['thermo equations']):
				thermoReactions[znum] = []
				thermoStoich[znum] = []
				tempValue = z.split('+')
				for j in tempValue:
					if j.find('*') != -1:
						thermoSplit = j.split('*')
						thermoSplit[1] = thermoSplit[1].replace('r','')
						thermoReactions[znum].append(int(thermoSplit[1]))
						thermoSplit[0] = thermoSplit[0].replace('(','')
						thermoSplit[0] = thermoSplit[0].replace(')','')
						thermoStoich[znum].append(float(thermoSplit[0]))
					else:
						j = j.replace('r','')
						thermoReactions[znum].append(int(j))
						thermoStoich[znum].append(1)
		
		path = './'+reac_input['Output Folder Name']+'_folder/'
		generateFolder(path)
			
		path_4 = './'+reac_input['Output Folder Name']+'_folder/surface_data/'
		generateFolder(path_4)

		path_3 = './'+reac_input['Output Folder Name']+'_folder/species_data/'
		generateFolder(path_3)
	
		if reac_input['Fit Parameters'].lower() == 'true' or fit_temperature == True:
			path_6 = './'+reac_input['Output Folder Name']+'_folder/fitting/'
			generateFolder(path_6)
	
		runge_kutta_approach = False

		# Declare and define the constants of interest
		r_const = constants_input
		r_Ao = Ao_in
		r_Ea = Ea_in
		r_fit = fitting_input
	
		for j in Ga_in:
			Ga_in[j] = Constant(Ga_in[j])
	
		for j in dG_in:
			dG_in[j] = Constant(dG_in[j])
	
		for j in r_const:
			r_const[j] = Constant(r_const[j])
	
		for j in r_Ao:
			r_Ao[j] = Constant(r_Ao[j])
	
	
		for j in r_Ea:
			r_Ea[j] = Constant(r_Ea[j])

		doe_form_pulse = True
		doe_form_surf = True
	
		if reac_input['Fit Parameters'].lower() == 'true' or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true'):
			
			controls = []
			legend_2 = []
			for j in r_fit:
				if j.find('Ga') > -1: # 'Ga' and 'dG'
					controls.append(Control(Ga_in[j]))
	
				elif j.find('dG') > -1: # 'Ga' and 'dG'			
					controls.append(Control(dG_in[j]))
				
				elif j.find('Ao') > -1:
					controls.append(Control(r_Ao[j]))
	
				elif j.find('Ea') > -1:
					controls.append(Control(r_Ea[j]))
				else:
					controls.append(Control(r_const[j]))
				legend_2.append(j)
			print('Fitting the parameters: ')
			print(controls)

		# Store input file in output folder	
		user_data = pd.read_csv(input_file,header=None)
		user_data.to_csv(path+input_file,header=None,index=False)
		original_input_structure = user_data.copy()
		
		#############################################################
		#### READ AND DEFINE INITIAL PARAMETERS FROM INPUT FILE #####
		#############################################################
		reac_input['Reference Pulse Size'] = 1
		set_log_level(30)
		tol = 1E-14
	
		kbt = 1.38064852e-23
		hb = 6.62607004e-34
		Rgas = 8.314
	
		# Read the initial composition of the catalyst
		if ',' in str(reac_input['Initial Surface Composition']):
			rangesurface_species = list(reversed(reac_input['Initial Surface Composition'].split(',')))
			reac_input['Initial Surface Composition'] = list(map(float, rangesurface_species))
		else:
			reac_input['Initial Surface Composition'] = [reac_input['Initial Surface Composition']]
			rangesurface_species = reac_input['Initial Surface Composition']
		
		# Initialize the grid system, time step size, pulse size
		dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
		
		constantTemp = reac_input['Reactor Temperature']
		constantTemp = Constant(constantTemp)
		testTemp = reac_input['Reactor Temperature']
		testTemp = Constant(testTemp)
	
		#if fit_temperature == True:
		#	controls = []
		#	new = Control(constantTemp)
		#	controls.append(new)
	
		# Construct the rate expression based on the microkinetic model
		necessary_values, rate_array, rev_irr = make_batch_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],arrForward,arrBackward,gForward,fit_temperature)
		
		if thermoConstraints == True:
			
			thermoReactants = {}
			thermoProducts = {}
			thermoReactantsValue = {}
			thermoProductsValue = {}

			#thermoReactants = []
			#thermoProducts = []
			#thermoReactantsValue = []
			#thermoProductsValue = []

			for jnum,jval in enumerate(thermoReactions):
				thermoReactants[jnum] = []
				thermoProducts[jnum] = []
				thermoReactantsValue[jnum] = []
				thermoProductsValue[jnum] = []

				for znum,zval in enumerate(thermoReactions[jval]):

					if rev_irr[jval-1] != 1:
						print('')
						print('To set a thermodynamic constraints, all reactions must be reversible.')
						print('Currently, elementary reaction '+str(jval)+' is irreversible and included in the thermo. equation.')
						sys.exit()

					else:
						for zen in range(0,necessary_values['molecules_in_gas_phase']):
							print(rate_array[zval-1])
							if rate_array[zval-1][zen] == 1:

								thermoProductsValue[jnum].append(necessary_values['reactants'][zen])						
								
								if thermoStoich[jnum][znum] != 1:
									thermoProducts[jnum].append('('+str(thermoStoich[jnum][znum])+')*'+necessary_values['reactants'][zen])
									
								else:
									thermoProducts[jnum].append(necessary_values['reactants'][zen])	

							elif rate_array[zval-1][zen] == -1:

								thermoReactantsValue[jnum].append(necessary_values['reactants'][zen])
								if thermoStoich[jnum][znum] != 1:
									thermoReactants[jnum].append('('+str(thermoStoich[jnum][znum])+')*'+necessary_values['reactants'][zen])
								else:
									thermoReactants[jnum].append(necessary_values['reactants'][zen])		
				print(thermoReactantsValue)
				print(thermoProductsValue)
				print(thermoStoich)

				overallReactionString = ''
			
				for znum in range(0,len(thermoReactants[jnum])):
					overallReactionString = overallReactionString+thermoReactants[jnum][znum]
					if znum != len(thermoReactants[jnum])-1:
						overallReactionString = overallReactionString + ' + '

				overallReactionString = overallReactionString + ' <-> '
				
				for znum in range(0,len(thermoProducts[jnum])):	
					overallReactionString = overallReactionString + thermoProducts[jnum][znum]
					if znum != len(thermoProducts[jnum])-1:
						overallReactionString = overallReactionString + ' + '
				print('Thermodynamic Constraint Equation:')
				print(overallReactionString)
				print('')

				if reactor_kinetics_input['thermo values'][jnum] == 'none':
					freeEnergyValue = 0

					try:
						for znum in range(0,len(thermoReactantsValue[jnum])):
							calcValue = (-1)*(thermoStoich[jnum][znum])*molecularProperties(thermoReactantsValue[jnum][znum],'freeEnergy',temperature=reac_input['Reactor Temperature'])
							freeEnergyValue += calcValue

						for znum in range(0,len(thermoProductsValue[jnum])):

							calcValue = (thermoStoich[jnum][znum])*molecularProperties(thermoProductsValue[jnum][znum],'freeEnergy',temperature=reac_input['Reactor Temperature'])
							
							freeEnergyValue += calcValue
					except:
						print('Failed to calculate free energy of reaction automatically. Please enter free energy manually.')
						sys.exit()

					reactor_kinetics_input['thermo values'][jnum] = freeEnergyValue
					print('Free energy calculated automatically: ')
					print(reactor_kinetics_input['thermo values'][jnum])
					print('')
					

				else:
					print('Free energy entered manually: ')
					print(reactor_kinetics_input['thermo values'][jnum])
					print('')

		if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or (reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' and sens_type == 'total') or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='1,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='1'

		else:
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='0,'
			
		reactor_kinetics_input['Objective Species'] = speciesString


		#############################################################
		######## INITIALIZATION OF FINITE ELEMENTS AND MESH #########
		#############################################################
	
		# Define the base mesh (the initial uniform mesh)
		reac_input['Mesh Size'] = 3
		mesh = UnitIntervalMesh(reac_input['Mesh Size'])
	
		# Generate element space (for all observables)
		P1 = FiniteElement('CG',mesh.ufl_cell(),1)
	
		if reac_input['reactions_test'] != ['INERT_ONLY']:
			test_new = eval(necessary_values['element'])
			element = MixedElement(test_new)
			V = FunctionSpace(mesh,element)
			V_du = FunctionSpace(mesh,P1)
		else:
			V = FunctionSpace(mesh,P1)
			V_du = FunctionSpace(mesh,P1)
	
		all_molecules = necessary_values['gas_num']
	
		u = Function(V)
		if runge_kutta_approach == True:
			W = Function(V)
			W.interpolate(Constant(0.0))
		u_n = Function(V)
		u_temp = Function(V)
		
		graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data = initializeVariableDictionaries(necessary_values,all_molecules,V,u,u_n)
	
		meshCells = reac_input['Mesh Size']
		
		monitored_gas = necessary_values['molecules_in_gas_phase']
		
		fig2,ax2,legend_label,header,colors = establishOutletGraph('batch',necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],None,reac_input['Scale Output'])
		ax2.set_xlim(0,reac_input['Pulse Duration'])
	
		if str(reac_input['Objective Species']).find(',') != -1:
			objSpecies = list(reac_input['Objective Species'].split(','))
		else:
			objSpecies = [str(int(reac_input['Objective Species']))]
	
		t = Constant(0)
	
		#############################################################
		######## EVALUATE AND INITIALIZE PDE FOR FEniCS #############
		#############################################################
	
		try:
			theta = 1
			Ftemp = eval(necessary_values['F'])
			theta = 0.5
			F = eval(necessary_values['F'])
		except NameError:
			errorOutput(reac_input['reactions_test'])
		
		constantT = Constant(0)
		if doe_form_pulse == True:

			reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))

			intensityFunctions ={}
			intensConst = {}
			for jk in range(0,len(reactant_feed)):
				intensConst['inten_'+str(jk)] = Constant(reactant_feed[jk])
				
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
			
			timeConst = {}

			for jk in range(0,len(reactant_time)):
				timeConst['sST_'+str(jk)] = Constant(reactant_time[jk]+0.001)

			fnew = eval(batch_functions(reac_input['reactions_test']))

			F += fnew

		if doe_form_surf == True:
			
			inComp = {}
			print(reac_input['Initial Surface Composition'])
			reac_input['Initial Surface Composition'].reverse()

			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):
				inComp[z_num] = Constant(z_sur)
			
			try:
				F += eval(batch_surface_functions(reac_input['reactions_test']))
			except:
				pass
		if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True or reac_input['Fitting Gif'].lower() == True:
			
			if reac_input['Uncertainty Quantification'].lower() == 'true':
				print("Uncertainty Quantification")
			try:
				if type(reac_input['Objective Points']) == float:
					output_fitting = pointFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
				elif reac_input['Objective Points'] == 'all':
					output_fitting = curveFitting(legend_label[:int(len(legend_label))],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
					
				else:
					print('Objective Points defined incorrectly')
					sys.exit()
			except TypeError:
				print('Objective Point Input Is Not Valid')
				sys.exit()
	
		sens_time_list = []
		
		if reac_input['Sensitivity Analysis'].lower() == 'true':
			if sens_type == 'trans':
				path_2 = reac_input['Output Folder Name']+'_folder/sensitivity/'
				generateFolder(path_2)
			
				path_molecules = path_2+reac_input['Sensitivity Parameter']
				generateFolder(path_molecules)
				sensFolder = path_molecules
	
			elif sens_type == 'total':
				path_2 = reac_input['Output Folder Name']+'_folder/sensitivity/'
				generateFolder(path_2)			
	
			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit() 
	
		if reac_input['Uncertainty Quantification'].lower() == 'true':
			path_2 = reac_input['Output Folder Name']+'_folder/UQ/'
			generateFolder(path_2)
			hessFolder = path_2
			
		#############################################################
		######## CONVERTING RAW CONCENTRATION TO OUTLET FLUX ########
		#############################################################
	
		storeDataFunc(reac_input['Store Outlet Flux'],reac_input['Output Folder Name'])
		
		if sens_type == 'trans':
			storeSens(reac_input['Sensitivity Analysis'],reac_input['Output Folder Name'],monitored_gas,legend_label)
		elif sens_type == 'total':
			pass
		else:
			pass
			#print('Sensitivity analysis is not properly defined.')
			#sys.exit()
	
	
		#############################################################
		############# DEFINE VARIATIONAL PROBLEM ####################
		#############################################################
	
		J = derivative(F,u)
		Jtemp = derivative(Ftemp,u)
	
		# Define a constrained variational problem solver
		reac_input['Variational Solver'] = 'newton'
	
		if reac_input['Variational Solver'].lower() == 'constrained':
			snes_solver_parameters = {"nonlinear_solver": "snes","snes_solver": {"linear_solver": "lu","line_search":'basic',"maximum_iterations": 10,"report": False,"error_on_nonconvergence": False}}
			
			lower = Function(V)
			upper = Function(V) 
			
			ninfty = Function(V); ninfty.vector()[:] = 0
			pinfty = Function(V); pinfty.vector()[:] =  np.infty
	
			
			problem = NonlinearVariationalProblem(F,u,None,J)
			
			problem.set_bounds(ninfty,pinfty)
	
			solver = NonlinearVariationalSolver(problem)
			solver.parameters.update(snes_solver_parameters)
		
		# Define a newton variational problem solver
		elif reac_input['Variational Solver'].lower() == 'newton':
			problem = NonlinearVariationalProblem(F,u,None,J)
			solver = NonlinearVariationalSolver(problem)
	
			problemtemp = NonlinearVariationalProblem(Ftemp,u,None,Jtemp)
			solvertemp = NonlinearVariationalSolver(problemtemp)
	
			solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
			solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10
	
		#############################################################
		############# PARSE INLET PULSE COMPOSITION #################
		#############################################################
		
		if '/' in str(reac_input['Pulse Size']):
			list_of_feed = reac_input['Pulse Size'].split('/')
			list_of_time = reac_input['Pulse Time'].split('/')
			pulse_variation = len(list_of_feed)
			pulse_time = len(list_of_time)
			reactant_feed = []
			reactant_time = []
			for k in range(0,len(reac_input['Pulse Size'].split('/'))):
				reactant_feed.append(list_of_feed[k].split(','))
				reactant_time.append(list_of_time[k].split(','))
		else:
			pulse_variation = 1
			reactant_feed = []
			reactant_time = []
			
			try:
				#reactant_feed.append(reac_input['Pulse Size'].split(','))
				reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))
				reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
			except AttributeError:
				#reactant_feed.append(reac_input['Pulse Size'])
				reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))
				reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
		
		#############################################################
		####### STORE DATA BASED ON USER SPECIFICATION ##############
		#############################################################
	
		current_reac_set = 1
		tot_objective = 0
		reac_input['Thin-Zone Analysis'] = 'true'
		if reac_input['Thin-Zone Analysis'].lower() == 'true':
			rateStrings = rateEqs(rate_array,rev_irr,gForward,arrForward,arrBackward)
	
		if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True:
			for k_fitting in range(0,len(legend_label[:int(len(legend_label))])):
				if objSpecies[k_fitting] == '1':
					for timeStep in range(0,len(output_fitting[legend_label[k_fitting]]['times'])):
						output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)
	
		
		if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':
	
			if reac_input['Uncertainty Quantification'].lower() == 'true' or (reac_input['Sensitivity Analysis'].lower() == 'true' and sens_type == 'trans'):
				
				if reac_input['Sensitivity Parameter'].find('Ga') > -1:
					c = Ga_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = Ga_in[reac_input['Sensitivity Parameter']]
				elif reac_input['Sensitivity Parameter'].find('dG') > -1:
					c = dG_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = dG_in[reac_input['Sensitivity Parameter']]
				elif reac_input['Sensitivity Parameter'].find('Ao') > -1:
					c = Ao_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = Ao_in[reac_input['Sensitivity Parameter']]
				elif reac_input['Sensitivity Parameter'].find('Ea') > -1:
					c = Ea_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = Ea_in[reac_input['Sensitivity Parameter']]
				else:
					c = r_const[reac_input['Sensitivity Parameter']]
					c.tlm_value = r_const[reac_input['Sensitivity Parameter']]
	
	
				#c = r_const[]
				#c.tlm_value = r_const[reac_input['Sensitivity Parameter']]
				#c2 = r_const['kf1']
				SV_du = FunctionSpace(mesh,P1)
				Sw_new = Expression('A',A=Constant(1),degree=0)
				Sw_new2 = interpolate(Sw_new,SV_du)
				Sw3 = project(Sw_new2,SV_du)
			
			elif sens_type == 'total':
				pass
	
			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit()
	
	
		if reac_input['Sensitivity Analysis'].lower() == 'true':
			
			if sens_type == 'trans':
				sensFuncs = {}
	
				for k_gasses in range(0,len(necessary_values['reactants'])):

					sensFuncs[str(k_gasses)] = []		

			elif sens_type == 'total':
				pass
			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit()
	
		sensAll = []
		sensAll2 = []
	
		if True == False:
			pass
	
		#############################################################
		############ RUN SIMULATION FOR EACH PULSE ##################
		#############################################################
	
		for k_pulse in range(0,1):
	
			start_time = time.time()
			
			species_pulse_list = reactant_feed#reactant_feed[current_reac_set-1]
			
			if type(reactant_time[current_reac_set-1]) == list:
				species_time = reactant_time[current_reac_set-1]
			else:
				species_time = reactant_time
	 
			dt = reac_input['Pulse Duration']/reac_input['Time Steps']
			
			for kTimeStep,kTime in enumerate(species_time.copy()):
				tNew = dt*round(float(kTime)/dt)
				species_time[kTimeStep] = round(tNew,6)		
			
			
			if sens_type == 'trans':
				sensitivity_output = {}
				RRM_der = {}
			
				for k_sens in range(monitored_gas):
					sensitivity_output[k_sens] = []
	
				for k_sens in range(len(necessary_values['reactants'])):
					RRM_der[k_sens] = []
			
			elif sens_type == 'total':
				pass
	
			else:
				#print('Sensitivity analysis is not properly defined.')
				pass
				#sys.exit()
	
			graph_data = {}
			u_graph_data = {}
			
			for k_gasses in range(0,necessary_values['molecules_in_gas_phase']):
				graph_data['conVtime_'+str(k_gasses)] = []
				u_graph_data['conVtime_'+str(k_gasses)] = []
			
			surf_data['timing'] = []
			
			for j_gasses in range(necessary_values['molecules_in_gas_phase'],len(necessary_values['reactants'])-1):
				surf_data['conVtime_'+str(j_gasses)] = []
			
			graph_data['timing'] = []
	
			cat_data['timing'] = []
			
			for z_gasses in range(0,all_molecules+1):
				cat_data['conVtime_'+str(z_gasses)] = []
	
			cat_data['rawData'] = []
			
			t = 0

			test_new = project(-u.dx(0),V)
			test_new_split = test_new.split(deepcopy=False)
	
			w_new = Expression("1", degree=0)
			w_new2 = interpolate(w_new,V_du)
	
			x_dim = list(range(0, int(reac_input['Mesh Size'])+1))
			
			cum_sum = 0
	
			meshCells = int(reac_input['Mesh Size']) 
	
			#############################################################
			############ SOLVE PDEs AT EACH TIME STEP ###################
			#############################################################

			# Design of experiment form
	
			while t <= reac_input['Pulse Duration']:
				graph_data['timing'].append(t)
				for k in range(0,monitored_gas):
					new_val = (( u_n.vector().get_local()[k+all_molecules]))

					graph_data['conVtime_'+str(k)].append((new_val))
	
				mesh_size = reac_input['Mesh Size']

				if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True:
	
					# Define the objective function 	
					objectiveAnalysis = True
					if objectiveAnalysis == True:
						for k_fitting in range(0,len(legend_label[:int(len(legend_label))])):
						
							if objSpecies[k_fitting] == '1':
									
								if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
									
									c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
									
									w_new = Expression('1',degree=0)
									w_new2 = interpolate(w_new,V_du)
									w3 = project(w_new2,V_du)
	
									try:
										if legend_label[k_fitting] != 'Inert':
											jfunc_2 += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dx())							
										else:
											pass
	
									except UnboundLocalError:
										if legend_label[k_fitting] != 'Inert':
											w_temp_2 = Expression('1',degree=0) # deltaG = sum(-R*T*ln(kf/kb))
											w_temp2_2 = interpolate(w_temp_2,V_du)
											w4_2 = project(w_temp2_2,V_du)		
	
											jfunc_2 = assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dx())
											
											if thermoConstraints == True:
												w_temp = {}
												w_temp2 = {}
												w4 = {}
												tempFunc = {}

												for kip in thermoReactions.keys():
													w_temp[kip] = Expression(str(reactor_kinetics_input['thermo values'][kip]),degree=0) # deltaG = sum(-R*T*ln(kf/kb))
													w_temp2[kip] = interpolate(w_temp[kip],V_du)
													w4[kip] = project(w_temp2[kip],V_du)		
													thermoWeight = 1e-2	

													for jnum,jval in enumerate(thermoReactions[kip]):
																				
														if jnum == 0:
															tempFunc[kip] = thermoStoich[kip][jnum]*(-0.008314*reac_input['Reactor Temperature'])*ln(r_const["kf"+str(jval-1)]/r_const["kb"+str(jval-1)])
														else:
															tempFunc[kip] += thermoStoich[kip][jnum]*(-0.008314*reac_input['Reactor Temperature'])*ln(r_const["kf"+str(jval-1)]/r_const["kb"+str(jval-1)])

													#tempFunc = (-w4)*thermoWeight
													#tempFunc = (tempFunc)*thermoWeight
													#tempFunc_2 = (w4_2)*thermoWeight
													tempFunc[kip] = (tempFunc[kip]-w4[kip])*thermoWeight
													jfunc_2 += assemble(inner(tempFunc[kip],tempFunc[kip])*dx())
													#jfunc_2 += assemble(inner(tempFunc,tempFunc)*dx())
													print('Simulated Thermo Value')
													print(jfunc_2)
													#sys.exit()
										else:
											pass
				
				if runge_kutta_approach == False:	
					if round(t,6) not in species_time:
						#timeDiff = round(t,6)
						try:
							if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
								# C.N. method
								if t > 0.0011+timeDiff:
									solver.solve()
								# B.E. method
								else:
									solvertemp.solve()
		
									if round(t,6) == 0.001+timeDiff:
										dt = reac_input['Pulse Duration']/reac_input['Time Steps']
										dk.assign(dt)
										u_n.assign(u)
										solver.solve()
							else:
								if t > 0.0011+timeDiff:
									solver.solve(annotate = False)
								else:
									solvertemp.solve(annotate=False)
		
									if round(t,6) == 0.001+timeDiff:
										dt = reac_input['Pulse Duration']/reac_input['Time Steps']
										dk.assign(dt)
										u_n.assign(u)
										solver.solve()
								
							#############################################################
							################### STORE THIN-ZONE DATA ####################
							#############################################################
		
							if reac_input['Thin-Zone Analysis'].lower() == 'true':
								arrayNew = np.vstack((cat_data['rawData'],u_n.vector().get_local()))
								cat_data['rawData'] = arrayNew

	
						except RuntimeError:
							print('Time Step Failure')
							sys.exit()
						
					else:
						if dt == 0.0001:
							pass
						else:
							dt = 0.0001
							dk.assign(dt)
		
						if round(t,6) in species_time:
							timeDiff = round(t,6)
		
						if reac_input['Reactor Type'] == 'tap':
							
							#if reac_input['Fit Parameters'].lower() == 'true':
							#	if t == 0:
							#		for k in range(0,int(reac_input['Number of Reactants'])):
							#			u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+additionalCells))+k-(all_molecules)] = -1e-10			
							
							if reac_input['Thin-Zone Analysis'].lower() == 'true':
								
								if round(t,6) != 0:
									arrayNew = np.vstack((cat_data['rawData'],u_n.vector().get_local()))
									cat_data['rawData'] = arrayNew
								else:
									cat_data['rawData'] = u_n.vector().get_local()
		
						try:
							if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
								solvertemp.solve()
		
							else:
								solvertemp.solve(annotate = False)
		
						except RuntimeError:
							print('Time Step Failure')
							sys.exit()
				
				elif runge_kutta_approach == True:
					timeDomain = [0.0, 5.0]
					#problem = NonlinearVariationalProblem(F,u,bcs,J)
					obj = ESDIRK(timeDomain, W, F, bcs=bcs)

				#############################################################
				################ STORE SENSITIVITY DATA #####################
				#############################################################
	
				if reac_input['Sensitivity Analysis'].lower() == 'true':
					
					if sens_type == 'trans':
	
						u_final = u.split(deepcopy=False)
						should_it = int(round(t*reac_input['Time Steps']/reac_input['Pulse Duration'],0))
					
						for k in range(0,monitored_gas):
							new_val = (( u.vector().get_local()[(all_molecules)+k]))
							u_graph_data['conVtime_'+str(k)].append((new_val))
	
						for kGasses in range(0,len(necessary_values['reactants'])):
							sensFuncs[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dP(1)))
	
					elif sens_type == 'total':
	
						pass
	
					else:
						print('Sensitivity analysis is not properly defined.')
						sys.exit()
						
				sens_time_list.append(t)
				progressBar(t, reac_input['Pulse Duration'])
				u_n.assign(u)
				
				t += dt
				constantT.assign(round(t,6))
	
			print()
			print(processTime(start_time))
	
			if reac_input['Sensitivity Analysis'].lower() == 'true':
				if sens_type == 'trans':
					print()
					start_time = time.time()
					print('Evaluating Tape with Tangent Linear Method. Could take some time.')
					tape2.evaluate_tlm()
	
				elif sens_type == 'total':
					dJdm = compute_gradient(jfunc_2, controls)
					djv = [v.values()[0] for v in dJdm]
					print(djv)
					with open('./'+reac_input['Output Folder Name']+'_folder/sensitivity/objSens.txt', 'w') as f:
						f.write("Parameters: "+str(legend_2))
						f.write("Change: "+str(djv))
						f.close
	
				else:
					print('Sensitivity analysis is not properly defined.')
					sys.exit()
	
	
			if reac_input['Sensitivity Analysis'].lower() == 'true':
	
	
				if sens_type == 'trans':
					for numEachSens,eachSens in enumerate(u_graph_data):

						#sys.exit()
						np.savetxt(sensFolder+'/c_'+legend_label[numEachSens]+'.csv',u_graph_data[eachSens],delimiter=",")
	
					for numEachSens,eachSens in enumerate(sensFuncs):
						newList = []
						for kSensNum, kSens in enumerate(sensFuncs[eachSens]):
							newValue = kSens.block_variable.tlm_value
							newList.append(newValue)
						np.savetxt(sensFolder+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
				
				elif sens_type == 'total':
					pass
	
				else:
					print('Sensitivity analysis is not properly defined.')
					sys.exit()
		
			if reac_input['Sensitivity Analysis'].lower() == 'true':
				print(processTime(start_time))
			
			current_reac_set += 1
			if current_reac_set >  pulse_variation:
				current_reac_set = 1
	
			x_values = []
			it_times = []
			j_values = []
			dj_values = []
	
			#############################################################
			############# OPTIMIZATION ITERATION CALLBACKS ##############
			#############################################################
	
			def derivCB(j,dj,m):
				it_times.append(time.time())
				j_values.append(j)
				djv = [v.values()[0] for v in dj]
				dj_values.append(djv)
				mv = [v.values()[0] for v in m]
				x_values.append(mv)
				print("Derivative Time: "+str(time.time() - start_time))
				with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt', 'w') as f:
					f.write("Contents: "+str(it_times))
					f.write('\n')
					f.write("Change: "+str(dj_values))
					f.write('\n')
					f.write("Constants: "+str(x_values))
	
					f.close
				print(j)
				print(djv)
				print(mv)
				#start_time = time.time()
	
			def hessCB(j,dj,m):
				it_times.append(time.time())
				j_values.append(j)
				djv = [v.values()[0] for v in dj]
				dj_values.append(djv)
				mv = [v.values()[0] for v in m]
				x_values.append(mv)
				print(time.time() - start_time)
				with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt', 'w') as f:
					f.write("Contents: "+str(it_times))
					f.write('\n')
					f.write("Change: "+str(dj_values))
					f.write('\n')
					f.write("Constants: "+str(x_values))
	
					f.close
				print(j)
				print(djv)
				print(mv)
	
			#############################################################
			############### SAVE RESULTS FROM THIN-ZONE #################
			#############################################################
	
			if reac_input['Thin-Zone Analysis'].lower() == 'true':
				#print('Storing Catalyst Zone Gas and Surface Data')

				cat_dataRate = {}
				#for j_species in range(0,all_molecules):
				for j_species in range(0,all_molecules):
					timeBetween = time.time()
					new_values = [] 
					mol_values = []
		
					for z in range(0, int(reac_input['Mesh Size'])-1):
	
						if z != 0:
							#value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
							new_values = np.vstack((mol_values,cat_data['rawData'][:,z*(all_molecules)+j_species]))
							mol_values = new_values
					
						else:					
							new_values = cat_data['rawData'][:,z]
							mol_values = cat_data['rawData'][:,z]
					
					cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values),axis=1)
					tempSurface = pd.DataFrame(np.transpose(mol_values[1:,11:]))
					tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/surface_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
						
				#print('Storing Catalyst Zone Gas and Surface Rate Data')
				for j_species in range(0,all_molecules):
						
					rateTest = eval(rateStrings[j_species])
					tempSurface = pd.DataFrame(rateTest)
					tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/surface_data/r_'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)						
						
	
			#############################################################
			################# OPTIMIZATION EXECUTION ####################
			#############################################################
	
			fitting_time = time.time()

			if (reac_input['Fit Parameters'].lower() == 'true' or fit_temperature == True):
	
				start_time = time.time()
				print()
				print('Fitting Kinetic Parameters. Will take some time!')
				

				#print(type(ln(r_const["kf0"]/r_const["kb0"])))
				#print(type(r_const["kf0"]))
				#print(jfunc_2)
				#print(AdjFloat(inner((ln(r_const["kf0"]/r_const["kb0"])+ln(r_const["kf1"]/r_const["kb1"])+ln(r_const["kf2"]/r_const["kb2"])),1)))
				#rf_2 = ReducedFunctional(jfunc_2+(r_const["kf0"]/r_const["kb0"])+(r_const["kf1"]/r_const["kb1"])+(r_const["kf2"]/r_const["kb2"]), controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
				#rf_2 = ReducedFunctional(jfunc_2+AdjFloat((ln(r_const["kf0"]/r_const["kb0"])+ln(r_const["kf1"]/r_const["kb1"])+ln(r_const["kf2"]/r_const["kb2"])) - 1), controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
				#rf_2 = ReducedFunctional(jfunc_2 + assemble(((Constant(8.314*350)*(ln(r_const["kf0"]/r_const["kb0"])+ln(r_const["kf1"]/r_const["kb1"])+ln(r_const["kf2"]/r_const["kb2"]))-w4)**2)*dx()), controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
				print(controls)
				rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
	
				taylorTest = False
	
				if taylorTest == True:
					h = Constant(0.0001)  # the direction of the perturbation
					for jContNum, jCont in enumerate(legend_2): #legend_2
						print('taylor test for parameter '+jCont)
						tTest = ReducedFunctional(jfunc_2, [controls[jContNum]],tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
						conv_rate = taylor_test(tTest, [r_const[jCont]], h)
						conv_rate = taylor_test(tTest, [r_const[jCont]], h, dJdm = 0)
					
					sys.exit()
					
				low_bounds = []
				up_bounds = []
	
				for gt in range(0,len(controls)):
					low_bounds.append(1e-144)
					up_bounds.append(np.inf)
	
				#dj = compute_gradient(jfunc_2, controls)
	
				if reac_input['Optimization Method'] == 'L-BFGS-B' or reac_input['Optimization Method'] == '':
					u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds), tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
					#u_opt_2 = minimize(rf_2, bounds = [low_bounds,up_bounds], tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
				elif reac_input['Optimization Method'] == 'Newton-CG':
					u_opt_2 = minimize(rf_2, method = 'Newton-CG',tol=1e-9, options={"xtol":1e-22})
				elif reac_input['Optimization Method'] == 'BFGS':
					u_opt_2 = minimize(rf_2, method = 'BFGS',tol=1e-13, options={"gtol":1e-13})# , "constraints":bounds
				elif reac_input['Optimization Method'] == 'SLSQP':
					u_opt_2 = minimize(rf_2, method = 'SLSQP', bounds = (low_bounds,up_bounds),tol=1e-9, options={"ftol":1e-9})
				elif reac_input['Optimization Method'] == 'CG':
					u_opt_2 = minimize(rf_2,bounds = (low_bounds,up_bounds), method = 'CG',tol=1e-9, options={"gtol":1e-9})
				elif reac_input['Optimization Method'] == 'basinhopping':
					u_opt_2 = minimize(rf_2, method = 'basinhopping', bounds = (low_bounds,up_bounds),tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
				elif reac_input['Optimization Method'] == 'nonlinear':
					print('non-linear optimization')
					print(low_bounds)
					print(up_bounds)
					problem = MinimizationProblem(rf_2) # ,bounds = (low_bounds,up_bounds)
					ipoptSolver = IPOPTSolver(problem)
					rf_2 = ipoptSolver.solve()
				else:
					print('Requested Optimization Method Does Not Exist')
					sys.exit()
				#sys.exit()
				print(processTime(start_time))
	
				optimization_success = True
	
			#############################################################
			################# CALCULATE HESSIAN FOR UQ ##################
			#############################################################
	
			if reac_input['Uncertainty Quantification'].lower() == 'true':
				start_time = time.time()
				print()
				print('Calculating hessian. Could take some time.')
				
				control = Control(c)
				rf_2 = ReducedFunctional(jfunc_2, control)
	
				dJdm = c._ad_dot(rf_2.derivative())
	
				testHessian = rf_2.hessian(c)
				#print(type(testHessian))
	
				Hm = c._ad_dot(testHessian)
				#print(type(Hm))
				
				#Hm = c._ad_dot(rf_2.hessian(c))
				#Hm = c._ad_dot(rf_2.hessian(c))
				print(processTime(start_time))
				np.savetxt(hessFolder+'/'+reactor_kinetics_input['Sensitivity Parameter']+'_sens.csv',[dJdm],delimiter=",")#
				np.savetxt(hessFolder+'/'+reactor_kinetics_input['Sensitivity Parameter']+'.csv',[Hm],delimiter=",")#
	
			#############################################################
			############# STORE OUTLET FLUX DATA per PULSE ##############
			#############################################################
	
			if k_pulse == 0:
				dictionary_of_numpy_data = {}
				for j_species in range(0,monitored_gas):
					dictionary_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
					new_data = np.asarray(graph_data['timing'])
					dictionary_of_numpy_data[legend_label[j_species]] = np.vstack((dictionary_of_numpy_data[legend_label[j_species]],new_data))
				
			
			#############################################################
			################### PLOT OUTLET FLUX DATA ###################
			#############################################################
	
			if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Display Experimental Data'].lower() == 'true':
	
				if reac_input['Display Objective Points'].lower() == 'true' or reac_input['Fit Parameters'].lower() == 'true':
					if k_pulse == 0:
						for k_fitting in range(0,len(legend_label[:int(len(legend_label))])):
							if objSpecies[k_fitting] == '1':
								ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting])
			
				if reac_input['Display Objective Points'].lower() == 'true':# and reac_input['Fit Inert'].lower() == 'true':
					for k_fitting in range(len(legend_label[:int(len(legend_label))]),len(legend_label)):
						if objSpecies[k_fitting] == '1':
							ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting], alpha=0.3)
	
	
				for k,j in enumerate(graph_data):
					if j != 'timing':
	
						if k_pulse > 0:
							pass
						else:
							dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/species_data/'+legend_label[k]+'.csv',header=None)
							ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
					else:
						pass
					
			#############################################################
			################## ADDITION OF WHITE NOISE ##################
			#############################################################
	
			sig = 0.1
			beta_2 = 0.00270
			w_2 = 2*3.14159*70
	
			for k,j in enumerate(graph_data):
				if j != 'timing':
					if noise == True:
						print('white Noise added')
						for z in range(0,int(reac_input['Time Steps'])):
							graph_data[j][z] += np.random.normal(0,1)*sig +beta_2*np.cos(w_2*(k*dt))
					
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)
					new_data = np.asarray(graph_data[j])
					dictionary_of_numpy_data[legend_label[k]] = np.vstack((dictionary_of_numpy_data[legend_label[k]],new_data))
				
				else:
					pass	
	
			name_list = necessary_values['reactants'][monitored_gas:]
	
			#############################################################
			################## SAVE OUTLET FLUX DATA ####################
			#############################################################
			
			for j_species in range(0,monitored_gas):
				tempDict = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
				
				np.savetxt('./'+reac_input['Output Folder Name']+'_folder/species_data/'+legend_label[j_species]+'.csv', tempDict, delimiter=",")
	
		ax2.legend(title="Gas Species")
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		
		if reac_input['Store Graph'].lower() == 'true':
			plt.savefig('./'+reac_input['Output Folder Name']+'_folder/graphs/species_data.png')
		
		if reac_input['Display Graph'].lower() == 'true':
			plt.show()
		
		for k,j in enumerate(graph_data):
			if j != 'timing':
				ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
			else:
				pass
	
		plt.clf()
		plt.close()
	
		#############################################################
		############## GENERATE OPTIMIZATION GIF ####################
		#############################################################
	
		if reac_input['Fitting Gif'].lower() == 'true':
			if os.path.exists('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt') == True:
				with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt', 'r') as f:
					lines = f.readlines()
				f.close
				lines = [x.strip() for x in lines]
				times = lines[0]
				times = times.replace('Contents: ','')
				times = eval(times)
				constants = lines[2]
				constants = constants.replace('Constants: ','')
				constants = eval(constants)
	
				things = len(times)
				print('things')
				for k_num in range(0,things):
					alter = pd.read_csv(input_file,header=None)
	
					variables_to_change = ['Display Graph','Fit Parameters','Sensitivity Analysis','Store Outlet Flux','Output Folder Name','Reaction_Information']
				
					for k in range(alter.shape[0]):
						if alter[0][k] in variables_to_change:
							if alter[0][k] == 'Store Outlet Flux':
								alter.iloc[k,1] = 'TRUE'
							elif alter[0][k] == 'Output Folder Name':
								alter.iloc[k,1] = reac_input['Output Folder Name']+'_folder/fitting/iter_'+str(k_num)
							elif alter[0][k] == 'Reaction_Information':
								value = 1
								current_rate = k+1
								kz = 0
								while kz < len(kVals):
									if value == 1:
										alter.iloc[current_rate,value] = constants[k_num][kz]
										
										value = 2
										kz+=1
									elif value == 2:
										if str(alter.iloc[current_rate,value]) != 'nan':
											alter.iloc[current_rate,value] = constants[k_num][kz]
											kz+=1
										else:
											pass
										current_rate += 1
										value = 1
	
							else:
								alter.iloc[k,1] = 'FALSE'
					
					alter.to_csv(input_file,header=None,index=False)	
				
					try:
						print()
						print('Iteration: '+str(k_num+1))
						call_sim()
						sys.exit()
					except:
						k_num = things
			
				generateGif(legend_label[:len(legend_label)], reac_input['Experimental Data Folder']+'/species_data', './'+reac_input['Output Folder Name']+'_folder/fitting', len(constants), constants, reactor_kinetics_input['reactions_test'], times)
			
				for k_num in range(0,things):
					shutil.rmtree('./'+reac_input['Output Folder Name']+'_folder/fitting/iter_'+str(k_num)+'_folder') 
	
			user_data = pd.read_csv('./'+reac_input['Output Folder Name']+'_folder/'+input_file,header=None)
			user_data.to_csv(input_file,header=None,index=False)
	
		else:
			return graph_data, legend_label, necessary_values['reactants']
	
	#############################################################
	############## INPUTFILE DICTIONARY READING #################
	#############################################################
	
	def call_sim():
		
		reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readBatchInput(input_file)
		
		reactor_kinetics_input['Fit Parameters'] = 'FALSE'
		reactor_kinetics_input['Display Experimental Data'] = 'FALSE'
		reactor_kinetics_input['Display Objective Points'] = 'FALSE'
		reactor_kinetics_input['Uncertainty Quantification'] = 'FALSE'
		reactor_kinetics_input['Sensitivity Analysis'] =  'TRUE'
		reactor_kinetics_input['Noise'] = 'FALSE'

		if fitting_gif != None:
			reactor_kinetics_input['Fitting Gif'] = 'TRUE'
			reactor_kinetics_input['Display Experimental Data'] = 'TRUE'
		else:
			reactor_kinetics_input['Fitting Gif'] = 'FALSE'

		if noise != 'FALSE':
			reactor_kinetics_input['Noise'] = 'TRUE'
		else:
			reactor_kinetics_input['Noise'] = 'FALSE'			

		if optimization != None:
			reactor_kinetics_input['Fit Parameters'] = 'TRUE'
			reactor_kinetics_input['Optimization Method'] = optimization
			reactor_kinetics_input['Objective Points'] = 'all'

		if uncertainty_quantificaiton != None:
			reactor_kinetics_input['Uncertainty Quantification'] = 'TRUE'
			reactor_kinetics_input['Fit Parameters'] = 'FALSE'
			
		if sensitivityType == None:
			reactor_kinetics_input['Sensitivity Analysis'] =  'FALSE'

		sens_type = sensitivityType
		if (sens_type != 'trans' and sens_type != 'total') and reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true':
			print('loc 1')
			print('Error: sensitivity analysis must be "trans" or "total"')
			sys.exit()

		reactor_kinetics_input['experiment_design'] = None

		if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
			if sens_type == 'trans' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
				for parameters in kinetic_parameters:
					reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readBatchInput(input_file)
					
					reactor_kinetics_input['Sensitivity Parameter'] = parameters
					reactor_kinetics_input['Fit Parameters'] = 'FALSE'
					reactor_kinetics_input['Display Experimental Data'] = 'FALSE'
					reactor_kinetics_input['Fitting Gif'] = 'FALSE'
					reactor_kinetics_input['Uncertainty Quantification'] = 'FALSE'
					reactor_kinetics_input['Display Objective Points'] = 'FALSE'
					reactor_kinetics_input['Sensitivity Analysis'] = 'FALSE'

					if noise != 'FALSE':
						reactor_kinetics_input['Noise'] = 'TRUE'
					else:
						reactor_kinetics_input['Noise'] = 'FALSE'

					if optimization != None:
						reactor_kinetics_input['Fit Parameters'] = 'TRUE'
						reactor_kinetics_input['Optimization Method'] = optimization
						reactor_kinetics_input['Objective Points'] = 'all'
			
					if uncertainty_quantificaiton != None:
						reactor_kinetics_input['Uncertainty Quantification'] = 'TRUE'
						reactor_kinetics_input['Fit Parameters'] = 'FALSE'
					
					if fitting_gif != None:
						reactor_kinetics_input['Fitting Gif'] = 'TRUE'
						reactor_kinetics_input['Display Experimental Data'] = 'TRUE'

					if sensitivityType != None:
						reactor_kinetics_input['Sensitivity Analysis'] =  'TRUE'

					sens_type = sensitivityType					

					print(parameters)
					reactor_kinetics_input['Display Graph'] = 'FALSE'
					
					if reactor_kinetics_input['Fit Parameters'].lower() == 'true':
						print('')
						print('')
					
						print('Running the Sensitivity and Parameter Fitting methods simultaniously is not possible due to conflicts between the tangent linear and adjoint methods.')
						print('Please run again with one of these methods excluded.')
						sys.exit()
					
					if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true':
						print('')
						print('Processing '+parameters)
	
					graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,Ao_in,Ea_in,Ga_in,dG_in,kin_fit,arrForward,arrBackward,gForward)
			
			elif sens_type == 'total':
				graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,Ao_in,Ea_in,Ga_in,dG_in,kin_fit,arrForward,arrBackward,gForward)
		else:
			graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,Ao_in,Ea_in,Ga_in,dG_in,kin_fit,arrForward,arrBackward,gForward)
		
	call_sim()


def run_batch_reactor(timeFunc,includeNoise=False,inputFile = './input_file_batch.csv'):
	
	general_run(timeFunc,input_file = inputFile,noise=includeNoise)	

def fit_batch(timeFunc,optim = 'L-BFGS-B',inputFile = './input_file.csv'):
	general_run(timeFunc,optimization=optim,input_file = inputFile)

def batchGraph(input_file = './input_file.csv',dispExper=False,objectivePoints=False,displayGraph=True,storeGraph=False,outputName='./flux.png'):
	timeFunc = 0.2

	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readBatchInput(input_file)
	#reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file)

	reac_input = reactor_kinetics_input
	reac_input['Pulse Duration'] = timeFunc
	reac_input['Infinite Inert'] = 'FALSE'
	reac_input['Display Experimental Data'] = 'False'

	if dispExper != False:
		reac_input['Display Experimental Data'] = 'TRUE'
		reac_input['Display Objective Points'] = 'TRUE'
		reac_input['Objective Points'] = 'all'

	if objectivePoints != False:
		reac_input['Display Objective Points'] = 'TRUE'
		reac_input['Objective Points'] = 'all'
	print(reac_input['Pulse Duration'])
	#sys.exit()
	reac_input['Time Steps'] = reac_input['Pulse Duration']*1000
	dt = 0.001

	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_batch_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],arrForward,arrBackward,gForward,fit_temperature)

	fig2,ax2,legend_label,header,colors = establishOutletGraph('batch',necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],None,'FALSE')

	monitored_gas = necessary_values['molecules_in_gas_phase']
	reac_input['Scale Output'] = 'FALSE'

	reac_input['Objective Species'] = '1,1,1,0,0'

	if '/' in str(reac_input['Pulse Size']):
		list_of_feed = reac_input['Pulse Size'].split('/')
		list_of_time = reac_input['Pulse Time'].split('/')
		pulse_variation = len(list_of_feed)
		pulse_time = len(list_of_time)
		reactant_feed = []
		reactant_time = []
		for k in range(0,len(reac_input['Pulse Size'].split('/'))):
			reactant_feed.append(list_of_feed[k].split(','))
			reactant_time.append(list_of_time[k].split(','))
	else:
		pulse_variation = 1
		reactant_feed = []
		reactant_time = []
		
		try:
			reactant_feed.append(reac_input['Pulse Size'].split(','))
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
		except AttributeError:
			reactant_feed.append(reac_input['Pulse Size'])
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
	
	species_pulse_list = reactant_feed[0]
	current_reac_set = 1
	if type(reactant_time[current_reac_set-1]) == list:
		species_time = reactant_time[current_reac_set-1]
	else:
		species_time = reactant_time

	if str(reac_input['Objective Species']).find(',') != -1:
		objSpecies = list(reac_input['Objective Species'].split(','))
	else:
		print(reac_input['Objective Species'])
		objSpecies = [str(int(reac_input['Objective Species']))]

	if dispExper != False:
		output_fitting = curveFitting(legend_label[:int(len(legend_label))],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)

	dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])

	constantTemp = reac_input['Reactor Temperature']
	constantTemp = Constant(constantTemp)
	testTemp = reac_input['Reactor Temperature']
	testTemp = Constant(testTemp)

	if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Display Experimental Data'].lower() == 'true':

		for k,j in enumerate(legend_label):
			if j != 'timing': #and j != 'conVtime_1' and j != 'conVtime_2' and j != 'conVtime_3' and j != 'conVtime_0' 
				dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/species_data/'+legend_label[k]+'.csv',header=None)
				ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
			else:
				pass

	for knum,k in enumerate(necessary_values['reactants']):
		if knum < necessary_values['molecules_in_gas_phase']:

			dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/species_data/'+k+'.csv',header=None)
			graphLimit = dfTemp[0][len(dfTemp[0]) - 1]
			
			legendInclude = True
			for j in range(1, len(dfTemp.columns)):
				if legendInclude == True:
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label=legend_label[knum], alpha=0.7)
					legendInclude = False
				else:
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label='_nolegend_', alpha=0.7)
	ax2.set_xlim(0,graphLimit)
	
	ax2.legend(title="Gas Species",loc = 1)

	if storeGraph == True:
		plt.savefig(outputName)

	if displayGraph == True:
		plt.show()
