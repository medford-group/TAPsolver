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

def general_run(timeFunc,uncertainty_quantificaiton=None,optimization=None,fitting_gif=None,pulseNumber=1,store_flux_func='TRUE',store_thin_func='FALSE',input_file = './input_file.csv',sensitivityType=None,noise='FALSE',fitInert=None,inputForm='old'):
	
	### How much concentration data do you want stored for 'thin-zone analysis'
	# 'point' = only center point in the catalyst zone
	# 'average' = average concentration values along the catalyst zone
	# 'cat' = all points in the catalyst zone
	# 'all' = all points in the reactor
	
	sampling = False

	thinSize = 'cat'
	#thinSize = 'point'
	#thinSize = 'average'
	#thinSize = 'all'
	
	### Fitting the temperature?
	
	# Fit = True
	# Don't fit = False
	
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
		reac_input['Thin-Zone Analysis'] = store_thin_func 
		reac_input['Time Steps'] = reac_input['Pulse Duration']*1000
	
		sens_type = sensitivityType
		if (sens_type != 'trans' and sens_type != 'total') and reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true':
			print('Error: sensitivity analysis must be "trans" or "total"')
			sys.exit()

		reac_input['Optimization Method'] = 'L-BFGS-B' #
		reac_input['Objective Points'] = 'all'
		reac_input['Store Outlet Flux'] = 'TRUE'
		reac_input['Store Graph'] = 'FALSE'
		reac_input['Display Graph'] = 'FALSE'
		reac_input['Scale Output'] = 'FALSE'
		reactor_kinetics_input['Reactor Type'] = 'tap'

		if reac_input['Advection Value'] > 0.0:
			reac_input['Advection'] = 'true'
		else:
			reac_input['Advection'] = 'false'

		if 'thermo equations' in reac_input.keys():
			if len(reac_input['thermo equations']) > 0:
				thermoConstraints = True

			else:
				thermoConstraints = False

		else:
			thermoConstraints = False

		if thermoConstraints == True:
			
			for z in reac_input['thermo equations']:
				thermoReactions = []
				thermoStoich = []
				tempValue = z.split('+')
				for j in tempValue:
					if j.find('*') != -1:
						thermoSplit = j.split('*')
						thermoSplit[1] = thermoSplit[1].replace('r','')
						thermoReactions.append(int(thermoSplit[1]))
						thermoSplit[0] = thermoSplit[0].replace('(','')
						thermoSplit[0] = thermoSplit[0].replace(')','')
						thermoStoich.append(float(thermoSplit[0]))
					else:
						j = j.replace('r','')
						thermoReactions.append(int(j))
						thermoStoich.append(1)
	
		path = './'+reac_input['Output Folder Name']+'_folder/'
		generateFolder(path)
			
		path_3 = './'+reac_input['Output Folder Name']+'_folder/flux_data/'
		generateFolder(path_3)
		
		if reac_input['Thin-Zone Analysis'].lower() == 'true':
			path_4 = './'+reac_input['Output Folder Name']+'_folder/thin_data/'
			generateFolder(path_4)
	
		path_5 = './'+reac_input['Output Folder Name']+'_folder/graphs/'
		generateFolder(path_5)
	
		if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or fit_temperature == True:
			path_6 = './'+reac_input['Output Folder Name']+'_folder/fitting/'
			generateFolder(path_6)
	
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
			rangesurface_species = reac_input['Initial Surface Composition']
		
		# Initialize the grid system, time step size, pulse size
		r_param, dx_r, dx2_r, frac_length, cat_location = establishMesh(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['Mesh Size'])
		cat_location = 1 - reac_input['Catalyst Location']
		dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
		eb = np.array((reac_input['Void Fraction Inert'],reac_input['Void Fraction Catalyst'],reac_input['Void Fraction Inert']))
		ca = (reac_input['Reactor Radius']**2)*3.14159 
		point_volume = dx_r * ca * eb[0]
		Inert_pulse_conc = reac_input['Reference Pulse Size']/(point_volume)
		
		# Define the diffusion constants
		ref_rate = np.append(reac_input['Reference Diffusion Inert'],reac_input['Reference Diffusion Catalyst'])
		ref_rate = np.append(ref_rate,ref_rate[0])  	
		D = np.empty((len(reac_input['Mass List'].split(',')),3))
		
		def diff_func(ref_mass,ref_T,mol_mass,ref_r):
			return ref_r*(mp.sqrt(ref_mass*reac_input['Reactor Temperature'])/mp.sqrt(ref_T*mol_mass))
		Dout = []
		Din = []
		constantTemp = reac_input['Reactor Temperature']
		constantTemp = Constant(constantTemp)
		testTemp = reac_input['Reactor Temperature']
		testTemp = Constant(testTemp)
	
		if fit_temperature == True:
			controls = []
			new = Control(constantTemp)
			controls.append(new)
	
		compMass_list = reac_input['Mass List'].split(',')
	
		for k,j in enumerate(reac_input['Mass List'].split(',')):
			for k_2,j_2 in enumerate(ref_rate):
				D[k,k_2] = diff_func(reac_input['Reference Mass'],reac_input['Reference Temperature'],float(j),j_2) ###??? should this (np.sum(r_param)**2) be here? This puts it in dimensional form!
							
				if k_2 == 0:
					Dout.append(Constant(D[k,k_2]))
				if k_2 == 1:
					Din.append(Constant(D[k,k_2]))
	
		# Construct the rate expression based on the microkinetic model
		necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,fit_temperature)
		
		if thermoConstraints == True:

			thermoReactants = []
			thermoProducts = []

			thermoReactantsValue = []
			thermoProductsValue = []

			for jnum,jval in enumerate(thermoReactions):
				if rev_irr[jval-1] != 1:
					print('')
					print('To set a thermodynamic constraints, all reactions must be reversible.')
					print('Currently, elementary reaction '+str(jval)+' is irreversible and included in the thermo. equation.')
					sys.exit()

				else:
					for zen in range(0,necessary_values['molecules_in_gas_phase']):
						if rate_array[jval-1][zen] == 1:
							thermoProductsValue.append(necessary_values['reactants'][zen])						
							if thermoStoich[jnum] != 1:
								thermoProducts.append('('+str(thermoStoich[jnum])+')*'+necessary_values['reactants'][zen])
							else:
								thermoProducts.append(necessary_values['reactants'][zen])	

						elif rate_array[jval-1][zen] == -1:
							thermoReactantsValue.append(necessary_values['reactants'][zen])
							if thermoStoich[jnum] != 1:
								thermoReactants.append('('+str(thermoStoich[jnum])+')*'+necessary_values['reactants'][zen])
							else:
								thermoReactants.append(necessary_values['reactants'][zen])		
			
			overallReactionString = ''
			
			for jnum in range(0,len(thermoReactants)):
				overallReactionString = overallReactionString+thermoReactants[jnum]
				if jnum != len(thermoReactants)-1:
					overallReactionString = overallReactionString + ' + '
			overallReactionString = overallReactionString + ' <-> '
			
			for jnum in range(0,len(thermoProducts)):	
				overallReactionString = overallReactionString + thermoProducts[jnum]
				if jnum == len(thermoReactants)-1:
					overallReactionString = overallReactionString + ' + '

			print('Thermodynamic Constraint Equation:')
			print(overallReactionString)
			print('')
			if len(reactor_kinetics_input['thermo values']) == 0: 
				freeEnergyValue = 0

				try:
					for jnum in range(0,len(thermoReactantsValue)):
						freeEnergyValue += (-1)*(thermoStoich[jnum])*molecularProperties(thermoReactantsValue[jnum],'freeEnergy',temperature=reac_input['Reactor Temperature'])

					for jnum in range(0,len(thermoProductsValue)):
						freeEnergyValue += (thermoStoich[jnum])*molecularProperties(thermoProductsValue[jnum],'freeEnergy',temperature=reac_input['Reactor Temperature'])
				except:
					print('Failed to calculate free energy of reaction automatically. Please enter free energy manually.')
					sys.exit()

				reactor_kinetics_input['thermo values'].append(freeEnergyValue)
				print('Free energy calculated automatically: ')
				print(reactor_kinetics_input['thermo values'][0])

			else:
				print('Free energy entered manually: ')
				print(reactor_kinetics_input['thermo values'][0])

		#print(necessary_values['molecules_in_gas_phase'])

		#print(thermoReactions)
		#print(thermoStoich)
		#print(rate_array)
		#print(rev_irr)
		#print(necessary_values['molecules_in_gas_phase'])
		#print(necessary_values['reactants'])

		#sys.exit()

		if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or (reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' and sens_type == 'total') or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='1,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='1'
			for k in range(0,int(reac_input['Number of Inerts'])):
				if k != int(reac_input['Number of Inerts']):
					speciesString+='0,'
				else:
					speciesString+='0'

		elif reactor_kinetics_input['Fit Inert'].lower() == 'true':
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='0,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='0'
			for k in range(0,int(reac_input['Number of Inerts'])):
				if k != int(reac_input['Number of Inerts']):
					speciesString+='1,'
				else:
					speciesString+='1'

		else:
			speciesString = ''
			for k in range(0,int(necessary_values['molecules_in_gas_phase'])):
				if k != int(necessary_values['molecules_in_gas_phase']):
					speciesString+='0,'
				elif int(reac_input['Number of Inerts']) != 0:
					speciesString+='0'
			for k in range(0,int(reac_input['Number of Inerts'])):
				if k != int(reac_input['Number of Inerts']):
					speciesString+='0,'
				else:
					speciesString+='0'

		reactor_kinetics_input['Objective Species'] = speciesString


		#############################################################
		######## INITIALIZATION OF FINITE ELEMENTS AND MESH #########
		#############################################################
	
		# Define the base mesh (the initial uniform mesh)
		mesh = UnitIntervalMesh(int(reac_input['Mesh Size']))
	
		# Refine the mesh (depending on user specification)
		catalystRefinement = int(reac_input['Catalyst Mesh Density'])
		cfDict = {}
	
		roundedMesh2 = ((1-cat_location) + 0.5*frac_length)*reac_input['Mesh Size']	
		Mesh2 = round(((1-cat_location) + 0.5*frac_length)*reac_input['Mesh Size'])
	
		roundedMesh1 = ((1-cat_location) - 0.5*frac_length)*reac_input['Mesh Size']
		Mesh1 = round(((1-cat_location) - 0.5*frac_length)*reac_input['Mesh Size'])
	
		if Mesh2 != roundedMesh2 or Mesh1 != roundedMesh1:
			print('Warning: Catalyst zone will be refined and rounded to the nearest whole mesh point!')
			trueMesh = (roundedMesh2 - roundedMesh1)/reac_input['Mesh Size']
			newMesh = (Mesh2 - Mesh1)/reac_input['Mesh Size']
			print()
			print('New Catalyst Fraction = '+str(newMesh))
			print('Old Catalyst Fraction = '+str(trueMesh))
			percentChange = abs(round(100*(trueMesh - newMesh)/trueMesh,2))
			print('Change = '+str(percentChange)+'%')
			print()
			if percentChange > 4:
				print('Consider refining the mesh to improve the accuracy of the simulation!')
				sys.exit()
			
		for jayz in range(0,catalystRefinement+1):
			class thin_zoneTest(SubDomain):
				def inside(self, x, on_boundary):
					return between(x[0], (((1-cat_location) - 0.5*frac_length), ((1-cat_location) + 0.5*frac_length)))
			
			thin_zoneTest = thin_zoneTest()
			cfDict[jayz] = MeshFunction("bool",mesh,1)
			
			thin_zoneTest.mark(cfDict[jayz],jayz)
			mesh = refine(mesh,cfDict[jayz],True)
	
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
		u_n = Function(V)
		u_temp = Function(V)
		
		if reac_input['Advection'].lower() == 'true':
			W = VectorFunctionSpace(mesh, 'P', 1)
			advTerm = Function(W)
			advTerm.vector()[:] = reac_input['Advection Value']
	
		graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data = initializeVariableDictionaries(necessary_values,all_molecules,V,u,u_n)
	
		meshCells = int((cat_location + 0.5*frac_length)*reac_input['Mesh Size']) - mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])
		frac_temp = (meshCells*2**(int(reac_input['Catalyst Mesh Density'])))/(int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells)
		totalNumCells = int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells
	
		# Define subdomains (thin-zone / final cell)
		class thin_zone(SubDomain):
			def inside(self, x, on_boundary):
				return between(x[0], (((1-cat_location) - 0.5*frac_length)-1/reac_input['Mesh Size'], ((1-cat_location) + 0.5*frac_length)+1/reac_input['Mesh Size']))
				
		class singlePoint(SubDomain):
			def inside(self,x,on_boundary):
				return between(x[0], (((1-cat_location) - 1/totalNumCells), ((1-cat_location) + 1/totalNumCells)))
	
		thin_zone = thin_zone()
		domains = MeshFunction("size_t", mesh,1)
		thin_zone.mark(domains,1)
	
		newBoundaries = domains.mesh().coordinates().transpose().tolist()[0]
		additionalCells = len(newBoundaries[int(reac_input['Mesh Size']+1):])
	
		centralPoint = singlePoint()
		boundary_parts0 = MeshFunction("size_t", mesh,0)
		boundary_parts0.set_all(0)
		centralPoint.mark(boundary_parts0, 1)
	
		dx = Measure("dx",subdomain_data=domains)
		dT = Measure("dx",subdomain_data=boundary_parts0)
		
		# Define inlet and outlet boundaries
		def boundary_L(x, on_boundary):
			return on_boundary and near(x[0],0,tol)
		
		def boundary_R(x, on_boundary):
			return on_boundary and near(x[0],1,tol)
		
		dz = 1/reac_input['Mesh Size']
	
		class integration_section(SubDomain):
			def inside(self, x, on_boundary):
				return between(x[0], (1-dz,1.0))
		
		right = CompiledSubDomain("near(x[0], 1.)")
	
		boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
		right.mark(boundary_parts, 1)
	
		monitored_gas = necessary_values['molecules_in_gas_phase']
		
		bcs = defineBCs(reac_input['Reactor Type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['Number of Reactants'],all_molecules,reac_input['Pulse Size'],boundary_L,boundary_R,reac_input['Number of Inerts'])
		
		fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])
		ax2.set_xlim(0,reac_input['Pulse Duration'])
	
		# Define controls for optimization or sensitivity analysis
		if reac_input['Fit Inert'].lower() == 'true':
			controls = []
			legend_2 = []
			for j in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
				print(j)
				controls.append(Control(Dout[j]))
				legend_2.append(j)
	
		if str(reac_input['Objective Species']).find(',') != -1:
			objSpecies = list(reac_input['Objective Species'].split(','))
		else:
			print(reac_input['Objective Species'])
			objSpecies = [str(int(reac_input['Objective Species']))]
	
		t = Constant(0)
		uv = Expression('(1/sqrt(3.14159))*exp(-x[0]*100*t)',t=0,degree=1)
	
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
	
		#############################################################
		##### DEFINE METHOD OF OPTIMIZATION (OBJECTIVE FUNC.) #######
		#############################################################
		
		if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Fit Inert'].lower() != 'true' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or sampling == True or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True or reac_input['Fitting Gif'].lower() == True:
			
			if reac_input['Uncertainty Quantification'].lower() == 'true':
				print("Uncertainty Quantification")
			try:
				if type(reac_input['Objective Points']) == float:
					output_fitting = pointFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
				elif reac_input['Objective Points'] == 'all':
					output_fitting = curveFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
					
				else:
					print('Objective Points defined incorrectly')
					sys.exit()
			except TypeError:
				print('Objective Point Input Is Not Valid')
				sys.exit()
	
		if reac_input['Fit Inert'].lower() == 'true':# or reac_input['Display Objective Points'].lower() == 'true':
			try:
				if type(reac_input['Objective Points']) == int:
					output_fitting = pointFitting(legend_label[int(reac_input['Number of Inerts'])+1:],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'])
				elif reac_input['Objective Points'] == 'all':
					
					output_fitting = curveFitting(legend_label[(len(legend_label) - int(reac_input['Number of Inerts']) ) :],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies[(len(legend_label) - int(reac_input['Number of Inerts']) ) :])
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
			
		to_flux = fluxGeneration(reac_input['Reactor Type'],len(legend_label),reac_input['Number of Reactants'],reac_input['Reference Pulse Size'],D,eb,dx_r,reac_input['Reactor Radius'],dx2_r,reac_input['Scale Output'])
	
	
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
	
			
			problem = NonlinearVariationalProblem(F,u,bcs,J)
			
			problem.set_bounds(ninfty,pinfty)
	
			solver = NonlinearVariationalSolver(problem)
			solver.parameters.update(snes_solver_parameters)
		
		# Define a newton variational problem solver
		elif reac_input['Variational Solver'].lower() == 'newton':
			problem = NonlinearVariationalProblem(F,u,bcs,J)
			solver = NonlinearVariationalSolver(problem)
	
			problemtemp = NonlinearVariationalProblem(Ftemp,u,bcs,Jtemp)
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
	
		if reac_input['Thin-Zone Analysis'].lower() == 'true':
			rateStrings = rateEqs(rate_array,rev_irr,gForward,arrForward,arrBackward)
	
		if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or sampling == True or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True:
			for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
				if objSpecies[k_fitting] == '1':
					for timeStep in range(0,len(output_fitting[legend_label[k_fitting]]['times'])):
						output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)
	
		if reac_input['Fit Inert'].lower() == 'true':
			for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts']+1)]),len(legend_label)):
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
	
		for k_pulse in range(0,int(reac_input['Number of Pulses'])):
	
			start_time = time.time()
			
			if sampling == False:
				print("")
				print("Simulation Status: "+"Pulse #"+str(k_pulse+1))
	
			species_pulse_list = reactant_feed#reactant_feed[current_reac_set-1]
			
			if type(reactant_time[current_reac_set-1]) == list:
				species_time = reactant_time[current_reac_set-1]
			else:
				species_time = reactant_time
	
			if reac_input['Knudsen Test'].lower() == 'true':
				knudsenTest(legend_label[(int(len(legend_label))-int(reac_input['Number of Inerts'])):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],reac_input['Reference Pulse Size'],species_pulse_list[(int(len(legend_label))-int(reac_input['Number of Inerts'])):])
	 
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
			
			for kjc in range(0,int(reac_input['Number of Inerts'])):
				graph_data['conVtime_'+str((necessary_values['gas_num']+1)+necessary_values['surf_num']+kjc)] = []
			
			surf_data['timing'] = []
			
			for j_gasses in range(necessary_values['molecules_in_gas_phase'],len(necessary_values['reactants'])-1):
				surf_data['conVtime_'+str(j_gasses)] = []
			
			graph_data['timing'] = []
	
			cat_data['timing'] = []
			
			for z_gasses in range(0,all_molecules+1):
				cat_data['conVtime_'+str(z_gasses)] = []
	
			cat_data['rawData'] = []
			
			t = 0
			if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or sampling == True or fit_temperature == True:
				osub = integration_section()
				domains = MeshFunction("size_t", mesh,0)
				domains.set_all(0)
				osub.mark(domains, 1)
				dP = Measure('vertex',domain = mesh, subdomain_data=domains)
			
			test_new = project(-u.dx(0),V)
			test_new_split = test_new.split(deepcopy=False)
	
			w_new = Expression("1", degree=0)
			w_new2 = interpolate(w_new,V_du)
	
			x_dim = list(range(0, int(reac_input['Mesh Size'])+1))
			
			cum_sum = 0
	
			meshCells = int((cat_location + 0.5*frac_length)*reac_input['Mesh Size']) - mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']) 
	
			#############################################################
			############ SOLVE PDEs AT EACH TIME STEP ###################
			#############################################################
	
			while t <= reac_input['Pulse Duration']:
				graph_data['timing'].append(t)
				for k in range(0,monitored_gas):
					new_val = (to_flux[k]*( u_n.vector().get_local()[(all_molecules)+k]))
					graph_data['conVtime_'+str(k)].append((new_val))
	
				for kjc in range(0,int(reac_input['Number of Inerts'])):
					new_val = ((to_flux[monitored_gas+kjc]*u_n.vector().get_local()[2*(all_molecules+1)-2-(int(reac_input['Number of Inerts'])-kjc)]))
					graph_data['conVtime_'+str(all_molecules-(int(reac_input['Number of Inerts'])-kjc))].append((new_val))
				
				mesh_size = reac_input['Mesh Size']
	
				#if reac_input['Maximization Process'].lower() == 'true':
	
				#	for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
				#		if objSpecies[k_fitting] == '1':
				#			#if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
	
				#			try:
				#				if legend_label[k_fitting] != 'Inert':
				#					jfunc_2 -= assemble(inner(u_n[k_fitting]*to_flux[k_fitting],u_n[k_fitting]*to_flux[k_fitting])*dP(1))							
				#				else:
				#					pass
	
				#			except UnboundLocalError:
				#				if legend_label[k_fitting] != 'Inert':
				#					jfunc_2 = -assemble(inner(u_n[k_fitting]*to_flux[k_fitting],u_n[k_fitting]*to_flux[k_fitting])*dP(1))
				#				else:
				#					pass
	
				#############################################################
				######## STEP FOR CONSTRUCTING OBJECTIVE FUNCTION ###########
				#############################################################
	
				if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or sampling == True or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True:
	
					# Define the objective function 	
					objectiveAnalysis = True
					selectivity = False
					conversion = False
					yieldValue = False
	
					if objectiveAnalysis == True:
						for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
						
							if objSpecies[k_fitting] == '1':
															
								if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
										
									c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
									slope = (-c_exp)/(1/mesh_size)
									intercept = c_exp - ((1-(1/mesh_size))*slope)
									w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)
									w_new2 = interpolate(w_new,V_du)
									w3 = project(w_new2,V_du)
	
									try:
										if legend_label[k_fitting] != 'Inert':
											jfunc_2 += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))							
										else:
											pass
	
									except UnboundLocalError:
										if legend_label[k_fitting] != 'Inert':
											w_temp = Expression(str(1),degree=0) # str(freeEnergyValue)
											w_temp2 = interpolate(w_temp,V_du)
											w4 = project(w_temp2,V_du)		
											thermoWeight = 1e-1#1e-2									
											jfunc_2 = assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
											if thermoConstraints == True:

												for jnum,jval in enumerate(thermoReactions):
													if jnum == 0:
														tempFunc = thermoStoich[jnum]*(r_const["kf"+str(jval-1)]/r_const["kb"+str(jval-1)])
													else:
														tempFunc += thermoStoich[jnum]*(-8.314*reac_input['Reactor Temperature'])*(r_const["kf"+str(jval-1)]/r_const["kb"+str(jval-1)])
												tempFunc = tempFunc*thermoWeight
												jfunc_2 += assemble(inner(tempFunc-w4,tempFunc-w4)*dx())

										else:
											pass
	
					if selectivity == True:
						w_new = Expression('1',degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)
	
						try:
							if legend_label[2] != 'Inert':
								selectDenom = 0
								for zap in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
									selectDenom += u_n[zap]*to_flux[zap]
								jfunc_2 += assemble(inner(u_n[2]*to_flux[2]/selectDenom,w3)*dP(1))							
							else:
								pass
	
						except UnboundLocalError:
							if legend_label[2] != 'Inert':
								selectDenom = 0
								for zap in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
									selectDenom += u_n[zap]*to_flux[zap]
								jfunc_2 = assemble(inner(u_n[2]*to_flux[2]/selectDenom,w3)*dP(1))
							else:
								pass
						
	
					if conversion == True:
											
						w_new = Expression('1',degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)
	
						try:
							if legend_label[2] != 'Inert':
								jfunc_2 += assemble(inner(dt*u_n[2]*to_flux[2],w3)*dP(1))							
							else:
								pass
	
						except UnboundLocalError:
							if legend_label[2] != 'Inert':
								jfunc_2 = assemble(inner(dt*u_n[2]*to_flux[2],w3)*dP(1))
							else:
								pass
	
					if yieldValue == True:
											
						w_new = Expression('1',degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)
	
						try:
							if legend_label[2] != 'Inert':
								jfunc_2 += assemble(inner( (1/(dt*u_n[2]*to_flux[2])) -  (u_n[0]*to_flux[0]/(u_n[2]*to_flux[2]))  ,w3)*dP(1))
							else:
								pass
	
						except UnboundLocalError:
							if legend_label[2] != 'Inert':
								jfunc_2 = assemble(inner( (1/(dt*u_n[2]*to_flux[2])) -  (u_n[0]*to_flux[0]/(u_n[2]*to_flux[2]))  ,w3)*dP(1))
							else:
								pass
	
				if reac_input['Fit Inert'].lower() == 'true':
	
					for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
						if objSpecies[k_fitting] == '1':
							if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
								c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
								slope = (-c_exp)/(1/mesh_size)
								intercept = c_exp - ((1-(1/mesh_size))*slope)
								w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)
								w_new2 = interpolate(w_new,V_du)
								w3 = project(w_new2,V_du)
	
								try:
									jfunc_2 += assemble(inner(u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3,u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3)*dP(1))								
									
								except UnboundLocalError:
									jfunc_2 = assemble(inner(u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3,u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3)*dP(1))
	
				if round(t,6) not in species_time:
					#timeDiff = round(t,6)
					try:
						if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
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
							
							#if len(cat_data['rawData'].shape) > 1:
							#	print('greater')
							arrayNew = np.vstack((cat_data['rawData'],u_n.vector().get_local()))
							#else:	
							#	print('less')
							#arrayNew = np.stack((cat_data['rawData'],u_n.vector().get_local()[:,0]),axis=1)
							cat_data['rawData'] = arrayNew
							#np.append(cat_data['rawData'],u_n.vector().get_local())
							#cat_data['rawData'].append(u_n.vector().get_local())
	
							#for jk in range(0,all_molecules+1):
							#	new_values = [] 
							#	cat_values = []
							#
							#	#for z in range(0,100):
							#	#### for z in range(4500,int((reac_input['Mesh Size'])+1)+additionalCells):
							#	#for z in range(0,100):
							#	#for z in range(4000,int((reac_input['Mesh Size'])-1)+additionalCells+3):
							#	#print(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1)
							#	#print((int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1)+additionalCells)
							#	#sys.exit()
							#	#for z in range(0,int((reac_input['Mesh Size'])+1)+additionalCells):
							#	
							#	#for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1,int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])):
							#	if additionalCells == 0:
							#		#print(mp.int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1)
							#		#print(int((cat_location + 0.5*frac_length)*reac_input['Mesh Size']))
							#		#sys.exit()
							#		for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
							#		#for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']),int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])):
							#			try:
							#				#value_cat = u_n.vector().get_local()[z*(all_molecules)-(all_molecules)+jk]
							#				value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
							#				cat_values.append(value_cat)
							#	
							#			except AttributeError:
							#				value_cat = u_n.vector().get_local()[z*(all_molecules)]
							#				sys.exit()
							#				#cat_values[0] += value_cat*dx_r
							#				cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))
							#	else:
							#		#print(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']))
							#		#print(int((cat_location + 0.5*frac_length)*reac_input['Mesh Size']))
							#		#sys.exit()
							#		for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
							#
							#		#for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']),int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))):
							#			#print(z)
							#			try:
							#				value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
							#				cat_values.append(value_cat)
							#	
							#			except AttributeError:
							#				value_cat = u_n.vector().get_local()[z*(all_molecules)]
							#				sys.exit()
							#				cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))
							#		###sys.exit()
							#
							#	cat_data['conVtime_'+str(jk)].append(cat_values)
	
					except RuntimeError:
						print('Time Step Failure')
						sys.exit()
					
				else:
					dt = 0.0001
					dk.assign(dt)
	
					if round(t,6) in species_time:
						timeDiff = round(t,6)
	
					if reac_input['Reactor Type'] == 'tap':
						
						if reac_input['Fit Parameters'].lower() == 'true':
							if t == 0:
								for k in range(0,int(reac_input['Number of Reactants'])):
									u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+additionalCells))+k-(all_molecules)] = -1e-10
	
						#if ',' in str(species_pulse_list):
						
						for k in range(0,int(reac_input['Number of Reactants'])):
								
							if species_time[k] == round(t,6):
								u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells)))+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]
									
						#else:
						#	u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']-10)-1)-(all_molecules)] = float(species_pulse_list)*Inert_pulse_conc
	
						for k_step in range(0,int(reac_input['Number of Inerts'])):
							if species_time[-1-k_step] == round(t,6):
								u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells))-1-k_step)] = float(species_pulse_list[-1-k_step])*Inert_pulse_conc###??? Added the porosity contribution
								
						#############################################################
						################### STORE THIN-ZONE DATA ####################
						#############################################################
	
						if reac_input['Thin-Zone Analysis'].lower() == 'true':
							
							if round(t,6) != 0:
								arrayNew = np.vstack((cat_data['rawData'],u_n.vector().get_local()))
								cat_data['rawData'] = arrayNew
							else:
								cat_data['rawData'] = u_n.vector().get_local()
	
							#for jk in range(0,all_molecules+1):
								#new_values = [] 
								#cat_values = []
	
								#if additionalCells == 0:
								#	for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
							
								#		try:
								#			value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
								#			cat_values.append(value_cat)
								#		except AttributeError:
								#			value_cat = u_n.vector().get_local()[z*(all_molecules)]
								#			sys.exit()
								#			
								#			cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))
								#else:
								#	for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
								#		try:
								#			value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
								#			cat_values.append(value_cat)
								#		except AttributeError:
								#			value_cat = u_n.vector().get_local()[z*(all_molecules)]
								#			sys.exit()
								#			cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))							
	
								#cat_data['conVtime_'+str(jk)].append(cat_values)
	
					#############################################################
					############### DEFINE INITIAL COMPOSITION ##################
					#############################################################
	
					if k_pulse == 0 and round(t,6) == 0:
						if additionalCells == 0:
	
							for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+1,int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])+2):
								if ',' in str(reac_input['Initial Surface Composition']):
									for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
										if int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1 <= 0:
											u_n.vector()[(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
										else:
											u_n.vector()[z*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
	
								else:
									if int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1 <= 0:
										u_n.vector()[(all_molecules)-(2)] = float(species_pulse_list)
									else:
										u_n.vector()[z*(all_molecules)-(2)] = float(species_pulse_list)
	
						else:
	
							transTest1 = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+1
							if ',' in str(reac_input['Initial Surface Composition']):
								for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
									u_n.vector()[transTest1*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)*0.5
	
							transTest2 = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))
							if ',' in str(reac_input['Initial Surface Composition']):
								for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
									u_n.vector()[transTest2*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)*0.5
	
	
							for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+2,int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))):
								if ',' in str(reac_input['Initial Surface Composition']):
									for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
										if int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1 <= 0:
											u_n.vector()[(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
										else:
											u_n.vector()[z*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
	
								else:
									if int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1 <= 0:
										u_n.vector()[(all_molecules)-(2)] = float(species_pulse_list)
									else:
										u_n.vector()[z*(all_molecules)-(2)] = float(species_pulse_list)
						
					try:
						if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
							solvertemp.solve()
	
						else:
							solvertemp.solve(annotate = False)
	
					except RuntimeError:
						print('Time Step Failure')
						sys.exit()
				
				#############################################################
				################ STORE SENSITIVITY DATA #####################
				#############################################################
	
				if reac_input['Sensitivity Analysis'].lower() == 'true':
					
					if sens_type == 'trans':
	
						u_final = u.split(deepcopy=False)
						should_it = int(round(t*reac_input['Time Steps']/reac_input['Pulse Duration'],0))
					
						for k in range(0,monitored_gas):
							new_val = (( to_flux[k]*u.vector().get_local()[(all_molecules)+k]))
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
						print(u_graph_data[eachSens])
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
				print('Storing Catalyst Zone Gas and Surface Data')
				if int(reac_input['Number of Pulses']) == 1:
					cat_dataRate = {}
					#for j_species in range(0,all_molecules):
					for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
						timeBetween = time.time()
						new_values = [] 
						mol_values = []
		
						for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
	
							if z != 0:
								#value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
								new_values = np.vstack((mol_values,cat_data['rawData'][:,z*(all_molecules)+j_species]))
								mol_values = new_values
					
							else:					
								new_values = cat_data['rawData'][:,z] 
								mol_values = cat_data['rawData'][:,z]
						
						top = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+1
						bottom = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-1
	
						## Change so that only thin-zone data is being stored
						if thinSize == 'cat':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values[top:bottom,:]),axis=1)
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							#tempSurface.to_csv(,header=None,index=False)
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)	
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
				
						elif thinSize == 'point':
							centerPoint = int((top + bottom)/2)
							cat_dataRate['convtime_'+str(j_species)] = np.transpose(mol_values[centerPoint,:])
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.transpose(mol_values[centerPoint,:]), delimiter=",")
				
						elif thinSize == 'average':
							cat_dataRate['convtime_'+str(j_species)] = np.mean(np.transpose(mol_values[top:bottom,:]),axis=1)
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.mean(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
	
						elif thinSize == 'all':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values),axis=1)
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values),axis=1), delimiter=",")
						
					print('Storing Catalyst Zone Gas and Surface Rate Data')
					for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
						
						rateTest = eval(rateStrings[j_species])
						tempSurface = pd.DataFrame(rateTest)
						tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)						
						#sys.exit()
						#timeBetween = time.time()

						#(np.asarray(parameters))
						#newFile.to_csv(reaction_file,header=None,index=False)	

						#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")			
						#print(time.time()-timeBetween)
						
				else:
					pulse_path = './'+reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(k_pulse+1)+'/'
					generateFolder(pulse_path)
					cat_dataRate = {}
					
					for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
						new_values = [] 
						mol_values = []
		
						for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
	
							if z != 0:
								#value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
								new_values = np.vstack((mol_values,cat_data['rawData'][:,z*(all_molecules)+j_species]))
								mol_values = new_values
					
							else:					
								new_values = cat_data['rawData'][:,z] 
								mol_values = cat_data['rawData'][:,z]
						
						top = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+1
						bottom = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-1				
	
						## Change so that only thin-zone data is being stored
						if thinSize == 'cat':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values[top:bottom,:]),axis=1)
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
				
						elif thinSize == 'point':
							centerPoint = int((top + bottom)/2)
							cat_dataRate['convtime_'+str(j_species)] = np.transpose(mol_values[centerPoint,:])
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.transpose(mol_values[centerPoint,:]), delimiter=",")
				
						elif thinSize == 'average':
							cat_dataRate['convtime_'+str(j_species)] = np.mean(np.transpose(mol_values[top:bottom,:]),axis=1)
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.mean(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
	
						elif thinSize == 'all':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values),axis=1)
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values),axis=1), delimiter=",")
	
					print('Storing Catalyst Zone Gas and Surface Rate Data')
					for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
						
						rateTest = eval(rateStrings[j_species])
						tempSurface = pd.DataFrame(rateTest)
						tempSurface.to_csv(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
						
						#np.savetxt(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")		
	
	
			#if reac_input['Thin-Zone Analysis'].lower() == 'true':
			#	if int(reac_input['Number of Pulses']) == 1:
			#		for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
			#			
			#			np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
			#			rateTest = eval(rateStrings[j_species])
			#			np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")
			#	else:
			#		pulse_path = './'+reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(k_pulse+1)+'/'
			#		generateFolder(pulse_path)
			#		for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
			#			np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
			#			rateTest = eval(rateStrings[j_species])
			#			np.savetxt(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")
			#	np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/Inert.csv', np.array(cat_data['conVtime_'+str(all_molecules-int(reac_input['Number of Inerts']))]), delimiter=",")
	
			#############################################################
			################# OPTIMIZATION EXECUTION ####################
			#############################################################
	
			fitting_time = time.time()
			
			if (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or fit_temperature == True):
	
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
					low_bounds.append(0)
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
					problem = MinimizationProblem(rf_2,bounds = (low_bounds,up_bounds))
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
				for j_species in range(0,monitored_gas+int(reac_input['Number of Inerts'])):
					dictionary_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
					new_data = np.asarray(graph_data['timing'])
					dictionary_of_numpy_data[legend_label[j_species]] = np.vstack((dictionary_of_numpy_data[legend_label[j_species]],new_data))
				
				for k_species in range(monitored_gas+1,len(necessary_values['reactants'])+1):
					dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.empty(shape=[0, len(graph_data['timing'])])
					new_data = np.asarray(graph_data['timing'])
					dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.vstack((dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]],new_data))
	
			#############################################################
			################### PLOT OUTLET FLUX DATA ###################
			#############################################################
	
			if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Display Experimental Data'].lower() == 'true':
	
				if reac_input['Display Objective Points'].lower() == 'true' or reac_input['Fit Parameters'].lower() == 'true':
					if k_pulse == 0:
						for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
							if objSpecies[k_fitting] == '1':
								ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting])
			
				if reac_input['Display Objective Points'].lower() == 'true':# and reac_input['Fit Inert'].lower() == 'true':
					#if k_pulse > 0:
					for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
					#for k_fitting in range( int(len(legend_label)-reac_input['Number of Inerts']) ,len(legend_label[int(len(legend_label)+reac_input['Number of Inerts'])])):
						if objSpecies[k_fitting] == '1':
							ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting], alpha=0.3)
	
	
				for k,j in enumerate(graph_data):
					if j != 'timing':
	
						if k_pulse > 0:
							pass
						else:
							dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/flux_data/'+legend_label[k]+'.csv',header=None)
							ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
					else:
						pass
	
			#############################################################
			############### PLOT ANALYTICAL INERT SOLUTION ##############
			#############################################################
	
			if reac_input['Infinite Inert'].lower() == 'true':
				if k_pulse == 0:
	
					# For reactant / product species
					analyticalTiming = np.arange(0, reac_input['Time Steps']*dt, dt).tolist()
					for kjc in range(0,monitored_gas):
						outlet = []
						
						if reac_input['Scale Output'].lower() == 'true':
							factor = 1
						else:
							factor = float(species_pulse_list[kjc])*reac_input['Reference Pulse Size']
						
						for k in range(0,int(reac_input['Time Steps'])+1):
							analyticalValue = 0
							if k*dt - species_time[kjc] > 0:
								for n in range(0,50):
									analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[kjc])*(D[kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
							else: 
								analyticalValue = -1
							if analyticalValue < 0:
								outlet.append(0)
							else:
								outlet.append(factor*(D[kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
						
						# Store analytical solution data
						np.savetxt('./analyticalCO19000.csv', outlet, delimiter=",")
	
						try:
							ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
							#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
						except ValueError:
							outlet = outlet[:-1]
							ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
							#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
					
					# For inert species
					for kjc in range(0,int(reac_input['Number of Inerts'])):
						outlet = []
						outlet.append(0)
					
						zAnalytical = 1
	
						while zAnalytical*dt < species_time[monitored_gas+kjc]:
							outlet.append(0)
							zAnalytical+=1	
						outlet.append(0)
	
						if reac_input['Scale Output'].lower() == 'true':
							factor = 1
						else:
							factor = float(species_pulse_list[monitored_gas+kjc])*reac_input['Reference Pulse Size']
	
						for k in range(zAnalytical+1,int(reac_input['Time Steps'])+1):
							analyticalValue = 0
							for n in range(0,50):
								analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[monitored_gas+kjc])*(D[monitored_gas+kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
							outlet.append(factor*(D[monitored_gas+kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
	
						try:
							ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
							#ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
						except ValueError:
							outlet = outlet[:-1]
							ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
							#ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
					
			#############################################################
			################## ADDITION OF WHITE NOISE ##################
			#############################################################
	
			sig = 0.1
			beta_2 = 0.00270
			w_2 = 2*3.14159*70
	
			for k,j in enumerate(graph_data):
				if j != 'timing':
					if reac_input['Noise'].lower() == 'true':
						for z in range(0,int(reac_input['Time Steps'])):
							graph_data[j][z] += np.random.normal(0,1)*sig +beta_2*np.cos(w_2*(k*dt))
					if k_pulse > 0:
						ax2.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
					else:
						ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)
					new_data = np.asarray(graph_data[j])
					dictionary_of_numpy_data[legend_label[k]] = np.vstack((dictionary_of_numpy_data[legend_label[k]],new_data))
				else:
					pass	
	
			name_list = necessary_values['reactants'][monitored_gas:]
	
			#############################################################
			################## SAVE OUTLET FLUX DATA ####################
			#############################################################
			
			for j_species in range(0,monitored_gas+int(reac_input['Number of Inerts'])):
				tempDict = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
				np.savetxt('./'+reac_input['Output Folder Name']+'_folder/flux_data/'+legend_label[j_species]+'.csv', tempDict, delimiter=",")
	
		ax2.legend(title="Gas Species")
		plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		
		if reac_input['Store Graph'].lower() == 'true':
			plt.savefig('./'+reac_input['Output Folder Name']+'_folder/graphs/flux_data.png')
		
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
							#alter.iloc[51,1] = 'FALSE'
					
					alter.to_csv(input_file,header=None,index=False)	
				
					try:
						print()
						print('Iteration: '+str(k_num+1))
						call_sim()
						sys.exit()
					except:
						k_num = things
			
				generateGif(legend_label[:len(legend_label)], reac_input['Experimental Data Folder']+'/flux_data', './'+reac_input['Output Folder Name']+'_folder/fitting', len(constants), constants, reactor_kinetics_input['reactions_test'], times)
			
				for k_num in range(0,things):
					shutil.rmtree('./'+reac_input['Output Folder Name']+'_folder/fitting/iter_'+str(k_num)+'_folder') 
	
			user_data = pd.read_csv('./'+reac_input['Output Folder Name']+'_folder/'+input_file,header=None)
			user_data.to_csv(input_file,header=None,index=False)
	
	
		if sampling == True:
			outValue = float(jfunc_2)
			tape2.clear_tape()		
			return outValue
	
		else:
			return graph_data, legend_label, necessary_values['reactants']
	
	#############################################################
	############## INPUTFILE DICTIONARY READING #################
	#############################################################
	
	def call_sim():
		
		reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm = inputForm)
		
		reactor_kinetics_input['Number of Pulses'] = pulseNumber
		reactor_kinetics_input['Fit Parameters'] = 'FALSE'
		reactor_kinetics_input['Display Experimental Data'] = 'FALSE'
		reactor_kinetics_input['Display Objective Points'] = 'FALSE'
		reactor_kinetics_input['Uncertainty Quantification'] = 'FALSE'
		reactor_kinetics_input['Sensitivity Analysis'] =  'TRUE'
		reactor_kinetics_input['Reactor Type'] = 'tap'
		reactor_kinetics_input['Knudsen Test'] = 'FALSE'
		reactor_kinetics_input['Noise'] = 'FALSE'
		reactor_kinetics_input['Fit Inert'] = 'FALSE'
		reactor_kinetics_input['Infinite Inert'] = 'FALSE'

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

		if fitInert != None:
			reactor_kinetics_input['Inert Fit'] =  'TRUE'

		sens_type = sensitivityType
		if (sens_type != 'trans' and sens_type != 'total') and reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true':
			print('loc 1')
			print('Error: sensitivity analysis must be "trans" or "total"')
			sys.exit()

		reactor_kinetics_input['Advection'] = 'FALSE'

		if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
			if sens_type == 'trans' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
				for parameters in kinetic_parameters:
					reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm=inputForm)
					
					reactor_kinetics_input['Sensitivity Parameter'] = parameters
					reactor_kinetics_input['Number of Pulses'] = pulseNumber
					reactor_kinetics_input['Reactor Type'] = 'tap'
					reactor_kinetics_input['Knudsen Test'] = 'FALSE'
					reactor_kinetics_input['Fit Inert'] = 'FALSE'
					reactor_kinetics_input['Fit Parameters'] = 'FALSE'
					reactor_kinetics_input['Display Experimental Data'] = 'FALSE'
					reactor_kinetics_input['Fitting Gif'] = 'FALSE'
					reactor_kinetics_input['Uncertainty Quantification'] = 'FALSE'
					reactor_kinetics_input['Infinite Inert'] = 'FALSE'
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

					if fitInert != None:
						reactor_kinetics_input['Inert Fit'] =  'TRUE'

					sens_type = sensitivityType					

					print(parameters)
					reactor_kinetics_input['Display Graph'] = 'FALSE'
					
					if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or reactor_kinetics_input['Fit Inert'].lower() == 'true':
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

def run_tapsolver(timeFunc,store_flux_func=True,includeNoise=False,pulseNumber = 1,store_thin_func=False,inputFile = './input_file.csv',input_form = 'new'):
	
	if store_thin_func == True:
		store_thin_func = 'TRUE'
	else:
		store_thin_func = 'FALSE'
	
	if store_flux_func == True:
		store_flux_func = 'TRUE'
	else:
		store_flux_func = 'FALSE'

	if includeNoise == True:
		includeNoise = 'TRUE'
	else:
		includeNoise = 'FALSE'


	general_run(timeFunc,store_flux_func='TRUE',store_thin_func=store_thin_func,pulseNumber = pulseNumber,input_file = inputFile,noise=includeNoise,inputForm = input_form)

def run_sensitivity(timeFunc,sensType=None,inputFile = './input_file.csv'):
	if sensType == None:
		sensType = 'total'
	general_run(timeFunc,store_flux_func='FALSE',store_thin_func='FALSE',input_file = inputFile,sensitivityType=sensType,input_form = 'old')

def fit_parameters(timeFunc,optim = 'L-BFGS-B',inputFile = './input_file.csv',input_form='new',inertFitting=None):
	if inertFitting == None:
		general_run(timeFunc,optimization=optim,input_file = inputFile,inputForm = input_form)
	else:
		general_run(timeFunc,optimization=optim,input_file = inputFile,fitInert=True,inputForm = input_form)

	sys.exit()

def run_uncertainty(timeFunc,inputFile = './input_file.csv'):
	general_run(timeFunc,uncertainty_quantificaiton=True,input_file = inputFile,inputForm=input_form)
	sys.exit()

def fitting_gif(timeFunc,outputName='./flux.png'):
	general_run(timeFunc,fitting_gif=True,input_file = './input_file.csv',inputForm=input_form)

def vary_Input(variableToChange, newValue, input_file='./input_file.csv'):
	df1 = pd.read_csv(input_file,header = None)
	#print(df1[0])
	#print(variableToChange)
	#print(df1[0].str.index(variableToChange))
	cellRow = df1[df1[0]==variableToChange].index.values[0]
	df1.iloc[cellRow,1] = newValue
	
	df1.to_csv('./input_file.csv',header=None,index=False)
	#print(df1[df1[0].str.match(variableToChange)])

	#to_csv(fileName)

def fluxGraph(input_file = './input_file.csv',pulse=None,dispExper=False,analytical=False,objectivePoints=False,displayGraph=True,storeGraph=False,outputName='./flux.png',inputForm='old'):
	timeFunc = 0.4

	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm = inputForm)
	#reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file)

	reac_input = reactor_kinetics_input
	reac_input['Pulse Duration'] = timeFunc
	reac_input['Infinite Inert'] = 'FALSE'

	if dispExper != False:
		reac_input['Display Experimental Data'] = 'TRUE'
		reac_input['Display Objective Points'] = 'TRUE'
		reac_input['Objective Points'] = 'all'

	if objectivePoints != False:
		reac_input['Display Objective Points'] = 'TRUE'
		reac_input['Objective Points'] = 'all'
#
	if analytical != False:
		reac_input['Infinite Inert'] = 'TRUE'

	reac_input['Advection'] = 'FALSE'
	reac_input['Time Steps'] = reac_input['Pulse Duration']*1000
	dt = 0.001
	reac_input['Reference Pulse Size'] = 1

	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,fit_temperature)

	reactor_kinetics_input['Reactor Type'] = 'tap'
	fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),'FALSE')

	monitored_gas = necessary_values['molecules_in_gas_phase']
	reac_input['Scale Output'] = 'FALSE'

	reac_input['Objective Species'] = '1,1,1,0'


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
		output_fitting = curveFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)


	r_param, dx_r, dx2_r, frac_length, cat_location = establishMesh(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['Mesh Size'])
	cat_location = 1 - reac_input['Catalyst Location']
	dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
	eb = np.array((reac_input['Void Fraction Inert'],reac_input['Void Fraction Catalyst'],reac_input['Void Fraction Inert']))
	ca = (reac_input['Reactor Radius']**2)*3.14159 
	point_volume = dx_r * ca * eb[0]
	Inert_pulse_conc = reac_input['Reference Pulse Size']/(point_volume)
		
	# Define the diffusion constants
	ref_rate = np.append(reac_input['Reference Diffusion Inert'],reac_input['Reference Diffusion Catalyst'])
	ref_rate = np.append(ref_rate,ref_rate[0])  	
	D = np.empty((len(reac_input['Mass List'].split(',')),3))
	
	def diff_func(ref_mass,ref_T,mol_mass,ref_r):
		return ref_r*(mp.sqrt(ref_mass*reac_input['Reactor Temperature'])/mp.sqrt(ref_T*mol_mass))
	Dout = []
	Din = []
	constantTemp = reac_input['Reactor Temperature']
	constantTemp = Constant(constantTemp)
	testTemp = reac_input['Reactor Temperature']
	testTemp = Constant(testTemp)

	if fit_temperature == True:
		controls = []
		new = Control(constantTemp)
		controls.append(new)

	compMass_list = reac_input['Mass List'].split(',')

	for k,j in enumerate(reac_input['Mass List'].split(',')):
		for k_2,j_2 in enumerate(ref_rate):
			D[k,k_2] = diff_func(reac_input['Reference Mass'],reac_input['Reference Temperature'],float(j),j_2) ###??? should this (np.sum(r_param)**2) be here? This puts it in dimensional form!
						
			if k_2 == 0:
				Dout.append(Constant(D[k,k_2]))
			if k_2 == 1:
				Din.append(Constant(D[k,k_2]))

	# Experimental Plot
#

	if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Display Experimental Data'].lower() == 'true':
#
		#if reac_input['Display Objective Points'].lower() == 'true':
		#	for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
		#		print(legend_label[k_fitting])
		#		if objSpecies[k_fitting] == '1':
		#			ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting])
#	
		#if reac_input['Display Objective Points'].lower() == 'true':# and reac_input['Fit Inert'].lower() == 'true':
		#	#if k_pulse > 0:
		#	for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
		#	#for k_fitting in range( int(len(legend_label)-reac_input['Number of Inerts']) ,len(legend_label[int(len(legend_label)+reac_input['Number of Inerts'])])):
		#		print(legend_label[k_fitting])
		#		if objSpecies[k_fitting] == '1':
		#			ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting], alpha=0.3)
#
		
		for k,j in enumerate(legend_label):
			if j != 'timing': #and j != 'conVtime_1' and j != 'conVtime_2' and j != 'conVtime_3' and j != 'conVtime_0' 
				dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/flux_data/'+legend_label[k]+'.csv',header=None)
				ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
			else:
				pass


	# Inert Plot
	if reac_input['Infinite Inert'].lower() == 'true':
		# For reactant / product species
		analyticalTiming = np.arange(0, reac_input['Time Steps']*dt, dt).tolist()
		for kjc in range(0,monitored_gas):
			outlet = []
			
			if reac_input['Scale Output'].lower() == 'true':
				factor = 1
			else:
				factor = float(species_pulse_list[kjc])*reac_input['Reference Pulse Size']
			
			for k in range(0,int(reac_input['Time Steps'])+1):
				analyticalValue = 0
				if k*dt - species_time[kjc] > 0:
					for n in range(0,50):
						analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[kjc])*(D[kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
				else: 
					analyticalValue = -1
				if analyticalValue < 0:
					outlet.append(0)
				else:
					outlet.append(factor*(D[kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
			
			# Store analytical solution data
			np.savetxt('./analyticalCO19000.csv', outlet, delimiter=",")
			try:
				ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
			except ValueError:
				outlet = outlet[:-1]
				ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
		
		# For inert species
		for kjc in range(0,int(reac_input['Number of Inerts'])):
			outlet = []
			outlet.append(0)
		
			zAnalytical = 1
			while zAnalytical*dt < species_time[monitored_gas+kjc]:
				outlet.append(0)
				zAnalytical+=1	
			outlet.append(0)
			if reac_input['Scale Output'].lower() == 'true':
				factor = 1
			else:
				factor = float(species_pulse_list[monitored_gas+kjc])*reac_input['Reference Pulse Size']
			for k in range(zAnalytical+1,int(reac_input['Time Steps'])+1):
				analyticalValue = 0
				for n in range(0,50):
					analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[monitored_gas+kjc])*(D[monitored_gas+kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
				outlet.append(factor*(D[monitored_gas+kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
			try:
				ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical '+legend_label[monitored_gas+kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
			except ValueError:
				outlet = outlet[:-1]
				ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical '+legend_label[monitored_gas+kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)

	for knum,k in enumerate(necessary_values['reactants']):
		if knum < necessary_values['molecules_in_gas_phase']:

			dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/'+k+'.csv',header=None)
			graphLimit = dfTemp[0][len(dfTemp[0]) - 1]
			
			if type(pulse) == int:
				ax2.plot(dfTemp[0],dfTemp[pulse],color=colors[knum], ls = '--',label=legend_label[knum], alpha=0.7)
			
			elif type(pulse) == list:
				legendInclude = True
				for j in pulse:
					if legendInclude == True:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label=legend_label[knum], alpha=0.7)
						legendInclude = False
					else:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label='_nolegend_', alpha=0.7)
			else:
				legendInclude = True
				for j in range(1, len(dfTemp.columns)):
					if legendInclude == True:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label=legend_label[knum], alpha=0.7)
						legendInclude = False
					else:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label='_nolegend_', alpha=0.7)
	ax2.set_xlim(0,graphLimit)
	for k in range(0,int(reac_input['Number of Inerts'])):
		
		dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/Inert-'+str(k+1)+'.csv',header=None)
		if type(pulse) == int:
			ax2.plot(dfTemp[0],dfTemp[pulse],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label=legend_label[necessary_values['molecules_in_gas_phase']+k], alpha=0.7)
			
		elif type(pulse) == list:
			legendInclude = True
			for j in pulse:
				if legendInclude == True:
					print(legend_label[necessary_values['molecules_in_gas_phase']+k-1])
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label=legend_label[necessary_values['molecules_in_gas_phase']+k], alpha=0.7)
					legendInclude = False
				else:
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label='_nolegend_', alpha=0.7)

		else:
			legendInclude = True
			for j in range(1, len(dfTemp.columns)):
				if legendInclude == True:
					print(legend_label[necessary_values['molecules_in_gas_phase']+k-1])
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label=legend_label[necessary_values['molecules_in_gas_phase']+k], alpha=0.7)
					legendInclude = False
				else:
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label='_nolegend_', alpha=0.7)

	ax2.legend(title="Gas Species",loc = 1)

	if storeGraph == True:
		plt.savefig(outputName)

	if displayGraph == True:
		plt.show()

def pulsingVisualization(input_file = './input_file.csv',fileName = './output.gif',inputForm='old'):
	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm=inputForm)
	reac_input = reactor_kinetics_input
	reac_input['Advection'] = 'FALSE'
	reac_input['Reactor Type'] = 'tap'
	reac_input['Scale Output'] = 'FALSE'
	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,fit_temperature)
	fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])

	dfColumns = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/'+necessary_values['reactants'][0]+'.csv',header=None)
	

	def plotPulse(pulseNumber):
		
		fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),'FALSE')

		for knum,k in enumerate(necessary_values['reactants']):
			if knum < necessary_values['molecules_in_gas_phase']:
				dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/'+k+'.csv',header=None)
				ax2.plot(dfTemp[0],dfTemp[pulseNumber],color=colors[knum], ls = '--',label=legend_label[knum], alpha=0.7)

		for k in range(1,int(reac_input['Number of Inerts'])+1):
			dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/Inert-'+str(k)+'.csv',header=None)
			ax2.plot(dfTemp[0],dfTemp[pulseNumber],color=colors[necessary_values['molecules_in_gas_phase']+k-1], ls = '--', label=legend_label[k], alpha=0.7)

		ax2.legend(title="Gas Species",loc = 1)

		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax2.text(0.4, 1.1, 'Pulse: '+str(pulseNumber), fontsize=12,verticalalignment='top', bbox=props)

		fig2.canvas.draw()
		image = np.frombuffer(fig2.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig2.canvas.get_width_height()[::-1] + (3,))
		
		return image


	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	
	imageio.mimsave(fileName, [plotPulse(i) for i in list(range(1,len(dfColumns.columns)))], fps=4)

def concnDistGif(input_file = './input_file.csv',dataType='surf',fileName='./cat.gif',pulse=1,inputForm='old'):
	
	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm = inputForm)
	reac_input = reactor_kinetics_input
	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,fit_temperature)
	#fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])

	try:
		dataLocation = reac_input['Output Folder Name']+'_folder/thin_data/'
		dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)
	except:
		dataLocation = reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(pulse)+'/'
		dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)

	dataDictionary = {}

	for knum,k in enumerate(necessary_values['reactants']):
		dataDictionary[k] = pd.read_csv(dataLocation+k+'.csv',header=None)
	
	def plotCat(time,type_of_data):

		fig,ax = plt.subplots()
		for knum,k in enumerate(necessary_values['reactants']):
			#dfTemp = pd.read_csv(dataLocation+k+'.csv',header=None)
			if dataType == 'surf':
				if knum > necessary_values['molecules_in_gas_phase']:
					ax.plot(dataDictionary[k].iloc[:,time],label = k)
			elif dataType == 'gas':
				if knum < necessary_values['molecules_in_gas_phase']:
					ax.plot(dataDictionary[k].iloc[:,time],label = k)
			else:
				if necessary_values['reactants'] in dataType:
					ax.plot(dataDictionary[k].iloc[:,time],label = k)

		ax.legend('Species')
		fig.canvas.draw()
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
		return image
	print(dfColumns[0].count())
	
	imageio.mimsave(fileName, [plotCat(i,dataType) for i in list(range(0,dfColumns[0].count()))], fps=4)

def concDistPlot(input_file = './input_file.csv',dataType='surf',fileName='./cat.gif',pulse=1,input_form='old'):
	
	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm = input_form)
	reac_input = reactor_kinetics_input
	fit_temperature = False
	reac_input['Advection'] = 'FALSE'
	reac_input['Reactor Type'] = 'tap'
	reac_input['Scale Output'] = 'FALSE'
	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,fit_temperature)
	#fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])

	reac_input['Advection'] = 'FALSE'
	reac_input['Reactor Type'] = 'tap'
	reac_input['Scale Output'] = 'FALSE'
	fit_temperature = False

	try:
		dataLocation = reac_input['Output Folder Name']+'_folder/thin_data/'
		dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)
	except:
		dataLocation = reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(pulse)+'/'
		dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)

	dataDictionary = {}

	for knum,k in enumerate(necessary_values['reactants']):
		dataDictionary[k] = pd.read_csv(dataLocation+k+'.csv',header=None)

	dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)

	fig,ax = plt.subplots()
	for knum,k in enumerate(necessary_values['reactants']):
		if dataType == 'surf':
			if knum > necessary_values['molecules_in_gas_phase']:
				ax.plot(dataDictionary[k].iloc[int(len(dfColumns.columns)/2),:],label = k)
		elif dataType == 'gas':
			if knum < necessary_values['molecules_in_gas_phase']:
				ax.plot(dataDictionary[k].iloc[dataDictionary[k].iloc[:,int(len(dfColumns.columns)/2)],:],label = k)
		else:
			if necessary_values['reactants'] in dataType:
				ax.plot(dataDictionary[k].iloc[dataDictionary[k].iloc[:,int(len(dfColumns.columns)/2)],:],label = k)

	plt.show()

def input_construction(reactor_name='./reactor_definition.csv',reaction_name='./reaction_definition.csv',new_name='./input_file.csv'):
	
	start = pd.read_csv(reactor_name,header=None)
	startRows = len(start)
	startCols = len(start.columns)

	reaction_info = pd.read_csv(reaction_name,header=None)	
	reactionRows = len(reaction_info)
	reactionCols = len(reaction_info.columns)

	reactions = reaction_info.iloc[:,0].tolist()
	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions)

	numberOfGasses = 0
	for j in reactants:
		if j.find('*') < 0:
			numberOfGasses+= 1
		else: 
			break

	gasses = reactants[:numberOfGasses]
	surface = reactants[numberOfGasses:]

	if len(gasses) > len(surface):
		N_cols = len(gasses)+1
	else:
		N_cols = len(surface)+1

	if N_cols < startCols:
		N_cols = startCols

	N_rows = 7

	d = pd.DataFrame(np.zeros((N_rows, N_cols)))

	for j_num,j in enumerate(gasses):
		d.iloc[0,j_num+1] = j
		d.iloc[1,j_num+1] = 1
		d.iloc[3,j_num+1] = molecularProperties(gasses[j_num],'mass')

	for j in range(len(gasses)+1,N_cols):
		d.iloc[0,j] = ''
		d.iloc[1,j] = ''
		d.iloc[2,j] = ''
		d.iloc[3,j] = ''

	for j in range(0,N_cols):
		d.iloc[4,j] = ''

	for j_num,j in enumerate(surface):
		d.iloc[5,j_num+1] = j
		if j_num+1 == len(surface):
			d.iloc[6,j_num+1] = 100
	for j in range(len(surface)+1,N_cols):
		d.iloc[5,j] = ''
		d.iloc[6,j] = ''
	
	d.iloc[0,0] = ''
	d.iloc[1,0] = 'Intensity'
	d.iloc[2,0] = 'Time'
	d.iloc[3,0] = 'Mass'
	d.iloc[4,0] = ''
	d.iloc[5,0] = ''
	d.iloc[6,0] = 'Initial Concentration'

	df = pd.read_csv(reactor_name,header=None)

	storageArray = np.zeros((N_rows+startRows+2+reactionRows+2, N_cols))
	newOutputFile = pd.DataFrame(storageArray)
	newOutputFile.iloc[:,:] = ''

	for j in range(0,len(start)):
		for k in range(0,len(start.columns)):
				newOutputFile.iloc[j,k] = start.iloc[j,k]

	for k in range(0,len(start.columns)):
		newOutputFile.iloc[len(start),k] = ''

	newOutputFile.iloc[len(start)+1,0] = 'Feed_&_Surface_Composition'

	for k in range(1,len(start.columns)):
		newOutputFile.iloc[len(start)+2,k] = ''

	for j in range(0,len(d)):
		for k in range(0,len(d.columns)):
			newOutputFile.iloc[startRows+j+2,k] = d.iloc[j,k]

	newOutputFile.iloc[len(start)+2+len(d)+1,0] = 'Reaction_Information'

	for j in range(0,len(reaction_info)):
		for k in range(0,len(reaction_info.columns)):
			newOutputFile.iloc[len(start)+2+len(d)+j+2,k] = reaction_info.iloc[j,k]

	newOutputFile.to_csv(new_name,header=None,index=False)

def exampleReaction(reaction_file = './reaction_example.csv'):

	parameters = ['CO + * <-> CO*',1,10]
	newFile = pd.DataFrame(np.asarray(parameters))
	newFile.to_csv(reaction_file,header=None,index=False)	
	
def define_reactor(reactor_name='./reactor_definition.csv',transport_type=['Knudsen','Advection']):
	
	parameters = ['Reactor Radius','Reactor Temperature','Mesh Size','Catalyst Mesh Density','Output Folder Name','Experimental Data Folder']
	initial_values = [0.2,700,400,2,'results','None']
	#parameters = ['Reactor_Information','Reactor Length','Mesh Size','Catalyst Fraction','Catalyst Mesh Density','Reactor Radius','Catalyst Location','Void Fraction Inert','Void Fraction Catalyst','Reactor Temperature','Reference Diffusion Inert','Reference Diffusion Catalyst','Reference Temperature','Reference Mass','Output Folder Name','Experimental Data Folder','Advection','','Feed_&_Surface_Composition']
	#initial_values = ['',3.832,400,0.02,5,0.2,0.5,0.4,0.4,700,13.5,13.5,700,40,'results','./probeTruncated/Mn',0,'','']


	if 'Knudsen' in transport_type:
		knudsenList = ['Reference Diffusion Inert','Reference Diffusion Catalyst','Reference Temperature','Reference Mass']
		knudsenValues = [13.5,13.5,700,40]
		for jnum,j in enumerate(knudsenList):
			parameters.append(j)
			initial_values.append(knudsenValues[jnum])

	if 'Advection' in transport_type:
		advectionList = ['Advection Value']
		advectionValues = [0.0]
		for jnum,j in enumerate(advectionList):
			parameters.append(j)
			initial_values.append(0)

	d = pd.DataFrame(np.zeros((len(parameters)+3, 4)))
	
	d.iloc[0,0] = 'Reactor_Information'
	d.iloc[0,1] = ''
	d.iloc[0,2] = ''
	d.iloc[0,3] = ''

	d.iloc[1,0] = 'Zone Length'
	d.iloc[1,1] = 2.95
	d.iloc[1,2] = 0.15
	d.iloc[1,3] = 2.95

	d.iloc[2,0] = 'Zone Void'
	d.iloc[2,1] = 0.4
	d.iloc[2,2] = 0.4
	d.iloc[2,3] = 0.4

	for j in range(0,len(parameters)):
		d.iloc[j+3,0] = parameters[j]
		d.iloc[j+3,1] = initial_values[j]
		d.iloc[j+3,2] = ''
		d.iloc[j+3,3] = ''

	d.to_csv(reactor_name,header=None,index=False)

def directoryGenerator(path_name):
	try:  
		os.mkdir(path_name)
	except OSError:  
		pass
	else:  
		pass

def bulk_structure(reactor_setup,iterDict,baseName=None):

	# generate csv submission calculation

	directoryGenerator('./'+baseName)
	copyfile(reactor_setup,'./'+baseName+'/'+reactor_setup)

	if baseName == None:
		print('No variables are changed. No generation performed.')
		sys.exit()

	#headers = []
	headers = ['directory']

	generalVariables = list(iterDict.keys())

	for j in iterDict.keys():
		if j == 'mechanisms':
			headers.append(j)
		else:
			for k in iterDict[j].keys():
				headers.append(k+'_'+j)

	def iterativeGeneration(newNpArray,directoryName,currentList,currentkey,currentSubKey):

		if list(iterDict.keys())[currentkey] == "mechanisms":
		
			for jnum,j in enumerate(list(iterDict[list(iterDict.keys())[currentkey]].keys())):
			
				newdirectoryName = directoryName+str(jnum)
				newList = [j]
			
				copyfile(iterDict["mechanisms"][j],'./'+baseName+'/'+iterDict["mechanisms"][j])
				# Construct new base input_file
				input_construction(reactor_name='./'+reactor_setup,reaction_name=iterDict["mechanisms"][j],new_name='./input_file.csv')				
				newNpArray = iterativeGeneration(newNpArray,newdirectoryName,newList,currentkey+1,0)

		else:
			for jnum,j in enumerate(iterDict[list(iterDict.keys())[currentkey]][list(iterDict[list(iterDict.keys())[currentkey]].keys())[currentSubKey]]):
				newdirectoryName = directoryName+str(jnum)
				newList = currentList[:]
				newList.append(j)
						
				tempInput = pd.read_csv('./input_file.csv',header=None)
				currentSpecies = list(iterDict[list(iterDict.keys())[currentkey]])[currentSubKey]
				if list(iterDict.keys())[currentkey] == "surface":
					currentSurface = list(tempInput.iloc[21,1:])
					tempInput.iloc[22,1+currentSurface.index(currentSpecies)] = float(j)
				elif list(iterDict.keys())[currentkey] == "intensity":
					currentGasses = list(tempInput.iloc[16,1:])
					tempInput.iloc[17,1+currentGasses.index(currentSpecies)] = float(j)
				elif list(iterDict.keys())[currentkey] == "time":
					currentGasses = list(tempInput.iloc[16,1:])
					tempInput.iloc[18,1+currentGasses.index(currentSpecies)] = float(j)
					tempInput.to_csv('./input_file.csv',header=None,index=False)

				if currentSubKey+1 != len(list(iterDict[list(iterDict.keys())[currentkey]].keys())):
					newNpArray = iterativeGeneration(newNpArray,newdirectoryName,newList,currentkey,currentSubKey+1)
				elif currentkey+1 != len(list(iterDict.keys())):
					newNpArray = iterativeGeneration(newNpArray,newdirectoryName,newList,currentkey+1,0)
				else:
					# Alter Parameters
					# Move new files to subdirectory
					# input_file, tap_sim, func_sim, vari_form, reac_odes, bash
					directoryNameList = [newdirectoryName]
					newList = directoryNameList+newList
					newNpArray = np.vstack((newNpArray,newList))
					directoryGenerator(newdirectoryName)
					tempInput.to_csv('./input_file.csv',header=None,index=False)
					copyfile('./tap_sim.py',newdirectoryName+'/tap_sim.py')
					copyfile('./run.sh',newdirectoryName+'/run.sh')
					copyfile('./func_sim.py',newdirectoryName+'/func_sim.py')
					copyfile('./reac_odes.py',newdirectoryName+'/reac_odes.py')
					copyfile('./vari_form.py',newdirectoryName+'/vari_form.py')
					copyfile('./input_file.csv',newdirectoryName+'/input_file.csv')
				
		return newNpArray
						
	dataStructure = np.asarray(headers)
	dataValues = np.asarray(headers[2:])

	dataStructurenew = iterativeGeneration(dataStructure,'./'+baseName+'/',[],0,0)
	df = pd.DataFrame(dataStructurenew)
	df.to_csv('./'+baseName+'/'+'directoryList.csv',header=None,index=False)

	valueList = []

	for j in list(iterDict.keys())[1:]:
		#print(j)
		for k in iterDict[j].keys():
			valueList.append(iterDict[j][k])

	dataValues = np.vstack((list(dataValues),valueList))

	df2 = pd.DataFrame(dataValues)

	df2.to_csv('./'+baseName+'/'+'changedVariables.csv',header=None,index=False)


def gen_help(parameter):
	sys.exit()

def folder_construcntion():
	sys.exit()

def bulk_analysis():
	sys.exit()
