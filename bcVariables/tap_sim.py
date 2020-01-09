from fenics import *
from fenics_adjoint import *
from func_sim import *
from vari_form import *
from reac_odes import *
import mpmath
import matplotlib.pyplot as plt
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
import ufl
import scipy
import pip
import pkg_resources
#from hippylib import *
import warnings

sampling = False
if sampling == True:
	from muq import pymuqModeling as mm # Needed for Gaussian distribution
	from muq import pymuqApproximation as ma # Needed for Gaussian processes
	from muq import pymuqSamplingAlgorithms as ms # Needed for MCMC
	from muq import pymuqUtilities as mu

import faulthandler


# New values for the objective function

objectiveAnalysis = False
selectivity = False
conversion = False
yieldValue = True

### How much concentration data do you want stored for 'thin-zone analysis'
# 'point' = only center point in the catalyst zone
# 'average' = average concentration values along the catalyst zone
# 'cat' = all points in the catalyst zone
# 'all' = all points in the reactor

thinSize = 'point'

### Sensitivity Type
# 'trans' = transient sensitivity analysis
# 'total' = sensitivity of summed objective function

sens_type = 'total'

warnings.simplefilter(action='ignore', category=FutureWarning)

#############################################################
################ TAPsolver FEniCS Function ##################
#############################################################

def tap_simulation_function(reactor_kinetics_input,constants_input,fitting_input):

	# Define or clear working Tape of FEniCS / Dolfin-Adjoint
	tape2 = Tape()
	tape2.clear_tape()
	set_working_tape(tape2)
	
	#############################################################
	#### GENERATE OUTPUT FOLDERS (depending on input values) ####
	#############################################################
	
	kVals = constants_input.copy()
	reac_input = reactor_kinetics_input

	reac_input['Optimization Method'] = 'BFGS'
	reac_input['Objective Points'] = 'all'

	if reac_input['Advection Value'] > 0.0:
		reac_input['Advection'] = 'true'
	else:
		reac_input['Advection'] = 'false'

	path = './'+reac_input['Output Folder Name']+'_folder/'
	generateFolder(path)
		
	path_3 = './'+reac_input['Output Folder Name']+'_folder/flux_data/'
	generateFolder(path_3)
	
	if reac_input['Thin-Zone Analysis'].lower() == 'true':
		path_4 = './'+reac_input['Output Folder Name']+'_folder/thin_data/'
		generateFolder(path_4)

	path_5 = './'+reac_input['Output Folder Name']+'_folder/graphs/'
	generateFolder(path_5)

	if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true':
		path_6 = './'+reac_input['Output Folder Name']+'_folder/fitting/'
		generateFolder(path_6)

	# Declare and define the constants of interest
	r_const = constants_input
	r_fit = fitting_input

	for j in r_const:
		r_const[j] = Constant(r_const[j])

	if reac_input['Fit Parameters'].lower() == 'true' or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true'):
		
		controls = []
		legend_2 = []
		for j in r_fit:
		#for j in r_const:
			controls.append(Control(r_const[j]))
			legend_2.append(j)

	# Store input file in output folder	
	user_data = pd.read_csv('./input_file.csv',header=None)
	user_data.to_csv(path+'input_file.csv',header=None,index=False)
	original_input_structure = user_data.copy()
	
	#############################################################
	#### READ AND DEFINE INITIAL PARAMETERS FROM INPUT FILE #####
	#############################################################

	set_log_level(30)
	tol = 1E-14

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
	for k,j in enumerate(reac_input['Mass List'].split(',')):
		for k_2,j_2 in enumerate(ref_rate):
			D[k,k_2] = Constant(diff_func(reac_input['Reference Mass'],reac_input['Reference Temperature'],float(j),j_2)) ###??? should this (np.sum(r_param)**2) be here? This puts it in dimensional form!
			testControl = Constant(D[k,k_2])
			if k_2 == 0:
				Dout.append(Constant(D[k,k_2]))
			if k_2 == 1:
				Din.append(Constant(D[k,k_2]))

	# Construct the rate expression based on the microkinetic model
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Active Sites'],reac_input['Number of Inerts'],reac_input['Advection'],False)

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
	u_1 = Function(V)
	u_2 = Function(V)
	u_3 = Function(V)
	
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
	
	bcs = defineBCs(reac_input['Reactor Type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['Number of Reactants'],all_molecules,reac_input['Pulse Ratio'],boundary_L,boundary_R,reac_input['Number of Inerts'])
	
	fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])
	ax2.set_xlim(0,reac_input['Pulse Duration'])

	# Define controls for optimization or sensitivity analysis
	if reac_input['Fit Inert'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
			controls.append(Control(Dout[j]))
			legend_2.append(j)

	if str(reac_input['Objective Species']).find(',') != -1:
		objSpecies = list(reac_input['Objective Species'].split(','))
	else:
		objSpecies = [str(int(reac_input['Objective Species']))]

	t = Constant(0)
	uv = Expression('(1/sqrt(3.14159))*exp(-x[0]*100*t)',t=0,degree=1)
	#uv = Expression('1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t', {'alpha': alpha, 'beta': beta})

	#############################################################
	######## EVALUATE AND INITIALIZE PDE FOR FEniCS #########
	#############################################################
	
	try:
		theta = 1
		Ftemp = eval(necessary_values['F'])
		theta = 0.5
		F = eval(necessary_values['F'])#-eval(uv)
	except NameError:
		errorOutput(reac_input['reactions_test'])
	
	#############################################################
	##### DEFINE METHOD OF OPTIMIZATION (OBJECTIVE FUNC.) #######
	#############################################################
	
	if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true') or sampling == True or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true'):
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

	if reac_input['Fit Inert'].lower() == 'true':
		try:
			if type(reac_input['Objective Points']) == int:
				output_fitting = pointFitting(legend_label[int(reac_input['Number of Inerts']):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'])
			elif reac_input['Objective Points'] == 'all':
				output_fitting = curveFitting(legend_label[int(reac_input['Number of Inerts']):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'])
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
		
	if reac_input['RRM Analysis'].lower() == 'true':
		path_2 = reac_input['Output Folder Name']+'_folder/RRM_analysis/'
		generateFolder(path_2)
		
		path_molecules = path_2+reac_input['Sensitivity Parameter']
		generateFolder(path_molecules)
		rrmFolder = path_molecules

		generateFolder(path_molecules+'/pointValue')
		generateFolder(path_molecules+'/thinValue')

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
		print('Sensitivity analysis is not properly defined.')
		sys.exit()


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
	
	if '/' in str(reac_input['Pulse Ratio']):
		list_of_feed = reac_input['Pulse Ratio'].split('/')
		list_of_time = reac_input['Pulse Time'].split('/')
		pulse_variation = len(list_of_feed)
		pulse_time = len(list_of_time)
		reactant_feed = []
		reactant_time = []
		for k in range(0,len(reac_input['Pulse Ratio'].split('/'))):
			reactant_feed.append(list_of_feed[k].split(','))
			reactant_time.append(list_of_time[k].split(','))
	else:
		pulse_variation = 1
		reactant_feed = []
		reactant_time = []
		try:
			reactant_feed.append(reac_input['Pulse Ratio'].split(','))
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
		except AttributeError:
			reactant_feed.append(reac_input['Pulse Ratio'])
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))

	#############################################################
	####### STORE DATA BASED ON USER SPECIFICATION ##############
	#############################################################

	current_reac_set = 1
	tot_objective = 0

	if reac_input['Thin-Zone Analysis'].lower() == 'true':
		rateStrings = rateEqs(rate_array,rev_irr)
	
	if reac_input['RRM Analysis'].lower() == 'true':
		rrmStringsThin = rrmEqs(rate_array,rev_irr,'dx(1)')	
		rrmStringsPoint = rrmEqs(rate_array,rev_irr,'dT(1)')

	if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or sampling == True or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true'):
		for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
			if objSpecies[k_fitting] == '1':

				for timeStep in range(0,len(output_fitting[legend_label[k_fitting]]['times'])):
					output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)

	if reac_input['Fit Inert'].lower() == 'true':
		for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
			for timeStep in range(0,len(output_fitting[legend_label[-1]]['times'])):
				output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)

	if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':

		if sens_type == 'trans' or reac_input['Uncertainty Quantification'].lower() == 'true':
			c = r_const[reac_input['Sensitivity Parameter']]
			c.tlm_value = r_const[reac_input['Sensitivity Parameter']]
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


	if reac_input['RRM Analysis'].lower() == 'true':
		thinSensFunc = {}
		pointSensFunc = {}
		
		thinSensFuncRate = {}
		pointSensFuncRate = {}

		for k_gasses in range(0,len(necessary_values['reactants'])):
			thinSensFunc[str(k_gasses)] = []
			pointSensFunc[str(k_gasses)] = []
			thinSensFuncRate[str(k_gasses)] = []
			pointSensFuncRate[str(k_gasses)] = []

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

		species_pulse_list = reactant_feed[current_reac_set-1]
		
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
			print('Sensitivity analysis is not properly defined.')
			sys.exit()

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
		if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or sampling == True:
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

			if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or sampling == True or (sens_type == 'total' and reac_input['Sensitivity Analysis'].lower() == 'true'):

				#print(legend_label)
				#sys.exit()

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
									if legend_label[k_fitting].find('Inert') == -1:
									#if legend_label[k_fitting] != 'Inert':
										jfunc_2 += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))							
									else:
										pass

								except UnboundLocalError:
									if legend_label[k_fitting].find('Inert') == -1:
									#if legend_label[k_fitting] != 'Inert':
										jfunc_2 = assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
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

								jfunc_2 = assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
							else:
								pass

			if round(t,6) not in species_time:

				# To improve accuracy / utilize crank-nicolson, small steps using backward euler must be utilized.
				try:
					if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
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
							#if round(t,6) == 0.001:
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

					if ',' in str(species_pulse_list):
						for k in range(0,int(reac_input['Number of Reactants'])):
							
							if species_time[k] == round(t,6):
								u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells)))+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]
								
					else:
						u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']-10)-1)-(all_molecules)] = float(species_pulse_list)*Inert_pulse_conc

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
							#cat_data['rawData'].append(u_n.vector().get_local())
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
					if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
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
						new_val = (( u.vector().get_local()[(all_molecules)+k]))
						u_graph_data['conVtime_'+str(k)].append((new_val))

					for kGasses in range(0,len(necessary_values['reactants'])):
						sensFuncs[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dP(1)))

				elif sens_type == 'total':
					pass

				else:
					print('Sensitivity analysis is not properly defined.')
					sys.exit()
					
			if reac_input['RRM Analysis'].lower() == 'true':
				for kGasses in range(0,len(necessary_values['reactants'])):
					thinSensFunc[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dx(1)))
					pointSensFunc[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dT(1)))

				for kGasses in range(0,len(necessary_values['reactants'])):
					thinSensFuncRate[str(kGasses)].append(eval(rrmStringsThin[kGasses]))
					pointSensFuncRate[str(kGasses)].append(eval(rrmStringsPoint[kGasses]))

			sens_time_list.append(t)

			progressBar(t, reac_input['Pulse Duration'])
			u_n.assign(u)
			t += dt

		print()
		print(processTime(start_time))

		if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true':
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
					f.write("Change: "+str(djv))
					f.close

			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit()


		if reac_input['Sensitivity Analysis'].lower() == 'true':


			if sens_type == 'trans':
				for numEachSens,eachSens in enumerate(u_graph_data):
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


		if reac_input['RRM Analysis'].lower() == 'true':

			for numEachSens,eachSens in enumerate(thinSensFunc):
				newList = []
				for kSensNum, kSens in enumerate(thinSensFunc[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)
				np.savetxt(rrmFolder+'/thinValue'+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")

			for numEachSens,eachSens in enumerate(pointSensFunc):
				newList = []
				for kSensNum, kSens in enumerate(pointSensFunc[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)

				np.savetxt(rrmFolder+'/pointValue'+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
			
			for numEachSens,eachSens in enumerate(thinSensFuncRate):
				newList = []
				for kSensNum, kSens in enumerate(thinSensFuncRate[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)

				np.savetxt(rrmFolder+'/thinValue'+'/dr_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
		
			for numEachSens,eachSens in enumerate(pointSensFuncRate):
				newList = []
				for kSensNum, kSens in enumerate(pointSensFuncRate[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)

				np.savetxt(rrmFolder+'/pointValue'+'/dr_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
			
		if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true':
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
			if int(reac_input['Number of Pulses']) == 1:
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
						np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
			
					elif thinSize == 'point':
						centerPoint = int((top + bottom)/2)
						cat_dataRate['convtime_'+str(j_species)] = np.transpose(mol_values[centerPoint,:])
						np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.transpose(mol_values[centerPoint,:]), delimiter=",")
			
					elif thinSize == 'average':
						cat_dataRate['convtime_'+str(j_species)] = np.mean(np.transpose(mol_values[top:bottom,:]),axis=1)
						np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.mean(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")

					elif thinSize == 'all':
						cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values),axis=1)
						np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values),axis=1), delimiter=",")


				for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					rateTest = eval(rateStrings[j_species])
					np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")			

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


				for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					rateTest = eval(rateStrings[j_species])
					np.savetxt(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")		


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
				
		if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true':

			start_time = time.time()
			print()
			print('Fitting Kinetic Parameters. Will take some time!')
			rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)

			low_bounds = []
			up_bounds = []

			for gt in range(0,len(controls)):
				low_bounds.append(0)
				up_bounds.append(np.inf)

			if reac_input['Optimization Method'] == 'L-BFGS-B' or reac_input['Optimization Method'] == '':
				u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds),tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
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
		
			if reac_input['Display Objective Points'].lower() == 'true' and reac_input['Fit Inert'].lower() == 'true':
				if k_pulse > 0:
					for k_fitting in range( int(len(legend_label)-reac_input['Number of Inerts']) ,len(legend_label[int(len(legend_label)+reac_input['Number of Inerts'])])):
						ax2.scatter(output_fitting[legend_label[monitored_gas+k_fitting]]['times'],output_fitting[legend_label[monitored_gas+k_fitting]]['values'],marker='^',color='r',label='Fitting'+legend_label[monitored_gas+k_fitting], alpha=0.3)


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

		sig = 0.05
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

			for k_num in range(0,things):
				alter = pd.read_csv('./input_file.csv',header=None)

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
						alter.iloc[51,1] = 'FALSE'
				
				alter.to_csv('./input_file.csv',header=None,index=False)	
			
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

		user_data = pd.read_csv('./'+reac_input['Output Folder Name']+'_folder/input_file.csv',header=None)
		user_data.to_csv('./input_file.csv',header=None,index=False)


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
	
	reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
	
	if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' or reactor_kinetics_input['RRM Analysis'].lower() == 'true' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
		
		if sens_type == 'trans':
			for parameters in kinetic_parameters:
				reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
				reactor_kinetics_input['Sensitivity Parameter'] = parameters
				reactor_kinetics_input['Display Graph'] = 'FALSE'
				if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or reactor_kinetics_input['Fit Inert'].lower() == 'true':
					print('')
					print('')
				
					print('Running the Sensitivity/RRM Analysis and Parameter Fitting methods simultaniously is not possible due to conflicts between the tangent linear and adjoint methods.')
					print('Please run again with one of these methods excluded.')
					sys.exit()
				if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true':
					print('')
					print('Processing '+parameters)
				graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)
		
		elif sens_type == 'total':
			graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)
	else:
		graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)


#############################################################
######## SAMPLING METHOD UNCERTAINTY PROPAGATION ############
#############################################################

if sampling == True:
	class fenicsModel(mm.PyModPiece):
		def __init__(self):
		#def __init__(self,a,b,c,d,e,f):
			mm.PyModPiece.__init__(self, [6], [1]) # One input containing 2 components [2]) # One output containing 2 components
			#self.a = a
			#self.b = b
			#self.c = c
			#self.d = d
			#self.e = e
			#self.f = f

		def EvaluateImpl(self, inputs):
			reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
			inputs = inputs[0].tolist()
			for j_num,j in enumerate(kinetic_parameters):
				kinetic_parameters[j] = inputs[j_num]

			print(kinetic_parameters)
			outputValue = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)
			
			#z = inputs[0]
			#m = np.zeros((6))

			print('done')
			return outputValue

	reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
	
	graph = mm.WorkGraph()

	zDist = mm.Gaussian(np.zeros((1)))

	#invf = fenicsModel(1,2,3,4,5,6)
	invf = fenicsModel()

	graph.AddNode(zDist.AsDensity(), "Gaussian Reference")
	graph.AddNode(invf, "TAP Evaluation")

	graph.AddEdge("TAP Evaluation", 0, "Gaussian Reference", 0)

	tgtDens = graph.CreateModPiece("Gaussian Reference")

	
	numSamps = 200
	paramDim = 6
	
	mSamps = np.zeros((paramDim, numSamps))

	samplingVersion = 'mcmc'

	#############################################################
	################ REJECTION SAMPLING METHOD ##################
	#############################################################

	if samplingVersion == 'rejection':
		propMu = np.array([29,0,3.6,11,8,0])
		propCov = np.array([ [2, 0, 0, 0, 0, 0],[ 0,0.25, 0, 0, 0, 0],[ 0, 0,1, 0, 0, 0],[ 0, 0, 0,2, 0, 0],[ 0, 0, 0, 0,3, 0],[ 0, 0, 0, 0, 0,0.25]])

		rejectSamps = np.zeros((paramDim, numSamps))
		rejectDens = np.zeros((numSamps))
		propDist = mm.Gaussian(propMu, propCov)
	
		numAccepts = 0
		numProposed = 0
		M = np.exp(3)
	
		while(numAccepts < numSamps):
			numProposed += 1
			#print(numAccepts)
			#print(numSamps)
			#print(numProposed)
		
			# Propose in the banana space
			mprop = propDist.Sample()
			
			# Evaluate the log target density
			#faulthandler.enable()
	   

			reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
			mprop = mprop.tolist()
			
			for j_num,j in enumerate(kinetic_parameters):
				if mprop[j_num] < 0:
					mprop[j_num] = 0
					kinetic_parameters[j] = 0
				else:
					kinetic_parameters[j] = mprop[j_num]
			print(mprop)
			#print(kinetic_parameters)
			logTgt = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)

			#logTgt = tgtDens.Evaluate([ mprop ])
			
			# Evaluate the log proposal density
			logProp = propDist.LogDensity(mprop)
			print(mprop)
			print(logProp)

			# Compute the acceptance ratio
			alpha = np.exp(logTgt - np.log(M) - logProp)
			
			print(logTgt)

			testRG = mu.RandomGenerator.GetUniform()
			print(testRG)
			print(alpha)
			sys.exit()
			if(testRG < alpha):
				print(alpha)
				rejectSamps[:,numAccepts] = mprop
				rejectDens[numAccepts] = logTgt
				numAccepts += 1

	#############################################################
	############# MARKOV CHAIN MONTE CARLO METHOD ###############
	#############################################################

	if samplingVersion == 'mcmc':
		propMu = np.array([29,0,3.6,11,8,0])
		#propMu = np.zeros((paramDim))
		propCov = 4.0*np.eye(paramDim)
		#propCov = np.array([ [2, 0, 0, 0, 0, 0],[ 0,0.25, 0, 0, 0, 0],[ 0, 0,1, 0, 0, 0],[ 0, 0, 0,2, 0, 0],[ 0, 0, 0, 0,3, 0],[ 0, 0, 0, 0, 0,0.25]])

		mcmcProp = mm.Gaussian(propMu, propCov)

		mcmcSamps = np.zeros((paramDim, numSamps))
		mcmcDens = np.zeros((numSamps))

		currPt = propMu#np.zeros((paramDim))

		reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
		currPt = currPt.tolist()
			
		for j_num,j in enumerate(kinetic_parameters):
			if currPt[j_num] < 0:
				currPt[j_num] = 0
				kinetic_parameters[j] = 0
			else:
				kinetic_parameters[j] = currPt[j_num]
		#print(currPt)

		currLogTgt = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)

		#currLogTgt = tgtDens.Evaluate([currPt])[0]

		numAccepts = 0

		for i in range(numSamps):

			print()
			print('Sample '+str(i+1)+' of '+str(numSamps))

			propSamp = mcmcProp.Sample()

			propSamp2 = propSamp
			reactor_kinetics_input,kinetic_parameters,kin_in,kin_fit = readInput()
			propSamp = propSamp.tolist()
			
			for j_num,j in enumerate(kinetic_parameters):
				if propSamp[j_num] < 0:
					propSamp[j_num] = 0
					kinetic_parameters[j] = 0
				else:
					kinetic_parameters[j] = propSamp[j_num]
			#print(propSamp)
			#print(kinetic_parameters)
			propLogTgt = tap_simulation_function(reactor_kinetics_input,kinetic_parameters,kin_fit)

		
			#propLogTgt = tgtDens.Evaluate([propSamp])[0]
	
			u = np.exp(propLogTgt - currLogTgt)
			#print(mu.RandomGenerator.GetUniform())
			#print(u)
			if(mu.RandomGenerator.GetUniform() < u):
				print('1')
				numAccepts += 1
		
				mcmcSamps[:,i] = propSamp2
		
				currPt = propSamp
				currLogTgt = propLogTgt
			else:
				print('2')
				mcmcSamps[:,i] = currPt
		
		print('Acceptance Rate = %f'%(float(numAccepts)/numSamps))

		plt.figure(figsize=(8,8))
		print(mcmcSamps[0,:])
		print(mcmcSamps[1,:])
		plt.scatter(mcmcSamps[0,:], mcmcSamps[1,:], alpha=0.5)
		plt.title('MCMC Samples1')
		plt.savefig('MCMC_Samples1.png')
		plt.show()

		plt.figure(figsize=(8,8))
		print(mcmcSamps[0,:])
		print(mcmcSamps[2,:])
		plt.scatter(mcmcSamps[0,:], mcmcSamps[2,:], alpha=0.5)
		plt.title('MCMC Samples2')
		plt.savefig('MCMC_Samples1.png')
		plt.show()

		plt.figure(figsize=(8,8))
		print(mcmcSamps[0,:])
		print(mcmcSamps[3,:])
		plt.scatter(mcmcSamps[0,:], mcmcSamps[3,:], alpha=0.5)
		plt.title('MCMC Samples3')
		plt.savefig('MCMC_Samples3.png')
		plt.show()

		plt.figure(figsize=(8,8))
		print(mcmcSamps[0,:])
		print(mcmcSamps[4,:])
		plt.scatter(mcmcSamps[0,:], mcmcSamps[4,:], alpha=0.5)
		plt.title('MCMC Samples4')
		plt.savefig('MCMC_Samples4.png')
		plt.show()

		plt.figure(figsize=(8,8))
		print(mcmcSamps[0,:])
		print(mcmcSamps[5,:])
		plt.scatter(mcmcSamps[0,:], mcmcSamps[5,:], alpha=0.5)
		plt.title('MCMC Samples5')
		plt.savefig('MCMC_Samples5.png')
		plt.show()

		plt.figure(figsize=(16,8))

		plt.plot(mcmcSamps[0,:],label='kf0')
		plt.plot(mcmcSamps[1,:],label='kb0')
		plt.plot(mcmcSamps[2,:],label='kf1')
		plt.plot(mcmcSamps[3,:],label='kb1')
		plt.plot(mcmcSamps[4,:],label='kf2')
		plt.plot(mcmcSamps[5,:],label='kb2')

		plt.legend()

		plt.title('MCMC Chain')
		plt.savefig('MCMC_Chain.png')
		plt.show()

	#############################################################
	############### IMPORTANCE SAMPLING METHOD ################
	#############################################################

	if samplingVersion == 'significance':
		propMu = np.array([29,0,3.6,11,8,0])
		propCov = np.array([ [2, 0, 0, 0, 0, 0],[ 0,0.25, 0, 0, 0, 0],[ 0, 0,1, 0, 0, 0],[ 0, 0, 0,2, 0, 0],[ 0, 0, 0, 0,3, 0],[ 0, 0, 0, 0, 0,0.25]])

		isProp = mm.Gaussian(propMu, propCov)

		isSamps   = np.zeros((paramDim, numSamps))
		isWeights = np.zeros((numSamps))

		for i in range(numSamps):
			isSamps[:,i] = isProp.Sample()
			isWeights[i] = np.exp( tgtDens.Evaluate([isSamps[:,i]])[0] - isProp.LogDensity(isSamps[:,i]) )

	#	logDensVals2 = np.zeros(numSamps)
	#	otherTest = []
	#	for i in range(numSamps):
	#		print(mSamps[:,i])
	#		otherTest.append(tgtDens.Evaluate([ mSamps[:,i] ]))
	#		logDensVals2[i] = tgtDens.Evaluate([ mSamps[:,i] ]) # NOTE: This only works because the jacobian is 1
	#		print('done')
	#	#graphsph.Visualize("EvaluationGraph.png")

	print('pass')
	sys.exit()

call_sim()

