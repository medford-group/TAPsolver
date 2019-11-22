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


warnings.simplefilter(action='ignore', category=FutureWarning)

def tap_simulation_function(reactor_kinetics_input,constants_input):

	"""Core of the TAPsolver simulator.

	Inputs to the function are taken from the input_file parsing script.

	"""

	# Specify the working Tape for optimization / sensitivity analysis
	tape2 = Tape()
	tape2.clear_tape()
	set_working_tape(tape2)
	
	# Generate the output folders (depending on input values)
	kVals = constants_input.copy()
	reac_input = reactor_kinetics_input

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
	for j in r_const:
		r_const[j] = Constant(r_const[j])

	if reac_input['Fit Parameters'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in r_const:
			controls.append(Control(r_const[j]))
			legend_2.append(j)


	# Store input file in output folder
	user_data = pd.read_csv('./input_file.csv',header=None)
	user_data.to_csv(path+'input_file.csv',header=None,index=False)
	original_input_structure = user_data.copy()
	
	# Control tolerence and error level
	set_log_level(30)
	tol = 1E-14

	# Define solving method
	#theta = 0.5

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
	#dt = 0.0001
	#dt = reac_input['Pulse Duration']/reac_input['Time Steps']
	
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

	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],reac_input['Reactor Type'],reac_input['Number of Active Sites'],reac_input['Number of Inerts'],reac_input['Advection'],False)

	####### step 1
	# Define the initial, non-refined mesh
	mesh = UnitIntervalMesh(int(reac_input['Mesh Size']))

	# Mesh refinement process
	catalystRefinement = int(reac_input['Catalyst Mesh Density'])
	cfDict = {}

	for jayz in range(0,catalystRefinement+1):
		class thin_zoneTest(SubDomain):
			def inside(self, x, on_boundary):
				return between(x[0], (((1-cat_location) - 0.5*frac_length), ((1-cat_location) + 0.5*frac_length)))

		thin_zoneTest = thin_zoneTest()
		cfDict[jayz] = MeshFunction("bool",mesh,1)
		
		thin_zoneTest.mark(cfDict[jayz],jayz)
		mesh = refine(mesh,cfDict[jayz],True)

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

	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (((1-cat_location) - 0.5*frac_length)-1/reac_input['Mesh Size'], ((1-cat_location) + 0.5*frac_length)+1/reac_input['Mesh Size']))
			
	class singlePoint(SubDomain):
		def inside(self,x,on_boundary):
			return between(x[0], (((1-cat_location) - 1/totalNumCells), ((1-cat_location) + 1/totalNumCells)))

	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,1)
	#domains.set_all(0)
	thin_zone.mark(domains,1)

	newBoundaries = domains.mesh().coordinates().transpose().tolist()[0]
	additionalCells = len(newBoundaries[int(reac_input['Mesh Size']+1):])

	centralPoint = singlePoint()
	boundary_parts0 = MeshFunction("size_t", mesh,0)
	boundary_parts0.set_all(0)
	centralPoint.mark(boundary_parts0, 1)

	dx = Measure("dx",subdomain_data=domains)
	dT = Measure("dx",subdomain_data=boundary_parts0)
	
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
	
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)
	
	dz = 1/reac_input['Mesh Size']

	class integration_section(SubDomain):
		def inside(self, x, on_boundary):
			#return between(x[0], (1-dz,1.0))
			return between(x[0], (1-dz,1.0))
	
	right = CompiledSubDomain("near(x[0], 1.)")

	boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
	right.mark(boundary_parts, 1)

	# Define the new graphs
	monitored_gas = necessary_values['molecules_in_gas_phase']
	
	bcs = defineBCs(reac_input['Reactor Type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['Number of Reactants'],all_molecules,reac_input['Pulse Ratio'],boundary_L,boundary_R,reac_input['Number of Inerts'])
	
	fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])
	ax2.set_xlim(0,reac_input['Pulse Duration'])


	# Define controls only if needed for differentiation based analysis
	if reac_input['Fit Inert'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
			controls.append(Control(Dout[j]))
			legend_2.append(j)

	objSpecies = list(reac_input['Objective Species'].split(','))

	try:
		theta = 1
		Ftemp = eval(necessary_values['F'])
		theta = 0.5
		F = eval(necessary_values['F'])
	except NameError:
		errorOutput(reac_input['reactions_test'])
	
	if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true'):
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
		path_2 = reac_input['Output Folder Name']+'_folder/sensitivity/'
		generateFolder(path_2)
		
		
		path_molecules = path_2+reac_input['Sensitivity Parameter']
		generateFolder(path_molecules)
		sensFolder = path_molecules

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

		
	to_flux = fluxGeneration(reac_input['Reactor Type'],len(legend_label),reac_input['Number of Reactants'],reac_input['Reference Pulse Size'],D,eb,dx_r,reac_input['Reactor Radius'],dx2_r,reac_input['Scale Output'])
	storeDataFunc(reac_input['Store Outlet Flux'],reac_input['Output Folder Name'])
	storeSens(reac_input['Sensitivity Analysis'],reac_input['Output Folder Name'],monitored_gas,legend_label)

	J = derivative(F,u)
	Jtemp = derivative(Ftemp,u)

	# Define the variational problem solver with newton or constrained
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
		
	elif reac_input['Variational Solver'].lower() == 'newton':
		problem = NonlinearVariationalProblem(F,u,bcs,J)
		solver = NonlinearVariationalSolver(problem)

		problemtemp = NonlinearVariationalProblem(Ftemp,u,bcs,Jtemp)
		solvertemp = NonlinearVariationalSolver(problemtemp)


		solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10

		solvertemp.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
		solvertemp.parameters["newton_solver"]["absolute_tolerance"] = 1e-10

		######################################################################
		snes_solver_parameters = {"nonlinear_solver": "snes","snes_solver": {"linear_solver": "lu","line_search":'basic',"maximum_iterations": 10,"report": False,"error_on_nonconvergence": False}}

		
		lower = Function(V)
		upper = Function(V) 
		
		ninfty = Function(V); ninfty.vector()[:] = 0
		pinfty = Function(V); pinfty.vector()[:] =  np.infty

		
		problem2 = NonlinearVariationalProblem(F,u,bcs,J)
		
		problem2.set_bounds(ninfty,pinfty)

		solver2 = NonlinearVariationalSolver(problem)
		solver2.parameters.update(snes_solver_parameters)
		########################################################################


	### Handle a pulse intensity that changes with each pulse ###
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

	### Begin simulation for each desired pulse ###
	current_reac_set = 1
	tot_objective = 0

	if reac_input['Thin-Zone Analysis'].lower() == 'true':
		rateStrings = rateEqs(rate_array,rev_irr)
	
	if reac_input['RRM Analysis'].lower() == 'true':
		rrmStringsThin = rrmEqs(rate_array,rev_irr,'dx(1)')	
		rrmStringsPoint = rrmEqs(rate_array,rev_irr,'dT(1)')

	if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true'):
		for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
			if objSpecies[k_fitting] == '1':
				for timeStep in range(0,len(output_fitting[legend_label[k_fitting]]['times'])):
					output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)

	if reac_input['Fit Inert'].lower() == 'true':
		for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
			for timeStep in range(0,len(output_fitting[legend_label[-1]]['times'])):
				output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)

	if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':

		c = r_const[reac_input['Sensitivity Parameter']]
		c.tlm_value = r_const[reac_input['Sensitivity Parameter']]

		SV_du = FunctionSpace(mesh,P1)
		Sw_new = Expression('A',A=Constant(1),degree=0)
		Sw_new2 = interpolate(Sw_new,SV_du)
		Sw3 = project(Sw_new2,SV_du)


	if reac_input['Sensitivity Analysis'].lower() == 'true':

		sensFuncs = {}

		for k_gasses in range(0,len(necessary_values['reactants'])):
			sensFuncs[str(k_gasses)] = []	
	

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


	for k_pulse in range(0,int(reac_input['Number of Pulses'])):

		start_time = time.time()
		print("")
		print("Simulation Status: "+"Pulse #"+str(k_pulse+1))

		species_pulse_list = reactant_feed[current_reac_set-1]
		
		if type(reactant_time[current_reac_set-1]) == list:
			species_time = reactant_time[current_reac_set-1]
		else:
			species_time = reactant_time

		if reac_input['Knudsen Test'].lower() == 'true':
			knudsenTest(legend_label[(int(len(legend_label))-int(reac_input['Number of Inerts'])):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],reac_input['Reference Pulse Size'],species_pulse_list[(int(len(legend_label))-int(reac_input['Number of Inerts'])):])
 
		#dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
		dt = reac_input['Pulse Duration']/reac_input['Time Steps']
		
		for kTimeStep,kTime in enumerate(species_time.copy()):
		#for kTimeStep,kTime in enumerate(reactant_time.copy()):
			tNew = dt*round(float(kTime)/dt)
			species_time[kTimeStep] = round(tNew,6)		
		
		### Redefine the lists tracking all the information ###
		sensitivity_output = {}
		RRM_der = {}
		
		for k_sens in range(monitored_gas):
			sensitivity_output[k_sens] = []

		for k_sens in range(len(necessary_values['reactants'])):
			RRM_der[k_sens] = []
				
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

		t = 0
		if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':
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




		transTest1 = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest1*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0
		transTest1 = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest1*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0
		transTest1 = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-2
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest1*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0
		transTest1 = mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-3
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest1*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0


		transTest2 = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))+1
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest2*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0
		transTest2 = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))+2
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest2*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0
		transTest2 = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))+3
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest2*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0
		transTest2 = int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))+4
		if ',' in str(reac_input['Initial Surface Composition']):
			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
				u_n.vector()[transTest2*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = 0




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

			if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':

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
									jfunc_2 = assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
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

				#dt = solverIteration(dt,reac_input['Solver Method'],solver,dk,1.5,1.1)
				#print(t)
				#time.sleep(0.1)

				try:
					if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
						solver.solve()
					else:
						
						if t > 0.0011+timeDiff:
						#if t > 0.0011:
							#if t > 0 and k_pulse < 4:
							solver.solve(annotate = False)
							#else:
							#	solver2.solve(annotate = False)
						else:
							solvertemp.solve(annotate=False)

							if round(t,6) == 0.001+timeDiff:
							#if round(t,6) == 0.001:
								dt = reac_input['Pulse Duration']/reac_input['Time Steps']
								dk.assign(dt)
								u_n.assign(u)
								solver.solve()
							


					if reac_input['Thin-Zone Analysis'].lower() == 'true':
						for jk in range(0,all_molecules+1):
							new_values = [] 
							cat_values = []

							#for z in range(0,100):
							#### for z in range(4500,int((reac_input['Mesh Size'])+1)+additionalCells):
							#for z in range(0,100):
							#for z in range(4000,int((reac_input['Mesh Size'])-1)+additionalCells+3):
							#print(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1)
							#print((int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1)+additionalCells)
							#sys.exit()
							#for z in range(0,int((reac_input['Mesh Size'])+1)+additionalCells):
							
							#for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1,int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])):
							if additionalCells == 0:
								#print(mp.int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1)
								#print(int((cat_location + 0.5*frac_length)*reac_input['Mesh Size']))
								#sys.exit()
								for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
								#for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']),int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])):
									try:
										#value_cat = u_n.vector().get_local()[z*(all_molecules)-(all_molecules)+jk]
										value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
										cat_values.append(value_cat)
							
									except AttributeError:
										value_cat = u_n.vector().get_local()[z*(all_molecules)]
										sys.exit()
										#cat_values[0] += value_cat*dx_r
										cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))
							else:
								#print(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']))
								#print(int((cat_location + 0.5*frac_length)*reac_input['Mesh Size']))
								#sys.exit()
								for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):

								#for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size']),int((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))):
									#print(z)
									try:
										#value_cat = u_n.vector().get_local()[z*(all_molecules)-(all_molecules)+jk]
										value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
										cat_values.append(value_cat)
							
									except AttributeError:
										value_cat = u_n.vector().get_local()[z*(all_molecules)]
										sys.exit()
										#cat_values[0] += value_cat*dx_r
										cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))
								###sys.exit()

							cat_data['conVtime_'+str(jk)].append(cat_values)

				except RuntimeError:
					######### 
					#print('runtimeError')
					if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
						solver2.solve()
					else:
						solver2.solve(annotate = False)
				
			else:
				dt = 0.0001
				dk.assign(dt)

				if round(t,6) in species_time:
					timeDiff = round(t,6)

				if reac_input['Reactor Type'] == 'tap':
					
					if reac_input['Fit Parameters'].lower() == 'true':
						if t == 0:
							for k in range(0,int(reac_input['Number of Reactants'])):
								#!#!#!#!##!#!#u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']+1))+k-(all_molecules)] = -1e-10#list_species_pulse[k]
								u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+additionalCells))+k-(all_molecules)] = -1e-10#list_species_pulse[k]

					if ',' in str(species_pulse_list):
						for k in range(0,int(reac_input['Number of Reactants'])):
							
							if species_time[k] == round(t,6):
								u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells)))+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]
								
					else:
						
						u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']-10)-1)-(all_molecules)] = float(species_pulse_list)*Inert_pulse_conc

					#sys.exit()
					for k_step in range(0,int(reac_input['Number of Inerts'])):
						if species_time[-1-k_step] == round(t,6):
							u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells))-1-k_step)] = float(species_pulse_list[-1-k_step])*Inert_pulse_conc###??? Added the porosity contribution
							
					if reac_input['Thin-Zone Analysis'].lower() == 'true':
						for jk in range(0,all_molecules+1):
							new_values = [] 
							cat_values = []

							if additionalCells == 0:
								for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
						
									try:
										value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
										cat_values.append(value_cat)
									except AttributeError:
										value_cat = u_n.vector().get_local()[z*(all_molecules)]
										sys.exit()
										
										cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))
							else:
								for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
									try:
										value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
										cat_values.append(value_cat)
									except AttributeError:
										value_cat = u_n.vector().get_local()[z*(all_molecules)]
										sys.exit()
										cat_values[0] += value_cat*dx_r/(2**(int(reac_input['Catalyst Mesh Density'])))							

							cat_data['conVtime_'+str(jk)].append(cat_values)

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
							#print(z)
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
						#if t > 0:
						solvertemp.solve()
						#else:
						#	solvertemp.solve(annotate = False)
						#	u_n.assign(u)
						#	solver.solve()

					else:
						#if t > 0:
						solvertemp.solve(annotate = False)
						#else:
						#	solvertemp.solve(annotate = False)
						#	u_n.assign(u)
						#	solver.solve()

				except RuntimeError:
					if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
						solver2.solve()
					else:
						solver2.solve(annotate = False)
					print('Time Step Failure')
				
			if reac_input['Sensitivity Analysis'].lower() == 'true':
				
				u_final = u.split(deepcopy=False)
				should_it = int(round(t*reac_input['Time Steps']/reac_input['Pulse Duration'],0))
				
				for k in range(0,monitored_gas):
					new_val = (( u.vector().get_local()[(all_molecules)+k]))
					u_graph_data['conVtime_'+str(k)].append((new_val))

				for kGasses in range(0,len(necessary_values['reactants'])):
					sensFuncs[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dP(1)))
					
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
			print()
			start_time = time.time()
			print('Evaluating Tape with Tangent Linear Method. Could take some time.')
			tape2.evaluate_tlm()

		if reac_input['Sensitivity Analysis'].lower() == 'true':
			
			for numEachSens,eachSens in enumerate(u_graph_data):
				np.savetxt(sensFolder+'/c_'+legend_label[numEachSens]+'.csv',u_graph_data[eachSens],delimiter=",")

			for numEachSens,eachSens in enumerate(sensFuncs):
				newList = []
				for kSensNum, kSens in enumerate(sensFuncs[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)
				np.savetxt(sensFolder+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
			
		if reac_input['RRM Analysis'].lower() == 'true':

			for numEachSens,eachSens in enumerate(thinSensFunc):
				newList = []
				for kSensNum, kSens in enumerate(thinSensFunc[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)
				np.savetxt(rrmFolder+'/thinValue'+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")#

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


		if reac_input['Thin-Zone Analysis'].lower() == 'true':
			if int(reac_input['Number of Pulses']) == 1:
				for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					
					np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
					rateTest = eval(rateStrings[j_species])
					np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")
			else:
				pulse_path = './'+reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(k_pulse+1)+'/'
				generateFolder(pulse_path)
				for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
					rateTest = eval(rateStrings[j_species])
					np.savetxt(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")
			np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/Inert.csv', np.array(cat_data['conVtime_'+str(all_molecules-int(reac_input['Number of Inerts']))]), delimiter=",")

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

		if reac_input['Uncertainty Quantification'].lower() == 'true':
			start_time = time.time()
			print()
			print('Calculating hessian. Could take some time.')
			
			control = Control(c)
			rf_2 = ReducedFunctional(jfunc_2, control)

			dJdm = c._ad_dot(rf_2.derivative())
			Hm = c._ad_dot(rf_2.hessian(c))
			print(processTime(start_time))
			np.savetxt(hessFolder+'/'+reactor_kinetics_input['Sensitivity Parameter']+'.csv',[Hm],delimiter=",")#

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
		
		if reac_input['Infinite Inert'].lower() == 'true':
			if k_pulse == 0:
				analyticalTiming = np.arange(0, reac_input['Time Steps']*dt, dt).tolist()
				for kjc in range(0,monitored_gas):

					outlet = []
					
					
					
					if reac_input['Scale Output'].lower() == 'true':
						factor = 1
					else:
						factor = float(species_pulse_list[kjc])*reac_input['Reference Pulse Size']
					
					for k in range(0,int(reac_input['Time Steps'])+1):
						analyticalValue = 0
						for n in range(0,50):
							analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[kjc])*(D[kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
						if analyticalValue < 0:
							outlet.append(0)
						else:
							outlet.append(factor*(D[kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
					np.savetxt('./analyticalCO19000.csv', outlet, delimiter=",")

					try:
						ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
						#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
					except ValueError:
						outlet = outlet[:-1]
						ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
						#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)


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
				#print('what?')
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

	return graph_data, legend_label, necessary_values['reactants']

def call_sim():
	
	reactor_kinetics_input,kinetic_parameters,kin_in = readInput()
	
	if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' or reactor_kinetics_input['RRM Analysis'].lower() == 'true' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
		for parameters in kinetic_parameters:
			reactor_kinetics_input,kinetic_parameters,kin_in = readInput()
			reactor_kinetics_input['Sensitivity Parameter'] = parameters
			reactor_kinetics_input['Display Graph'] = 'FALSE'
			if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or reactor_kinetics_input['Fit Inert'].lower() == 'true':
				print('')
				print('')
				
				print('Running the Sensitivity/RRM Analysis and Parameter Fitting methods simultaniously is not possible due to conflicts between the tangent linear and adjoint methods.')
				print('Please run again with one of these methods excluded.')
				sys.exit()
			print('')
			print('Processing '+parameters)
			graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)
	else:
		graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)
	
call_sim()
