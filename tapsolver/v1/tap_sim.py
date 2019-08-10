from fenics import *
import dolfin
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

fenics_version = dolfin.__version__

rrmProof = True

if fenics_version == '2017.2.0':
	fen_17 = False#True

else:
	fen_17 = False
	print('You are working with a newer version of FEniCS (beyond 2017). Some methods of analysis could be limited. Dolfin methods could also be limited to a smaller number of time steps.')
	print('If "fit gif" has been called, it might not properly generate a gif')
	time.sleep(2)

def tap_simulation_function(reactor_kinetics_input,constants_input):
	tape = Tape()
	#tape.clear_tape()
	set_working_tape(tape)
	### Initial forms for storing results and reading in the information ###
	simulation_time_list = []
	kVals = constants_input.copy()
	reac_input = reactor_kinetics_input
	r_const = constants_input

	### Generate the folders for information storage ###
	path = './'+reac_input['Output Folder Name']+'_folder/'
	generate_folder(path)
		
	path_3 = './'+reac_input['Output Folder Name']+'_folder/flux_data/'
	generate_folder(path_3)
	
	if reac_input['Thin-Zone Analysis'].lower() == 'true':
		path_4 = './'+reac_input['Output Folder Name']+'_folder/thin_data/'
		generate_folder(path_4)

	path_5 = './'+reac_input['Output Folder Name']+'_folder/graphs/'
	generate_folder(path_5)

	if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true':
		path_6 = './'+reac_input['Output Folder Name']+'_folder/fitting/'
		generate_folder(path_6)

	### Declare and define the constants ###
	for j in r_const:
		r_const[j] = Constant(r_const[j])

	### Define controls only if needed for differentiation based analysis ###
	if reac_input['Fit Parameters'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in r_const:
			controls.append(Control(r_const[j]))
			legend_2.append(j)


	### Copy and store the input_csv file in the new folder (to keep organized) ###
	user_data = pd.read_csv('./input_file.csv',header=None)
	user_data.to_csv(path+'input_file.csv',header=None,index=False)
	original_input_structure = user_data.copy()
	### Control the output information from FEniCS ###
	parameters["std_out_all_processes"] = False
	cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}
	
	if fen_17 == True:
		set_log_active(False)
	import warnings
	warnings.filterwarnings("ignore", category=DeprecationWarning)
	tol = 1E-14
	theta = reac_input['Theta']

	### Used to define the initial surface composition over the simulation ###
	if ',' in str(reac_input['Initial Surface Composition']):
		rangesurface_species = list(reversed(reac_input['Initial Surface Composition'].split(',')))
		reac_input['Initial Surface Composition'] = list(map(float, rangesurface_species))

	else:
		rangesurface_species = reac_input['Initial Surface Composition']
	
	### Initialize the grid system, time step size, pulse size and diffusion coefficients ###
	r_param, dx_r, dx2_r, frac_length, cat_location = establish_grid_system(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['Mesh Size'])
	
	cat_location = reac_input['Catalyst Location']
	
	dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
	eb = np.array((reac_input['Void Fraction Inert'],reac_input['Void Fraction Catalyst'],reac_input['Void Fraction Inert']))

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

	### Define the dimensions of the reactor ###
	ca = (reac_input['Reactor Radius']**2)*3.14159 

	point_volume = dx_r * ca * eb[0] ###??? Should this include the voidage?
	

	Inert_pulse_conc = reac_input['Reference Pulse Size']/(point_volume) ###??? Similar issue to the volume previously used

	dt = reac_input['Pulse Duration']/reac_input['Time Steps']

	### Construct the reaction equation in a form that Fenics can understand ###
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],reac_input['Reactor Type'],reac_input['Number of Active Sites'],reac_input['Number of Inerts'],reac_input['Advection'],False)
	
	### Declaring the trial / test functions in fenics and defining the finite elements ###
	mesh = UnitIntervalMesh(int(reac_input['Mesh Size'])) 
	P1 = FiniteElement('P',mesh.ufl_cell(),1)

	if reac_input['reactions_test'] != ['INERT_ONLY']:
		test_new = eval(necessary_values['element'])
		element = MixedElement(test_new)
		V = FunctionSpace(mesh,element)
		V_du = FunctionSpace(mesh,P1)
	else:
		V = FunctionSpace(mesh,P1)
		V_du = FunctionSpace(mesh,P1)
	#all_molecules = necessary_values['gas_num']+necessary_values['surf_num']+reac_input['Number of Inerts']
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

	graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data = initialize_variable_dictionaries(necessary_values,all_molecules,V,u,u_n)

	### Define the subdomain for the catalyst region
	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (((1-cat_location) - 0.5*frac_length), ((1-cat_location) + 0.5*frac_length)))

	class singlePoint(SubDomain):
		def inside(self,x,on_boundary):
			return between(x[0], (((1-cat_location) - 1/reac_input['Mesh Size']), ((1-cat_location) + 1/reac_input['Mesh Size'])))


	
	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	thin_zone.mark(domains,1)
	
	centralPoint = singlePoint()
	boundary_parts0 = MeshFunction("size_t", mesh,0)
	boundary_parts0.set_all(0)
	centralPoint.mark(boundary_parts0, 1)


	dx = Measure("dx")[domains]
	dT = Measure("dx")[boundary_parts0]

	### Define the boundaries at each end of the reactor ###
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
	#ds = Measure("ds", subdomain_data=boundary_parts)
	
	monitored_gas = necessary_values['molecules_in_gas_phase']
	
	bcs = define_boundary_conditions(reac_input['Reactor Type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['Number of Reactants'],all_molecules,reac_input['Pulse Ratio'],boundary_L,boundary_R,reac_input['Number of Inerts'])
	
	### Initialize the graphs and set the graphical parameters ###
	fig2,ax2,legend_label,header,colors = establish_output_graph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])
	#ax2.set_ylim(0,1.75)
	ax2.set_xlim(0,reac_input['Pulse Duration'])
	### Evaluate the reaction/diffusion expression for FEniCS to use ###


	### Define controls only if needed for differentiation based analysis ###
	if reac_input['Fit Inert'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
			controls.append(Control(Dout[j]))
			legend_2.append(j)

	objSpecies = list(reac_input['Objective Species'].split(','))

	try:
		F = eval(necessary_values['F'])
	except NameError:
		error_output(reac_input['reactions_test'])
	
	if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true'):
		try:
			if type(reac_input['Objective Points']) == float:
				output_fitting = exp_data_fitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
			elif reac_input['Objective Points'] == 'all':
				output_fitting = every_point_fitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)

			else:
				print('Objective Points defined incorrectly')
				sys.exit()
		except TypeError:
			print('Objective Point Input Is Not Valid')
			sys.exit()

	if reac_input['Fit Inert'].lower() == 'true':
		try:
			if type(reac_input['Objective Points']) == float:
				output_fitting = exp_data_fitting(legend_label[int(reac_input['Number of Inerts']):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'])
			elif reac_input['Objective Points'] == 'all':
				output_fitting = every_point_fitting(legend_label[int(reac_input['Number of Inerts']):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'])
			else:
				print('Objective Points defined incorrectly')
				sys.exit()
		except TypeError:
			print('Objective Point Input Is Not Valid')
			sys.exit()
	sens_time_list = []
	

	if reac_input['Sensitivity Analysis'].lower() == 'true':
		path_2 = reac_input['Output Folder Name']+'_folder/sensitivity/'
		generate_folder(path_2)
		
		
		path_molecules = path_2+reac_input['Sensitivity Parameter']
		generate_folder(path_molecules)
		sensFolder = path_molecules

	if reac_input['Uncertainty Quantification'].lower() == 'true':
		path_2 = reac_input['Output Folder Name']+'_folder/UQ/'
		generate_folder(path_2)
		
	if reac_input['RRM Analysis'].lower() == 'true':
		path_2 = reac_input['Output Folder Name']+'_folder/RRM_analysis/'
		generate_folder(path_2)
		
		path_molecules = path_2+reac_input['Sensitivity Parameter']
		generate_folder(path_molecules)
		rrmFolder = path_molecules

		generate_folder(path_molecules+'/pointValue')
		generate_folder(path_molecules+'/thinValue')

		
	to_flux = flux_generation(reac_input['Reactor Type'],len(legend_label),reac_input['Number of Reactants'],reac_input['Reference Pulse Size'],D,eb,dx_r,reac_input['Reactor Radius'],dx2_r,reac_input['Scale Output'])
	store_data_func(reac_input['Store Outlet Flux'],reac_input['Output Folder Name'])
	store_sens_analysis(reac_input['Sensitivity Analysis'],reac_input['Output Folder Name'],monitored_gas,legend_label)

	#class Plane(SubDomain):
	#	def inside(self, x, on_boundary):
	#		return x[0] > 1.0 - DOLFIN_EPS
	#solver = PETScTAOSolver()
	J = derivative(F,u)

	#### Might want to consider adding this later (the hessian, that is!)
	###H = derivative(J,u)

	if reac_input['Variational Solver'].lower() == 'constrained':
		snes_solver_parameters = {"nonlinear_solver": "snes","snes_solver": {"linear_solver": "lu","maximum_iterations": 40,"report": False,"error_on_nonconvergence": False}}
		W = FunctionSpace(mesh, "P", 1)
		#solver.parameters.update(snes_solver_parameters)
		testV_du = FunctionSpace(mesh,P1)
		testw_new = Expression('A',A=Constant(1),degree=0)#Expression("1", degree=0)
		#testw_new2 = 
		a_min = interpolate(testw_new,testV_du)
		
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

		#solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
		#solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10



	####solver.parameters["nonlinear_solver"] = 'snes'
	####solver.parameters["snes_solver"]["maximum_iterations"] = 100
	####solver.parameters["snes_solver"]["report"] = False
	#constraint_l = Expression('A',A=Constant(1),degree=0) 
	#umin = interpolate(constraint_l, FunctionSpace(mesh,P1))
	#sys.exit()
	# Create the PETScTAOSolver
	####solver = PETScTAOSolver(problem)

	# Set some parameters
	####solver.parameters["method"] = "tron"
	####solver.parameters["monitor_convergence"] = True
	####solver.parameters["report"] = True
	
	
	# Uncomment this line to see the available parameters
	####info(parameters, True)
	
	# Parse (PETSc) parameters
	####parameters.parse()
	####sys.exit()

	#dolfin.parameters["nonlinear_solver"]["test_gradient"] = True
	#dolfin.parameters["nonlinear_solver"] = 'snes'
	#solver.parameters["snes_solver"]["maximum_iterations"] = 50
	#solver.parameters["snes_solver"]["report"] = False

	#solver  = NonlinearVariationalSolver(problem)
	#solver.parameters["nonlinear_solver"] = 'snes'

	#a_max = Function(interpolate(Constant(1.0), V))


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

	### legend_label[int(reac_input['Number of Inerts']):]
	if reac_input['Fit Inert'].lower() == 'true':
		for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
			for timeStep in range(0,len(output_fitting[legend_label[-1]]['times'])):
				output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)

	if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':

		c = r_const[reac_input['Sensitivity Parameter']]
		c.tlm_value = r_const[reac_input['Sensitivity Parameter']]

		SV_du = FunctionSpace(mesh,P1)
		Sw_new = Expression('A',A=Constant(1),degree=0)#Expression("1", degree=0)
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

		#thinSensFuncRate['CO'] = []
		#pointSensFuncRate['CO'] = []
		#thinSensFuncRate['CO2'] = []
		#pointSensFuncRate['CO2'] = []

	sensAll = []
	sensAll2 = []
	#c = r_const["kf0"]
	for k_pulse in range(0,int(reac_input['Number of Pulses'])):

		### Clear the tape for the adjoint (will lead to challenges with fitting parameters over multiple pulses) ###
		
		print("")
		print("Simulation Status: "+"Pulse #"+str(k_pulse+1))

		species_pulse_list = reactant_feed[current_reac_set-1]
		species_time = reactant_time[current_reac_set-1]

		if reac_input['Knudsen Test'].lower() == 'true':
			knudsenTest(legend_label[(int(len(legend_label))-int(reac_input['Number of Inerts'])):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],reac_input['Reference Pulse Size'],species_pulse_list[(int(len(legend_label))-int(reac_input['Number of Inerts'])):])
			time.sleep(3)


		### Incase the time step was altered during a previous pulse, just reset to the original values ### 
		dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
		dt = reac_input['Pulse Duration']/reac_input['Time Steps']
		#!#!print(reactant_time)
		start_time = time.time()
		for kTimeStep,kTime in enumerate(reactant_time.copy()):
			tNew = dt*round(kTime/dt)
			reactant_time[kTimeStep] = round(tNew,6)
		### Redefine the lists tracking all the information ###
		sensitivity_output = {}
		RRM_der = {}
		
		for k_sens in range(monitored_gas):
			sensitivity_output[k_sens] = []

		for k_sens in range(len(necessary_values['reactants'])):
			RRM_der[k_sens] = []

		#for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
		#	#print(necessary_values['reactants'][j_species])
		#	np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
		#for k_rrm in range()
				
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
		
		for z_gasses in range(0,all_molecules):
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

		### Used to define the function to be minimized between the experimental points and synthetic points ###
		w_new = Expression("1", degree=0)
		w_new2 = interpolate(w_new,V_du)
		#W = VectorFunctionSpace(mesh, 'P', 1)
		#w = Function(W)
		#w.vector()[:] = 10

		x_dim = list(range(0, int(reac_input['Mesh Size'])+1))
		
		cum_sum = 0

		while t <= reac_input['Pulse Duration']:
			### Store the results of each iteration ###
			graph_data['timing'].append(t)
			for k in range(0,monitored_gas):
				new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k]))
				graph_data['conVtime_'+str(k)].append((new_val))
			#for kjc in range(all_molecules-int(reac_input['Number of Inerts']),all_molecules):
			for kjc in range(0,int(reac_input['Number of Inerts'])):
				new_val = ((to_flux[monitored_gas+kjc]*u.vector().get_local()[2*(all_molecules+1)-2-(int(reac_input['Number of Inerts'])-kjc)]))
				graph_data['conVtime_'+str(all_molecules-(int(reac_input['Number of Inerts'])-kjc))].append((new_val))
			

			mesh_size = reac_input['Mesh Size']

			### If you are trying to fit parameters, then check if the point must be defined as a part of the objective function ###
			if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true':

				for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
					if objSpecies[k_fitting] == '1':
						if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
						
							c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
							slope = (-c_exp)/(1/mesh_size)
							intercept = c_exp - ((1-(1/mesh_size))*slope)
							w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)#Expression("1", degree=0)
							w_new2 = interpolate(w_new,V_du)
							w3 = project(w_new2,V_du)#output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]*
							#test_meaning = assemble(inner(u_n[k_fitting] + output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))] ,u_n[k_fitting] + output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))])*ds(1))

							try:
								if legend_label[k_fitting] != 'Inert':
									jfunc_2 += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
									#tot_objective += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
									#tot_objective += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
								
								else:
									pass

							except UnboundLocalError:
								
								if legend_label[k_fitting] != 'Inert':
									jfunc_2 = assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
									#tot_objective += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
									#tot_objective += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
								else:
									pass
							#print(tot_objective)

			### len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label[int(reac_input['Number of Inerts']):]


			if reac_input['Fit Inert'].lower() == 'true':

				for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
					if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
						c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
						slope = (-c_exp)/(1/mesh_size)
						intercept = c_exp - ((1-(1/mesh_size))*slope)
						w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)#Expression("1", degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)#output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]*
						#test_meaning = assemble(inner(u_n[k_fitting] + output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))] ,u_n[k_fitting] + output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))])*ds(1))

						try:
							if legend_label[k_fitting] != 'Inert':
								jfunc_2 += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
								#tot_objective += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
								#tot_objective += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
								
							else:
								pass

						except UnboundLocalError:
							
							if legend_label[k_fitting] != 'Inert':
								jfunc_2 = assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
								#tot_objective += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))
								#tot_objective += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
							else:
								pass

					
			if reac_input['Thin-Zone Analysis'].lower() == 'true':
				for jk in range(0,all_molecules):
					new_values = [] 
					cat_values = []

					for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1,int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])):
						
						try:
							value_cat = u_n.vector().get_local()[z*(all_molecules)-(all_molecules)+jk]
							cat_values.append(value_cat)
							
						except AttributeError:
							value_cat = u_n.vector().get_local()[z*(all_molecules)]
							sys.exit()
							cat_values[0] += value_cat*dx_r

					cat_data['conVtime_'+str(jk)].append(cat_values)

			if round(t,6) not in reactant_time:

				dt = solver_iteration(dt,reac_input['Solver Method'],solver,dk,1.5,1.1)
				
			else:
				if reac_input['Reactor Type'] == 'tap':
					
					if reac_input['Fit Parameters'].lower() == 'true':
						if t == 0:
							for k in range(0,int(reac_input['Number of Reactants'])):
									u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']+1))+k-(all_molecules)] = -1e-10#list_species_pulse[k]

					if ',' in str(species_pulse_list):
						for k in range(0,int(reac_input['Number of Reactants'])):
							
							if reactant_time[k] == round(t,6):

								u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']+1))+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]

					else:
						u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']+1)-1)-(all_molecules)] = float(species_pulse_list)*Inert_pulse_conc
					
					for k_step in range(0,int(reac_input['Number of Inerts'])):
						if reactant_time[-1-k_step] == round(t,6):
							u_n.vector()[int((all_molecules)*(reac_input['Mesh Size']+1)-1-k_step)] = float(species_pulse_list[-1-k_step])*Inert_pulse_conc###??? Added the porosity contribution
					
				if k_pulse == 0:
					for z in range(mp.ceil((cat_location - 0.5*frac_length)*reac_input['Mesh Size'])-1,int((cat_location + 0.5*frac_length)*reac_input['Mesh Size'])):
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

				
				dt = solver_iteration(dt,reac_input['Solver Method'],solver,dk,1.5,1.1)
				
			if reac_input['Sensitivity Analysis'].lower() == 'true':
				time_size = time.time()
				u_final = u.split(deepcopy=False)
				should_it = int(round(t*reac_input['Time Steps']/reac_input['Pulse Duration'],0))
				#sensAll.append(assemble(inner(u_final[k_step],u_final[k_step]))*dx)
				### #Js.append(assemble( ( inner(u[0], w3) )* dx))		
				#sensAll.append(assemble(inner(u, u) * dx))
				
				#Store the actual value of the concentration
				for k in range(0,monitored_gas):
					new_val = (( u.vector().get_local()[(all_molecules)+k]))
					u_graph_data['conVtime_'+str(k)].append((new_val))

				#Store the 
				for kGasses in range(0,len(necessary_values['reactants'])):
					sensFuncs[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dP(1)))
					
				#for kGasses in range(0,necessary_values['molecules_in_gas_phase']):
			if reac_input['RRM Analysis'].lower() == 'true':
				## Store the values of the thin zone analysis
				for kGasses in range(0,len(necessary_values['reactants'])):
					thinSensFunc[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dx(1)))
					pointSensFunc[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3) )* dT(1)))

				#Js.append(assemble((inner(ufl.ln(ufl.exp(-dG)*u_n[0] - ufl.exp(-dG)*u_n[1]), w3) )* dx))
				#sys.exit()
				for kGasses in range(0,len(necessary_values['reactants'])):
					thinSensFuncRate[str(kGasses)].append(eval(rrmStringsThin[kGasses]))
					pointSensFuncRate[str(kGasses)].append(eval(rrmStringsPoint[kGasses]))
				
				#thinSensFuncRate['CO'].append(assemble( ( inner(-r_const["kf0"]*u[0]*u[4] + r_const["kb0"]*u[2], Sw3) )* dx(1)))
				#pointSensFuncRate['CO'].append(assemble( ( inner(-r_const["kf0"]*u[0]*u[4] + r_const["kb0"]*u[2], Sw3) )* dT(1)))
				#thinSensFuncRate['CO2'].append(assemble( ( inner(r_const["kf1"]*u[2]*u[3], Sw3) )* dx(1)))
				#pointSensFuncRate['CO2'].append(assemble( ( inner(r_const["kf1"]*u[2]*u[3], Sw3) )* dT(1)))

				#print(type(assemble((to_flux[k_step]*u_final[k_step])*(to_flux[k_step]*u_final[k_step])*dP(1))))
				#if should_it <= 6 and should_it%5 == 0:
				#	direction = project(Constant(1),V_du)
				#	all_together = [direction,direction,direction,direction]
				#	for k_step in range(0,monitored_gas):
				#		testAgain = [Control(Constant(1)),Control(Constant(1)),Control(Constant(1)),Control(Constant(1))]
				#		temp = call_sens_analysis(to_flux[k_step]*u_final[k_step],controls,dP(1),all_together)
				#		sensitivity_output[k_step].append(temp)
				#	simulation_time_list.append(time.time() - time_size)
			sens_time_list.append(t)

			progressBar(t, reac_input['Pulse Duration'])
			u_n.assign(u)
			t += dt

		if reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['RRM Analysis'].lower() == 'true':
			print()
			print()
			print('Evaluating Tape with Tangent Linear Method. Could take some time.')
			tape.evaluate_tlm()

		if reac_input['Sensitivity Analysis'].lower() == 'true':
			
			for numEachSens,eachSens in enumerate(u_graph_data):
				np.savetxt(sensFolder+'/c_'+legend_label[numEachSens]+'.csv',u_graph_data[eachSens],delimiter=",")#

			for numEachSens,eachSens in enumerate(sensFuncs):
				newList = []
				for kSensNum, kSens in enumerate(sensFuncs[eachSens]):
					newValue = kSens.block_variable.tlm_value
					newList.append(newValue)
				np.savetxt(sensFolder+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")#
			
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
			
			#sys.exit()
		
		current_reac_set += 1
		if current_reac_set >  pulse_variation:
			current_reac_set = 1

		x_values = []
		it_times = []
		j_values = []
		dj_values = []

		if reac_input['Thin-Zone Analysis'].lower() == 'true':
			if int(reac_input['Number of Pulses']) == 1:
				for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					#print(necessary_values['reactants'][j_species])
					np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
					rateTest = eval(rateStrings[j_species])
					np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")
			else:
				pulse_path = './'+reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(k_pulse+1)+'/'
				generate_folder(pulse_path)
				for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					#print(necessary_values['reactants'][j_species])
					np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")
					rateTest = eval(rateStrings[j_species])
					np.savetxt(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv', np.array(rateTest), delimiter=",")

		def eval_cb(j, m):
			print('eval')
			print(j)
			print(m)

		def deriv_cb(j,dj,m):
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
		
		if fen_17 == True:
			set_log_active(False)
		fitting_time = time.time()
				
		if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true':
			#if reac_input['Fit Parameters'].lower() == 'true':

			rf_2 = ReducedFunctional(jfunc_2, controls,derivative_cb_post=deriv_cb)
			#else:
			#	rf_2 = ReducedFunctional(jfunc_2, controls,derivative_cb_post=deriv_inert)
			low_bounds = []
			up_bounds = []
			#try:
			for gt in range(0,len(controls)):
				low_bounds.append(0)
				up_bounds.append(np.inf)
			if reac_input['Optimization Method'] == 'L-BFGS-B' or reac_input['Optimization Method'] == '':
				u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds),tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
			elif reac_input['Optimization Method'] == 'Newton-CG':
				u_opt_2 = minimize(rf_2, method = 'Newton-CG',tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
			elif reac_input['Optimization Method'] == 'BFGS':
				u_opt_2 = minimize(rf_2, method = 'BFGS',tol=1e-13, options={"gtol":1e-13})
			elif reac_input['Optimization Method'] == 'SLSQP':
				u_opt_2 = minimize(rf_2, method = 'SLSQP', bounds = (low_bounds,up_bounds),tol=1e-9, options={"ftol":1e-9})
			elif reac_input['Optimization Method'] == 'CG':
				u_opt_2 = minimize(rf_2,bounds = (low_bounds,up_bounds), method = 'CG',tol=1e-9, options={"gtol":1e-9})
			elif reac_input['Optimization Method'] == 'basinhopping':
				u_opt_2 = minimize(rf_2, method = 'basinhopping', bounds = (low_bounds,up_bounds),tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
			
			else:
				print('Requested Optimization Method Does Not Exist')
				sys.exit()
			optimization_success = True
			#except RuntimeError:
			#	print('Optimization Failed or Entered Too Stiff Region')
			#	time.sleep(1.5)
			#	print('Will shortly begin generating the optimization gif.')
			#	time.sleep(5)
			#	optimization_success = False
		if reac_input['Uncertainty Quantification'].lower() == 'true':
			h = Constant(1)
			print(type(jfunc_2))
			#Jfunc_2.adj_value = 1.0
			c.tlm_value = h
			tape.evaluate_adj()
			tape.evaluate_tlm()
			Jfunc_2.block_variable.hessian_value = 0
			tape.evaluate_hessian()
			Hm = c.originial_block_variable.hessian_value
			print(Hm)
			np.savetxt(sensFolder+'/k1.csv',Hm,delimiter=",")#
			sys.exit()



	############# Store the output data #######################################################
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

	############# Visualize/graph the data #######################################################
		if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Display Experimental Data'].lower() == 'true':

			#if reac_input['Fit Parameters'].lower() == 'true':
			if reac_input['Display Objective Points'].lower() == 'true' or reac_input['Fit Parameters'].lower() == 'true':
				for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
					if objSpecies[k_fitting] == '1':
						ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting])
		
			#if reac_input['Fit Parameters'].lower() == 'true':
			if reac_input['Display Objective Points'].lower() == 'true' and reac_input['Fit Inert'].lower() == 'true':
				for k_fitting in range( int(len(legend_label)-reac_input['Number of Inerts']) ,len(legend_label[int(len(legend_label)+reac_input['Number of Inerts'])])):
					ax2.scatter(output_fitting[legend_label[monitored_gas+k_fitting]]['times'],output_fitting[legend_label[monitored_gas+k_fitting]]['values'],marker='^',color='r',label='Fitting'+legend_label[monitored_gas+k_fitting], alpha=0.3)



			for k,j in enumerate(graph_data):
				if j != 'timing': #and j != 'conVtime_1' and j != 'conVtime_2' and j != 'conVtime_3' and j != 'conVtime_0' 

					if k_pulse > 0:
						pass
					else:
						dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/flux_data/'+legend_label[k]+'.csv',header=None)
						ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
				else:
					pass




		
		if reac_input['Infinite Inert'].lower() == 'true':
			
			for kjc in range(0,monitored_gas):

				outlet = []
				outlet.append(0)
				
				zAnalytical = 1
				while zAnalytical*dt < reactant_time[kjc]:#!!!!!!!!!!!!!!Need the actual time at this point
					#print(zAnalytical*dt) 
					outlet.append(0)
					zAnalytical+=1	
				outlet.append(0)
				
				if reac_input['Scale Output'].lower() == 'true':
					factor = 1
				else:
					factor = float(species_pulse_list[kjc])*reac_input['Reference Pulse Size']

				for k in range(zAnalytical+1,int(reac_input['Time Steps'])+1):
					analyticalValue = 0
					for n in range(0,50):
						analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - reactant_time[kjc])*(D[kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
						#outlet.append((Dout[3]*3.14159/(eb[0]*(np.sum(r_param)**2)))*nsum(lambda x: ((-1)**x)*(2*x+1)*exp((-(x+0.5)**2)*(3.14159**2)*((k*dt)*(Dout[3]/(eb[0]*(np.sum(r_param)**2))))) ,[0, inf]))
					outlet.append(factor*(D[kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
				#print(dfNew[0][:].tolist())
				#print(outlet)
				#sys.exit()
				try:
					ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
				except ValueError:
					outlet = outlet[:-1]
					ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)


			#for kjc in range(all_molecules-int(reac_input['Number of Inerts']),all_molecules):
			for kjc in range(0,int(reac_input['Number of Inerts'])):
				outlet = []
				outlet.append(0)
				
				zAnalytical = 1

				while zAnalytical*dt < reactant_time[monitored_gas+kjc]:#!!!!!!!!!!!!!!Need the actual time at this point
					#print(zAnalytical*dt) 
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
						analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - reactant_time[monitored_gas+kjc])*(D[monitored_gas+kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
						#outlet.append((Dout[3]*3.14159/(eb[0]*(np.sum(r_param)**2)))*nsum(lambda x: ((-1)**x)*(2*x+1)*exp((-(x+0.5)**2)*(3.14159**2)*((k*dt)*(Dout[3]/(eb[0]*(np.sum(r_param)**2))))) ,[0, inf]))
					outlet.append(factor*(D[monitored_gas+kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
				#print(dfNew[0][:].tolist())
				#print(outlet)
				#sys.exit()
				try:
					ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
				except ValueError:
					outlet = outlet[:-1]
					ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)

			###outlet = []
			###outlet.append(0)
		
		
			###zAnalytical = 1
			###while zAnalytical*dt < 0.01:
			###	#print(zAnalytical*dt) 
			###	outlet.append(0)
			###	zAnalytical+=1	
			###outlet.append(0)
			###for k in range(zAnalytical+1,int(reac_input['Time Steps'])+1):
			###	analyticalValue = 0
			###	#print(k*dt-0.01)
			###	#time.sleep(1)
			###	for n in range(0,50):
			###		analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*(((k*dt-0.01))*(D[3][0]/(eb[0]*(np.sum(r_param)**2)))))
			###		#outlet.append((Dout[3]*3.14159/(eb[0]*(np.sum(r_param)**2)))*nsum(lambda x: ((-1)**x)*(2*x+1)*exp((-(x+0.5)**2)*(3.14159**2)*((k*dt)*(Dout[3]/(eb[0]*(np.sum(r_param)**2))))) ,[0, inf]))
			###	outlet.append((D[3][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
			####print(dfNew[0][:].tolist())
			####print(outlet)
			####sys.exit()
			###ax2.plot(graph_data['timing'],outlet,color='r',label='Analytical Inert-1', alpha=0.7)
				
			###outlet = []
			###outlet.append(0)
			####print(np.sum(r_param))
			####sys.exit()
			###for k in range(1,int(reac_input['Time Steps'])+1):
			###	analyticalValue = 0
			###	for n in range(0,50):
			###		analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt)*(D[4][0]/(eb[0]*(np.sum(r_param)**2)))))
			###		#outlet.append((Dout[3]*3.14159/(eb[0]*(np.sum(r_param)**2)))*nsum(lambda x: ((-1)**x)*(2*x+1)*exp((-(x+0.5)**2)*(3.14159**2)*((k*dt)*(Dout[3]/(eb[0]*(np.sum(r_param)**2))))) ,[0, inf]))
			###	outlet.append((D[4][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
			####print(dfNew[0][:].tolist())
			####print(outlet)
			####sys.exit()
			###ax2.plot(graph_data['timing'],outlet,color='k',label='Analytical Inert-1', alpha=0.7)

		sig = 0.05
		beta_2 = 0.00270
		w_2 = 2*3.14159*70

		for k,j in enumerate(graph_data):
			if j != 'timing': #and j != 'conVtime_1' and j != 'conVtime_2' and j != 'conVtime_3' and j != 'conVtime_0' 
				if reac_input['Noise'].lower() == 'true':
					
					for z in range(0,int(reac_input['Time Steps'])):
						graph_data[j][z] += np.random.normal(0,1)*sig +beta_2*np.cos(w_2*(k*dt))
				if k_pulse > 0:
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
				else:
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)#/reac_input['Reference Pulse Size']*(reac_input['Number of Reactants']+1)
				new_data = np.asarray(graph_data[j])
				dictionary_of_numpy_data[legend_label[k]] = np.vstack((dictionary_of_numpy_data[legend_label[k]],new_data))
			else:
				pass
		
#		if reac_input['Sensitivity Analysis'].lower() == 0:#'true':
#			for k_sens_step in range(monitored_gas):
#				sens_time = np.asarray(sens_time_list)
#				sens_time = sens_time.T#np.transpose(sens_time)
#				sensitivity_output_2 = np.asarray(sensitivity_output[k_sens_step])
#				sensitivity_output_2 = np.append(sens_time[:,None],sensitivity_output_2,axis=1)	
#				np.savetxt(reac_input['Output Folder Name']+'_folder/sensitivity/'+legend_label[k_sens_step]+'/pulse_'+str(k_pulse+1)+'.csv',sensitivity_output_2,delimiter=",",header='t,'+','.join(legend_2))#
#			df_sens_time = pd.DataFrame(simulation_time_list)
#			df_sens_time.to_csv(reac_input['Output Folder Name']+'_folder/sensitivity/time.csv',header=None)			
#
#		if reac_input['RRM Analysis'].lower() == 'true':
#			for k_sens_step in range(len(necessary_values['reactants'])):#necessary_values['reactants']
#				RRM_time = np.asarray(RRM_time_list)
#				#sens_time = np.asarray(graph_data['timing'][0:])
#				RRM_time = RRM_time.T#np.transpose(sens_time)
#				RRM_der_2 = np.asarray(RRM_der[k_sens_step])
#				RRM_der_2 = np.append(RRM_time[:,None],RRM_der_2,axis=1)	
#				np.savetxt(reac_input['Output Folder Name']+'_folder/RRM_derivatives/'+necessary_values['reactants'][k_sens_step]+'/pulse_'+str(k_pulse+1)+'.csv',RRM_der_2,delimiter=",",header='t,'+','.join(legend_2))#
#			df_sens_time = pd.DataFrame(simulation_time_list)
#			df_sens_time.to_csv(reac_input['Output Folder Name']+'_folder/RRM_derivatives/time.csv',header=None)			

		name_list = necessary_values['reactants'][monitored_gas:]

	ax2.legend(title="Gas Species")#!#!#!#!#
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	
	if reac_input['Store Graph'].lower() == 'true':
		plt.savefig('./'+reac_input['Output Folder Name']+'_folder/graphs/flux_data.png')
	
	if reac_input['Display Graph'].lower() == 'true':
		plt.show()
	
	for j_species in range(0,monitored_gas+int(reac_input['Number of Inerts'])):
		dictionary_of_numpy_data[legend_label[j_species]] = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
		np.savetxt('./'+reac_input['Output Folder Name']+'_folder/flux_data/'+legend_label[j_species]+'.csv', dictionary_of_numpy_data[legend_label[j_species]], delimiter=",")
	
	
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
			#if optimization_success == False:
			#things -= 1

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
				#print(alter)
				#sys.exit()
				alter.to_csv('./input_file.csv',header=None,index=False)	
			
				try:
					call_sim()
					sys.exit()
				except:
					k_num = things
		
			generate_gif(legend_label[:len(legend_label)], reac_input['Experimental Data Folder']+'/flux_data', './'+reac_input['Output Folder Name']+'_folder/fitting', len(constants), constants, reactor_kinetics_input['reactions_test'], times)
		
			for k_num in range(0,things):
				shutil.rmtree('./'+reac_input['Output Folder Name']+'_folder/fitting/iter_'+str(k_num)+'_folder') 

		user_data = pd.read_csv('./'+reac_input['Output Folder Name']+'_folder/input_file.csv',header=None)
		user_data.to_csv('./input_file.csv',header=None,index=False)

	return graph_data, legend_label, necessary_values['reactants']

def call_sim():
	
	reactor_kinetics_input,kinetic_parameters,kin_in = read_input()

	if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' or reactor_kinetics_input['RRM Analysis'].lower() == 'true' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
		for parameters in kinetic_parameters:
			reactor_kinetics_input,kinetic_parameters,kin_in = read_input()
			reactor_kinetics_input['Sensitivity Parameter'] = parameters
			reactor_kinetics_input['Display Graph'] = 'FALSE'
			if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or reactor_kinetics_input['Fit Inert'].lower() == 'true':
				print('')
				print('')
				
				print('Running the Sensitivity/RRM Analysis and Parameter Fitting methods simultaniously is not possible due to conflicts between the tangent linear and adjoint methods.')
				print('Please run again with one of these methods excluded.')
				sys.exit()

			graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)
	else:
		graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)

call_sim()
