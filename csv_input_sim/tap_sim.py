from fenics import *
from fenics_adjoint import *
#import pyipopt
from funs_tap_sim import *
from variational_form_constructor import make_f_equation
from read_exp import exp_data_fitting_3
from mpmath import nsum, exp, inf
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import csv
import sys
import os
import ufl

"""
The core function behind the tap simulator. When the 'user_input.py' file is run, the user defined data
from the input file will be parsed and fed to this function. 
"""

def tap_simulation_function(reactor_kinetics_input,constants_input):


	simulation_time_list = []
	#######################################################################################
	############# Read in the data and redefine information dictionary names ##############
	#######################################################################################

	reac_input = reactor_kinetics_input
	r_const = constants_input
	path = './'+reac_input['output_file_name']+'_folder/'
	try:  
		os.mkdir(path)
	except OSError:  
		print ("Creation of the directory %s failed" % path)
		pass
	else:  
		print ("Successfully created the directory %s " % path)
	
	path_3 = './'+reac_input['output_file_name']+'_folder/flux_data/'
	try:  
		os.mkdir(path_3)
	except OSError:  
		print ("Creation of the directory %s failed" % path_3)
		pass
	else:  
		print ("Successfully created the directory %s " % path_3)

	path_4 = './'+reac_input['output_file_name']+'_folder/thin_data/'
	try:  
		os.mkdir(path_4)
	except OSError:  
		print ("Creation of the directory %s failed" % path_4)
		pass
	else:  
		print ("Successfully created the directory %s " % path_4)

	path_5 = './'+reac_input['output_file_name']+'_folder/graphs/'
	try:  
		os.mkdir(path_5)
	except OSError:  
		print ("Creation of the directory %s failed" % path_5)
		pass
	else:  
		print ("Successfully created the directory %s " % path_5)

	for j in r_const:
		r_const[j] = Constant(r_const[j])
	

	#new_info
	refs = np.loadtxt("recorded.txt")
	gamma = 1.e-5
	if gamma > 0:
		noise = np.random.normal(0, gamma, refs.shape[0])

	# add noise to the refs
	refs += noise
	refs = list(map(Constant, refs))

	controls = []
	legend_2 = []
	for j in r_const:
		controls.append(Control(r_const[j]))
		legend_2.append(j)

	user_data = pd.read_csv('./input_file.csv',header=None)
	user_data.to_csv(path+'input_file.csv',header=None)
	
	##########################################################################
	############# Used to control the output/warnings of Fenics ##############
	##########################################################################

	parameters["std_out_all_processes"] = False																							
	cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
	set_log_active(False)
	import warnings
	warnings.filterwarnings("ignore", category=DeprecationWarning)
	tol = 1E-14
	theta = reac_input['theta']

	if ',' in str(reac_input['rangesurface_species']):
		rangesurface_species = list(reversed(reac_input['rangesurface_species'].split(',')))
		rangesurface_species = rangesurface_species[::-1]
	else:
		rangesurface_species = reac_input['rangesurface_species']

	##########################################################################
	############# Establish transport properties in FEniCS ###################
	##########################################################################

	def vect_vol(x,y):
		return x*(reac_input['reac_radius']**(2))*np.pi*y
	r_param, dx_r, dx2_r, frac_length, cat_location = establish_grid_system(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['mesh_size'])

	original_pulse_size = reac_input['Inert_pulse_size']*(reac_input['reactants_num']+1)#1#
	dk = Constant(reac_input['Time']/reac_input['Time_steps'])
	eb = np.array((reac_input['void_inert'],reac_input['void_cat'],reac_input['void_inert']))
	
	ref_rate = np.append(reac_input['ref_rate_inert'],reac_input['ref_rate_cat'])
	ref_rate = np.append(ref_rate,ref_rate[0])  	
	D = np.empty((len(reac_input['mass_list'].split(',')),3))#,dtype=Constant
	def diff_func(ref_m,mol_mass,ref_r):     	
		return ref_r*(mp.sqrt(40*reac_input['T'])/mp.sqrt(423.15*mol_mass))
	#Constant(D[j_num][0]/(np.sum(r_param)**2)))
	for k,j in enumerate(reac_input['mass_list'].split(',')):
		for k_2,j_2 in enumerate(ref_rate):
			#D[k,k_2] = Constant(float(j))
			D[k,k_2] = Constant(diff_func(reac_input['ref_mass'],float(j),j_2)/(np.sum(r_param)**2))
			#D[k,k_2] = Constant(D[k,k_2])
			#print(type(D[0,0]))
	#D = Constant(D)
	#print(type(D[0,0]))
	#print(type(D))
	#sys.exit()
	#print(D)

	##########################################################################
	############# Establish reactor parameters in FEniCS #####################
	##########################################################################


	ca = (reac_input['reac_radius']**2)*3.14159				
	point_volume = dx_r * ca 
	#print(dx_r)
	#print(dx_r*ca)
	#sys.exit()
	specific_area = 1e8
	surface_area_at_point = point_volume*specific_area
	sdensity =1
	sarea=1
	open_sites_per_point = sdensity*surface_area_at_point
	
	vect_vol_imp = np.vectorize(vect_vol)
	vr = np.sum(vect_vol_imp(r_param,eb))		### Volume of the reactor
	vc = r_param[1]*(reac_input['reac_radius']**(2))*np.pi*eb[1]		### Volume of the catalyst zone
	Q = sarea/vc 			### surface to gas scaling coefficient (1/cm)
	Inert_pulse_size = reac_input['Inert_pulse_size']/(dx_r*ca)
	
	necessary_values = make_f_equation(reac_input['reactions_test'],reac_input['reactants_num'],reac_input['reactor_type'],False)
	
	############# Establish the mesh and boundaries/domains in FEniCS #######################################################
	
	mesh = UnitIntervalMesh(int(reac_input['mesh_size']))
	P1 = FiniteElement('P',mesh.ufl_cell(),1)
	
	if reac_input['reactions_test'] != ['INERT_ONLY']:
		test_new = eval(necessary_values['element'])
		element = MixedElement(test_new)
		V = FunctionSpace(mesh,element)
		V_du = FunctionSpace(mesh,P1)
	else:
		V = FunctionSpace(mesh,P1)
		V_du = FunctionSpace(mesh,P1)
	
	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], ((cat_location - 0.5*frac_length), (cat_location + 0.5*frac_length)))

	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	thin_zone.mark(domains,1)
	
	dx = Measure("dx")[domains]
	
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
	
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)
	
	vcu = TestFunction(V)
	
	dz = 1/reac_input['mesh_size']

	class integration_section(SubDomain):
			def inside(self, x, on_boundary):
				return between(x[0], (1-dz,1.0))

	right = CompiledSubDomain("near(x[0], 1.)")

	boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
	#left.mark(boundary_parts, 0)
	right.mark(boundary_parts, 1)
	ds = Measure("ds", subdomain_data=boundary_parts)

	all_molecules = necessary_values['gas_num']+necessary_values['surf_num']+1
	u = Function(V)
	u_n = Function(V)
	u_temp = Function(V)
	u_1 = Function(V)
	u_2 = Function(V)
	u_3 = Function(V)
	
	graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data = initialize_variable_dictionaries(necessary_values,all_molecules,V,u,u_n)
	
	dt = reac_input['Time']/reac_input['Time_steps']
	Length = np.sum(r_param)**2
	
	try:
		F = eval(necessary_values['F'])
	except NameError:
		error_output(reac_input['reactions_test'])
	#### Just an equations for reeces code (different from others since it requires the use of a scaling parameter ( Q ))
	Q = Constant(5.5423e3)

	#F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0) + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1) + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0) + dk*Constant(D[1][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1) + dk*Constant(D[1][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + ((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + ((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + ((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + ((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + ((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + ((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1)+ ((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + ((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1)	- dk*(1.0* Q*Kd0*(u_d['u_3']**1.0)*v_d['v_1']*dx(1)) + dk*(1.0* Ke0*(u_d['u_1']**1.0)*v_d['v_1']*dx(1)) + dk*(1.0* Kd0*(u_d['u_3']**1.0)*v_d['v_3']*dx(1)) - dk*(1.0* (Ke0/Q)*(u_d['u_1']**1.0)*v_d['v_3']*dx(1)) + dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_3']*dx(1)) + dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_4']*dx(1)) - dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_5']*dx(1)) + dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_6']*dx(1))- dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke3*(u_d['u_5']**1.0)*v_d['v_5']*dx(1))	- dk*(1.0* Q*Ke3*(u_d['u_5']**1.0)*v_d['v_2']*dx(1))	+ Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)	+ dk * Constant(D[2][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0) + Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1) + dk * Constant(D[2][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)

	#F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)+ dk*D[0][0]*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)+ dk*D[0][1]*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1)	+ Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0) + dk*D[1][0]*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1) + dk*D[1][1]*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + ((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + ((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + ((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + ((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + ((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + ((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1) + ((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + ((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1) - dk*(1.0* Q*r_const['Kd0']*(u_d['u_3']**1.0)*v_d['v_1']*dx(1)) + dk*(1.0* r_const['Ke0']*(u_d['u_1']**1.0)*v_d['v_1']*dx(1)) + dk*(1.0* r_const['Kd0']*(u_d['u_3']**1.0)*v_d['v_3']*dx(1)) - dk*(1.0* (r_const['Ke0']/Q)*(u_d['u_1']**1.0)*v_d['v_3']*dx(1)) + dk*(1.0* r_const['Ke1']*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_3']*dx(1)) + dk*(1.0* r_const['Ke1']*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_4']*dx(1)) - dk*(1.0* r_const['Ke1']*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_5']*dx(1)) + dk*(1.0* r_const['Ke2']*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_3']*dx(1)) + dk*(1.0* r_const['Ke2']*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_6']*dx(1)) - dk*(1.0* r_const['Ke2']*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_5']*dx(1)) + dk*(1.0* r_const['Ke3']*(u_d['u_5']**1.0)*v_d['v_5']*dx(1)) - dk*(1.0* Q*r_const['Ke3']*(u_d['u_5']**1.0)*v_d['v_2']*dx(1)) + Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)+ dk * D[2][0] * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0)	+ Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1)	+ dk * D[2][1] * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)
	
	monitored_gas = necessary_values['molecules_in_gas_phase']
	
	bcs = define_boundary_conditions(reac_input['reactor_type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['reactants_num'],all_molecules,reac_input['reactant_ratio'],boundary_L,boundary_R)
	
	############# Initialize the graphs and set the graphical parameters #######################################################
	
	fig2,ax2,legend_label,header,colors = establish_output_graph(reac_input['reactor_type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'])
	

	if reac_input['fit_parameters'].lower() == 'true':
		output_fitting = exp_data_fitting_3(legend_label,reac_input['Time_steps'],reac_input['exp_data_folder'])

	if reac_input['sensitivity_analysis'].lower() == 'true':
		path_2 = './'+reac_input['output_file_name']+'_folder/sensitivity_'+reac_input['output_file_name']+'/'
		for k_gas in legend_label[:-1]:
			try:  
				os.mkdir(path_2+k_gas+'/')
			except OSError:  
				print ("Creation of the directory %s failed" % path_2+k_gas+'/')
				pass
			else:  
				print ("Successfully created the directory %s " % path_2+k_gas+'/')
		
		for j in r_const:
			r_const[j] = Constant(r_const[j])


	to_flux = flux_generation(reac_input['reactor_type'],monitored_gas,reac_input['reactants_num'],original_pulse_size,D,eb,dx_r,reac_input['reac_radius'],dx2_r)
	store_data_func(reac_input['store_data'],reac_input['output_file_name'])
	store_sens_analysis(reac_input['sensitivity_analysis'],reac_input['output_file_name'],monitored_gas,legend_label)

	############# Declare the solver to be used #######################################################
	
	class Plane(SubDomain):
		def inside(self, x, on_boundary):
			return x[0] > 1.0 - DOLFIN_EPS

	J = derivative(F,u)
	problem = NonlinearVariationalProblem(F,u,bcs,J)
	solver = NonlinearVariationalSolver(problem)
	
	solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
	#solver.parameters["newton_solver"]["absolute_tolerance"] = 1e3
	solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10
	#solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-18
	############# Define a method to iteratively choose time-step (will change based on meeting with tobin isaac) ##
	#%#%#%#%#% TASK #3	
	############# Solve the system #######################################################
	if '/' in str(reac_input['reactant_ratio']):
		list_of_feed = reac_input['reactant_ratio'].split('/')
		pulse_variation = len(list_of_feed)
		reactant_feed = []
		for k in range(0,len(reac_input['reactant_ratio'].split('/'))):
			reactant_feed.append(list_of_feed[k].split(','))

	else:
		pulse_variation = 1
		reactant_feed = []
		try:
			reactant_feed.append(reac_input['reactant_ratio'].split(','))
		except AttributeError:
			reactant_feed.append(reac_input['reactant_ratio'])

	current_reac_set = 1

	for k_pulse in range(0,int(reac_input['number_of_pulses'])):
		tape.clear_tape()
		
		print("")
		print("Simulation Status: "+"Pulse #"+str(k_pulse+1))

		species_pulse_list = reactant_feed[current_reac_set-1]
		
		dk = Constant(reac_input['Time']/reac_input['Time_steps'])
		dt = reac_input['Time']/reac_input['Time_steps']
		
		start_time = time.time()
		
		sensitivity_output = {}
		for k_sens in range(monitored_gas):
			sensitivity_output[k_sens] = []
		
		graph_data['timing'] = []
		for k_gasses in range(0,necessary_values['molecules_in_gas_phase']):
			graph_data['conVtime_'+str(k_gasses)] = []
		graph_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
		
		surf_data['timing'] = []
		for j_gasses in range(necessary_values['molecules_in_gas_phase'],len(necessary_values['reactants'])-1):
			surf_data['conVtime_'+str(j_gasses)] = []

		cat_data['timing'] = []

		for z_gasses in range(0,necessary_values['molecules_in_gas_phase']):
			
			cat_data['conVtime_'+str(z_gasses)] = []

		t = 0
		if reac_input['fit_parameters'].lower() == 'true':
			osub = integration_section()
			domains = CellFunction("size_t", mesh)
			domains.set_all(0)
			osub.mark(domains, 1)
			dP = Measure('vertex',domain = mesh, subdomain_data=domains)
		
		test_new = project(-u.dx(0),V)
		test_new_split = test_new.split(deepcopy=True)

		w_new = Expression("1", degree=0)
		w_new2 = interpolate(w_new,V_du)


		W = VectorFunctionSpace(mesh, 'P', 1)
		w = Function(W)
		w.vector()[:] = 10

		x_dim = list(range(0, int(reac_input['mesh_size'])+1))
		
		cum_sum =0
		while t <= reac_input['Time']:
			#####STORE THE RESULTS FOR EACH ITERATION
			graph_data['timing'].append(t)
			for k in range(0,monitored_gas):
				new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k+1]))#/(2*Inert_pulse_size)
				graph_data['conVtime_'+str(k)].append((new_val))
			new_val = ((to_flux[monitored_gas]*u.vector().get_local()[2*(all_molecules+1)-1]))#/(2*Inert_pulse_size)
			graph_data['conVtime_'+str(all_molecules-1)].append((new_val))
			
			new_values = [] 
			cat_values = []

			mesh_size = reac_input['mesh_size']
			if reac_input['fit_parameters'].lower() == 'true':
				for k_fitting in range(0,len(legend_label[:len(legend_label)])):
					if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:

						c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
						slope = (-c_exp)/(1/mesh_size)
						intercept = c_exp - ((1-(1/mesh_size))*slope)
						w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)#Expression("1", degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)#output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]*
						test_meaning = assemble(inner(u_n[k_fitting] + output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))] ,u_n[k_fitting] + output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))])*ds(1))

						print(all_molecules)
						try:
							if legend_label[k_fitting] != 'Inert':
								jfunc += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
							else:
								jfunc += assemble(inner(u_n[all_molecules] - w3,u_n[all_molecules] - w3)*dP(1))

						except UnboundLocalError:
							if legend_label[k_fitting] != 'Inert':
								jfunc = assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
							else:
								jfunc = assemble(inner(u_n[all_molecules] - w3,u_n[all_molecules] - w3)*dP(1))
						#print(type(test_new_split[k]))
						#sys.exit()


			for cj in range(0,monitored_gas+1):
				cat_values.append(0)

			for kj in range(monitored_gas+1,len(necessary_values['reactants'])+1):
				new_values.append(0)
			
			for z in range(int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1,int((cat_location + 0.5*frac_length)*reac_input['mesh_size'])):
				#print(z)
				#print()
				try:
					for z_num in range(0,1+necessary_values['molecules_in_gas_phase']):#.split(',')
						#print(z*(all_molecules+1)-(z_num))
						value_cat = u_n.vector().get_local()[z*(all_molecules+1)-(z_num)+1]# = float(z_sur)
						cat_values[z_num] += value_cat*dx_r
						#surf_data['conVtime_'+str(z_num+monitored_gas)] += value_new*dx_r
				except AttributeError:
					value_cat = u_n.vector().get_local()[z*(all_molecules+1)]# = float(z_sur)
					cat_values[0] += value_cat*dx_r
				#print(z*(all_molecules+1)-all_molecules)
				#value_cat = u_n.vector().get_local()[z*(all_molecules+1)-all_molecules]
				#cat_values[necessary_values['molecules_in_gas_phase']] += value_cat*dx_r

				try:
					swap_for_2 = necessary_values['molecules_in_gas_phase']
					for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
						#print(z*(all_molecules+1)-(swap_for_2+z_num))
						value_new = u_n.vector().get_local()[z*(all_molecules+1)-(swap_for_2+z_num)]# = float(z_sur)
						new_values[z_num] += value_new*dx_r
						#surf_data['conVtime_'+str(z_num+monitored_gas)] += value_new*dx_r
				except AttributeError:
					value_new = u_n.vector().get_local()[z*(all_molecules+1)]# = float(z_sur)
					new_values[0] += value_new*dx_r
			#print(new_values)
			cat_length = reac_input['length_reac']*reac_input['factor']
			#if ',' in str(reac_input['rangesurface_species']):
			#	print(reac_input['rangesurface_species'].split(','))
			#	for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
			#		cat_data['conVtime_'+str(z_num)].append(float(z_sur)/(cat_length))
			#else:
			#	cat_data['conVtime_'+str(0)].append(cat_values[0])

			#if ',' in str(reac_input['rangesurface_species']):
			#	for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
			#		surf_data['conVtime_'+str(z_num+monitored_gas)].append(new_values[z_num])
			#else:
			#	surf_data['conVtime_'+str(monitored_gas)].append(new_values[0])
			#new_val = u.vector().get_local()[(all_molecules)+jk+1]
			#surf_data['conVtime_'+str(all_molecules-1)]
			#####Solve for remaining time steps
			
			#if ',' in str(reac_input['rangesurface_species']):
			#	rangesurface_species = list(reversed(reac_input['rangesurface_species'].split(',')))
			#else:
			#	rangesurface_species = reac_input['rangesurface_species']


			if t > 0:
				dt = solver_iteration(dt,reac_input['solver_method_options'],solver,dk,1.5,1.1)
			else:
				if reac_input['reactor_type'] == 'tap':
					if ',' in str(species_pulse_list):
						for k in range(0,int(reac_input['reactants_num'])):
							u_n.vector()[int((all_molecules+1)*(reac_input['mesh_size']+1)-1)+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_size#/(point_volume)#list_species_pulse[k]
					else:
						u_n.vector()[int((all_molecules+1)*(reac_input['mesh_size']+1)-1)-(all_molecules)] = float(species_pulse_list)*Inert_pulse_size#/(point_volume)
					u_n.vector()[int((all_molecules+1)*(reac_input['mesh_size']+1)-1)] = Inert_pulse_size#/(point_volume)
				
				if k_pulse == 0:
					for z in range(int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1,int((cat_location + 0.5*frac_length)*reac_input['mesh_size'])):
						if ',' in str(reac_input['rangesurface_species']):
							for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
								if int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1 <= 0:
									u_n.vector()[(all_molecules+1)-(2+z_num)] = float(z_sur)
								else:
									u_n.vector()[z*(all_molecules+1)-(2+z_num)] = float(z_sur)
						else:
							if int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1 <= 0:
								u_n.vector()[(all_molecules+1)-(2)] = float(species_pulse_list)
							else:
								u_n.vector()[z*(all_molecules+1)-(2)] = float(species_pulse_list)

				dt = solver_iteration(dt,reac_input['solver_method_options'],solver,dk,1.5,1.1)
				
			if reac_input['sensitivity_analysis'].lower() == 'true':
				time_size = time.time()
				u_final = u.split(deepcopy=False)
				for k_step in range(0,monitored_gas):
					temp = call_sens_analysis(u_final[k_step],controls,dP(1))
					sensitivity_output[k_step].append(temp)
				simulation_time_list.append(time.time() - time_size)
			progressBar(t, reac_input['Time'])
			u_n.assign(u)
			t += dt
		#print(cum_sum)
		#sys.exit()
		current_reac_set += 1
		if current_reac_set >  pulse_variation:
			current_reac_set = 1
		def print_fun(x):
			print(x)
			
		set_log_active(False)
		fitting_time = time.time()

		#controls_2 = []
		print(legend_label)
		#for j_num,j in enumerate(legend_label):
		#	controls_2.append(Control(D[j_num][0]))
		if reac_input['fit_parameters'].lower() == 'true':

			rf = ReducedFunctional(jfunc, controls[1:])
			#u_opt = minimize(rf, method="COBYLA", bounds = ((0,0),(np.inf,np.inf)), callback=print_fun,tol = 1e-5)#, options={"disp": True} # , method="COBYLA", bounds = ((0,0),(np.inf,np.inf), callback=print_fun,tol = 1e-5
			#u_opt = minimize(rf, method="TNC", callback=print_fun,tol = 1e-5, options={"disp": True})#,options={'gtol': 1e-05})#  , bounds=(0,1000)  #method="L-BFGS-B" ###, , constraints={'type':'ineq', 'fun': rf}
			u_opt = minimize(rf, callback=print_fun, bounds = ((0,0),(np.inf,np.inf)), tol=1e-12, options={"disp": True,'maxiter': 15000, "gtol": 1.0e-12, "ftol": 1e-100})
		print('fitting done in '+str(time.time() - fitting_time))
		print("")
		total_simulation_time = time.time() - start_time
		print("Complete in "+str(time.time() - start_time)+" seconds")
		
	############# Store the output data #######################################################
		if k_pulse == 0:
			dictionary_of_numpy_data = {}
			for j_species in range(0,monitored_gas+1):
				dictionary_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionary_of_numpy_data[legend_label[j_species]] = np.vstack((dictionary_of_numpy_data[legend_label[j_species]],new_data))
			
			for k_species in range(monitored_gas+1,len(necessary_values['reactants'])+1):
				dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.vstack((dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]],new_data))
			dictionary_of_gas_cat_data = {}
			for z_species in range(0,monitored_gas):
				dictionary_of_gas_cat_data[legend_label[z_species]] = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionary_of_gas_cat_data[legend_label[z_species]] = np.vstack((dictionary_of_gas_cat_data[legend_label[z_species]],new_data))

#for z in range(int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1,int((cat_location + 0.5*frac_length)*reac_input['mesh_size'])):
#	for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
#		value_new = u_n.vector().get_local()[z*(all_molecules+1)-(2+z_num)]# = float(z_sur)
#		surf_data['conVtime_'+str(z_num+monitored_gas)] += value_new*dx_r

	############# Visualize/graph the data #######################################################
	
		for k,j in enumerate(graph_data):
			if j != 'timing':
				if k_pulse > 0:
					#new_set = [x / (reac_input['Inert_pulse_size']*(reac_input['reactants_num']+1)) for x in graph_data[j] ]
					#ax2.plot(graph_data['timing'],new_set,color=colors[k], ls = '--', alpha=0.7)
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
				else:
					#new_set = [x / (reac_input['Inert_pulse_size']*(reac_input['reactants_num']+1)) for x in graph_data[j] ]
					#print(type(graph_data[j]))
					#ax2.plot(graph_data['timing'],new_set,color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)#/reac_input['Inert_pulse_size']*(reac_input['reactants_num']+1)
				new_data = np.asarray(graph_data[j])
				dictionary_of_numpy_data[legend_label[k]] = np.vstack((dictionary_of_numpy_data[legend_label[k]],new_data))
			else:
				pass
		
		if reac_input['sensitivity_analysis'].lower() == 'true':
			for k_sens_step in range(monitored_gas):
				sens_time = np.asarray(graph_data['timing'][0:])
				sens_time = sens_time.T#np.transpose(sens_time)
				sensitivity_output_2 = np.asarray(sensitivity_output[k_sens_step])
				sensitivity_output_2 = np.append(sens_time[:,None],sensitivity_output_2,axis=1)	
				np.savetxt('./'+reac_input['output_file_name']+'_folder/sensitivity_'+reac_input['output_file_name']+'/'+legend_label[k_sens_step]+'/pulse_'+str(k_pulse+1)+'.csv',sensitivity_output_2,delimiter=",",header='t,'+','.join(legend_2))#
			df_sens_time = pd.DataFrame(simulation_time_list)
			df_sens_time.to_csv('./'+reac_input['output_file_name']+'_folder/sensitivity_'+reac_input['output_file_name']+'/time.csv',header=None)			

		name_list = necessary_values['reactants'][monitored_gas:]
		#for k,j in enumerate(cat_data):
		#	#print(surf_data)
		#	if j != 'timing':
		#		#print(necessary_values['reactants'][k+monitored_gas])
		#		#print('conVtime_'+str(k+monitored_gas))
		#		#print(j)#necessary_values['reactants'][k+monitored_gas]
#
		#		new_data = np.asarray(cat_data['conVtime_'+str(k)])
		#		#print(type(new_data))
		#		#print(dictionary_of_gas_cat_data)
		#		#sys.exit()
		#		dictionary_of_gas_cat_data[legend_label[k]] = np.vstack((dictionary_of_gas_cat_data[legend_label[k]][0],new_data))
		#	else:
		#		pass

		#for k,j in enumerate(surf_data):
		#	#print(surf_data)
		#	if j != 'timing':
		#		#print(necessary_values['reactants'][k+monitored_gas])
		#		#print('conVtime_'+str(k+monitored_gas))
		#		#print(j)#necessary_values['reactants'][k+monitored_gas]

		#		new_data = np.asarray(surf_data['conVtime_'+str(k+monitored_gas)])
		#		print(name_list[k])
		#		print(new_data)
				#sys.exit()
		#		dictionary_of_numpy_data[name_list[k-1]] = np.vstack((dictionary_of_numpy_data[name_list[k-1]][0],new_data))
		#	else:
		#		pass

	#reece = pd.read_csv('./reece_output_data.csv',header=None)
	#ax2.plot(reece[0],reece[1])
	#ax2.plot(reece[0],reece[2])
	#ax2.plot(reece[0],reece[3])
	ax2.legend(title="Gas Species")
	#ax2.set_ylim(0,0)
	if reac_input['save_figure'].lower() == 'true':
		plt.savefig('./'+reac_input['output_file_name']+'_folder/graphs/'+reac_input['output_file_name']+'.png')
	
	if reac_input['Display_figure'].lower() == 'true':
		plt.show()
	
	for j_species in range(0,monitored_gas+1):
		dictionary_of_numpy_data[legend_label[j_species]] = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
		np.savetxt('./'+reac_input['output_file_name']+'_folder/flux_data/'+legend_label[j_species]+'.csv', dictionary_of_numpy_data[legend_label[j_species]], delimiter=",")
	
	#for j_species in range(0,monitored_gas):
	#	dictionary_of_gas_cat_data[legend_label[j_species]] = np.transpose(dictionary_of_gas_cat_data[legend_label[j_species]])
	#	np.savetxt('./'+reac_input['output_file_name']+'_folder/thin_data/'+legend_label[j_species]+'.csv', dictionary_of_gas_cat_data[legend_label[j_species]], delimiter=",")

	#for k_species in range(monitored_gas+1,len(necessary_values['reactants'])+1):
	#	dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.transpose(dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]])
	#	np.savetxt('./'+reac_input['output_file_name']+'_folder/thin_data/'+necessary_values['reactants'][k_species-1]+'.csv', dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]], delimiter=",")

	for k,j in enumerate(graph_data):
		if j != 'timing':
			ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
		else:
			pass

	plt.clf()
	plt.close()
	#sys.exit()
	return total_simulation_time, graph_data, legend_label, necessary_values['reactants']