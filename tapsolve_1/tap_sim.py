from fenics import *
from fenics_adjoint import *
from func_sim import *
from vari_form import *
from reac_odes import *
from mpmath import nsum, exp, inf
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

#exp_data_fitting_2

def tap_simulation_function(reactor_kinetics_input,constants_input):

	### Initial forms for storing results and reading in the information ###
	simulation_time_list = []
	float_constants = constants_input.copy()
	reac_input = reactor_kinetics_input
	r_const = constants_input
	
	### Generate the folders for information storage ###
	path = './'+reac_input['output_folder_name']+'_folder/'
	generate_folder(path)
		
	path_3 = './'+reac_input['output_folder_name']+'_folder/flux_data/'
	generate_folder(path_3)
	
	if reac_input['mkm_analysis'].lower() == 'true' or reac_input['RRM_analysis'].lower() == 'true' or reac_input['petal_plots'].lower() == 'true':
		path_4 = './'+reac_input['output_folder_name']+'_folder/thin_data/'
		generate_folder(path_4)

	path_5 = './'+reac_input['output_folder_name']+'_folder/graphs/'
	generate_folder(path_5)

	if reac_input['fit_parameters'].lower() == 'true':
		path_6 = './'+reac_input['output_folder_name']+'_folder/fitting/'
		generate_folder(path_6)

	### Declare and define the constants ###
	for j in r_const:
		r_const[j] = Constant(r_const[j])

	### Define controls only if needed for differentiation based analysis ###
	if reac_input['fit_parameters'].lower() == 'true' or reac_input['sensitivity_analysis'].lower() == 'true' or reac_input['ad_rrm'].lower() == 'true':
		controls = []
		legend_2 = []
		for j in r_const:
			controls.append(Control(r_const[j]))
			legend_2.append(j)

	### Copy and store the input_csv file in the new folder (to keep organized) ###
	user_data = pd.read_csv('./input_file.csv',header=None)
	user_data.to_csv(path+'input_file.csv',header=None)
	
	### Control the output information from FEniCS ###
	parameters["std_out_all_processes"] = False
	cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}
	set_log_active(False)
	import warnings
	warnings.filterwarnings("ignore", category=DeprecationWarning)
	tol = 1E-14
	theta = reac_input['theta']

	### Used to define the initial surface composition over the simulation ###
	if ',' in str(reac_input['rangesurface_species']):
		rangesurface_species = list(reversed(reac_input['rangesurface_species'].split(',')))
		#rangesurface_species = rangesurface_species[::-1]
		reac_input['rangesurface_species'] = list(map(float, rangesurface_species))
	else:
		rangesurface_species = reac_input['rangesurface_species']
	
	### Initialize the grid system, time step size, pulse size and diffusion coefficients ###
	r_param, dx_r, dx2_r, frac_length, cat_location = establish_grid_system(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['mesh_size'])

	dk = Constant(reac_input['Time']/reac_input['Time_steps'])
	eb = np.array((reac_input['void_inert'],reac_input['void_cat'],reac_input['void_inert']))

	ref_rate = np.append(reac_input['ref_diff_inert'],reac_input['ref_diff_cat'])
	ref_rate = np.append(ref_rate,ref_rate[0])  	

	D = np.empty((len(reac_input['mass_list'].split(',')),3))#,dtype=Constant
	def diff_func(ref_mass,ref_T,mol_mass,ref_r):     	
		return ref_r*(mp.sqrt(ref_mass*reac_input['temp'])/mp.sqrt(ref_T*mol_mass))

	for k,j in enumerate(reac_input['mass_list'].split(',')):
		for k_2,j_2 in enumerate(ref_rate):
			D[k,k_2] = Constant(diff_func(reac_input['ref_mass'],reac_input['ref_temp'],float(j),j_2)) ###??? should this (np.sum(r_param)**2) be here? This puts it in dimensional form!

	### Define the dimensions of the reactor ###
	ca = (reac_input['reac_radius']**2)*3.14159 
	
	point_volume = dx_r * ca * eb[0] ###??? Should this include the voidage?
	#specific_area = 1e8
	#surface_area_at_point = point_volume*specific_area
	#sdensity =1
	#sarea=1
	#open_sites_per_point = sdensity*surface_area_at_point
	
	#def vect_vol(x,y):
	#	return x*(reac_input['reac_radius']**(2))*np.pi*y
	#vect_vol_imp = np.vectorize(vect_vol)
	#vr = np.sum(vect_vol_imp(r_param,eb))		### Volume of the reactor
	#vc = r_param[1]*(reac_input['reac_radius']**(2))*np.pi*eb[1]		### Volume of the catalyst zone
	#Q = sarea/vc 			### surface to gas scaling coefficient (1/cm)
	Inert_pulse_conc = reac_input['Inert_pulse_size']/(point_volume) ###??? Similar issue to the volume previously used

	dt = reac_input['Time']/reac_input['Time_steps']
	#Length = np.sum(r_param)**2

	### Construct the reaction equation in a form that Fenics can understand ###
	necessary_values = make_f_equation(reac_input['reactions_test'],reac_input['reactants_num'],reac_input['reactor_type'],reac_input['number_of_active_sites'],reac_input['number_of_inerts'],False)
	
	
	### Declaring the trial / test functions in fenics and defining the finite elements ###
	mesh = UnitIntervalMesh(int(reac_input['mesh_size']))  ###??? Swap the mesh size here
	P1 = FiniteElement('P',mesh.ufl_cell(),1)

	if reac_input['reactions_test'] != ['INERT_ONLY']:
		test_new = eval(necessary_values['element'])
		element = MixedElement(test_new)
		V = FunctionSpace(mesh,element)
		V_du = FunctionSpace(mesh,P1)
	else:
		V = FunctionSpace(mesh,P1)
		V_du = FunctionSpace(mesh,P1)
	#all_molecules = necessary_values['gas_num']+necessary_values['surf_num']+reac_input['number_of_inerts']
	all_molecules = necessary_values['gas_num']

	u = Function(V)
	u_n = Function(V)
	u_temp = Function(V)
	u_1 = Function(V)
	u_2 = Function(V)
	u_3 = Function(V)
	
	graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data = initialize_variable_dictionaries(necessary_values,all_molecules,V,u,u_n)
	
	### Define the subdomain for the catalyst region
	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], ((cat_location - 0.5*frac_length), (cat_location + 0.5*frac_length)))

	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	thin_zone.mark(domains,1)
	
	dx = Measure("dx")[domains]
	
	### Define the boundaries at each end of the reactor ###
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
	
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)
	
	dz = 1/reac_input['mesh_size']

	class integration_section(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (1-dz,1.0))
	
	right = CompiledSubDomain("near(x[0], 1.)")

	boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
	right.mark(boundary_parts, 1)
	#ds = Measure("ds", subdomain_data=boundary_parts)
	
	### Evaluate the reaction/diffusion expression for FEniCS to use ###
	try:
		F = eval(necessary_values['F'])
	except NameError:
		error_output(reac_input['reactions_test'])
	
	monitored_gas = necessary_values['molecules_in_gas_phase']
	
	bcs = define_boundary_conditions(reac_input['reactor_type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['reactants_num'],all_molecules,reac_input['reactant_ratio'],boundary_L,boundary_R,reac_input['number_of_inerts'])
	
	### Initialize the graphs and set the graphical parameters ###
	fig2,ax2,legend_label,header,colors = establish_output_graph(reac_input['reactor_type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['number_of_inerts']))

	if reac_input['fit_parameters'].lower() == 'true':
		if reac_input['objective_points'] == 1:
			output_fitting = exp_data_fitting_1(legend_label,reac_input['Time_steps'],reac_input['exp_data_folder'],reac_input['Time'])
		
		elif reac_input['objective_points'] == 3:
			output_fitting = exp_data_fitting_3(legend_label,reac_input['Time_steps'],reac_input['exp_data_folder'])

		elif reac_input['objective_points'] == 4:
			output_fitting = exp_data_fitting_4(legend_label,reac_input['Time_steps'],reac_input['exp_data_folder'])

		elif reac_input['objective_points'] == 5:
			output_fitting = exp_data_fitting_5(legend_label,reac_input['exp_data_folder'],reac_input['Time_steps'],reac_input['Time'])
		
		elif type(reac_input['objective_points']) == str:
			points = reac_input['objective_points'].split(',')
			output_fitting = exp_data_fitting_user(legend_label,reac_input['exp_data_folder'],reac_input['Time_steps'],reac_input['Time'])

		else:
			print('Objective Point Input Is Not Valid')
			sys.exit()

	if reac_input['sensitivity_analysis'].lower() == 'true':
		path_2 = reac_input['output_folder_name']+'_folder/sensitivity/'
		generate_folder(path_2)
		sens_time_list = []
		
		for k_gas in legend_label[:-1]:
			path_molecules = path_2+k_gas+'/'
			generate_folder(path_molecules)
		
		for j in r_const:
			r_const[j] = Constant(r_const[j])
	
	to_flux = flux_generation(reac_input['reactor_type'],monitored_gas,reac_input['reactants_num'],reac_input['Inert_pulse_size'],D,eb,dx_r,reac_input['reac_radius'],dx2_r)
	store_data_func(reac_input['store_data'],reac_input['output_folder_name'])
	store_sens_analysis(reac_input['sensitivity_analysis'],reac_input['output_folder_name'],monitored_gas,legend_label)
	
	### Declare the solver to be used ###
	
	#class Plane(SubDomain):
	#	def inside(self, x, on_boundary):
	#		return x[0] > 1.0 - DOLFIN_EPS

	J = derivative(F,u)
	problem = NonlinearVariationalProblem(F,u,bcs,J)
	solver = NonlinearVariationalSolver(problem)
	
	solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
	solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10

	### Handle a pulse intensity that changes with each pulse ###
	if '/' in str(reac_input['reactant_ratio']):
		list_of_feed = reac_input['reactant_ratio'].split('/')
		list_of_time = reac_input['pulse_time'].split('/')
		pulse_variation = len(list_of_feed)
		pulse_time = len(list_of_time)
		reactant_feed = []
		reactant_time = []
		for k in range(0,len(reac_input['reactant_ratio'].split('/'))):
			reactant_feed.append(list_of_feed[k].split(','))
			reactant_time.append(list_of_time[k].split(','))
	else:
		pulse_variation = 1
		reactant_feed = []
		reactant_time = []
		try:
			reactant_feed.append(reac_input['reactant_ratio'].split(','))
			reactant_time = list(map(float, reac_input['pulse_time'].split(',')))
		except AttributeError:
			reactant_feed.append(reac_input['reactant_ratio'])
			reactant_time = list(map(float, reac_input['pulse_time'].split(',')))

	### Begin simulation for each desired pulse ###
	current_reac_set = 1
	tot_objective = 0

	for k_pulse in range(0,int(reac_input['number_of_pulses'])):

		### Clear the tape for the adjoint (will lead to challenges with fitting parameters over multiple pulses) ###
		tape.clear_tape()
		
		print("")
		print("Simulation Status: "+"Pulse #"+str(k_pulse+1))

		species_pulse_list = reactant_feed[current_reac_set-1]
		species_time = reactant_time[current_reac_set-1]
		### Incase the time step was altered during a previous pulse, just reset to the original values ### 
		dk = Constant(reac_input['Time']/reac_input['Time_steps'])
		dt = reac_input['Time']/reac_input['Time_steps']
		
		start_time = time.time()
		
		### Redefine the lists tracking all the information ###
		sensitivity_output = {}
		for k_sens in range(monitored_gas):
			sensitivity_output[k_sens] = []
		
		graph_data = {}
		for k_gasses in range(0,necessary_values['molecules_in_gas_phase']):
			graph_data['conVtime_'+str(k_gasses)] = []
		for kjc in range(0,int(reac_input['number_of_inerts'])):
			graph_data['conVtime_'+str((necessary_values['gas_num']+1)+necessary_values['surf_num']+kjc)] = []
		
		surf_data['timing'] = []
		
		for j_gasses in range(necessary_values['molecules_in_gas_phase'],len(necessary_values['reactants'])-1):
			surf_data['conVtime_'+str(j_gasses)] = []
		
		graph_data['timing'] = []

		cat_data['timing'] = []
		
		for z_gasses in range(0,all_molecules):
			cat_data['conVtime_'+str(z_gasses)] = []

		t = 0
		if reac_input['fit_parameters'].lower() == 'true' or reac_input['sensitivity_analysis'].lower() == 'true' or reac_input['ad_rrm'].lower() == 'true':
			osub = integration_section()
			domains = CellFunction("size_t", mesh)
			domains.set_all(0)
			osub.mark(domains, 1)
			dP = Measure('vertex',domain = mesh, subdomain_data=domains)
		
		test_new = project(-u.dx(0),V)
		test_new_split = test_new.split(deepcopy=True)

		### Used to define the function to be minimized between the experimental points and synthetic points ###
		w_new = Expression("1", degree=0)
		w_new2 = interpolate(w_new,V_du)
		W = VectorFunctionSpace(mesh, 'P', 1)
		w = Function(W)
		w.vector()[:] = 10

		x_dim = list(range(0, int(reac_input['mesh_size'])+1))
		
		cum_sum = 0
		
		while t <= reac_input['Time']:

			### Store the results of each iteration ###
			graph_data['timing'].append(t)
			for k in range(0,monitored_gas):
				new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k]))
				graph_data['conVtime_'+str(k)].append((new_val))
			#for kjc in range(all_molecules-int(reac_input['number_of_inerts']),all_molecules):
			for kjc in range(0,int(reac_input['number_of_inerts'])):
				new_val = ((to_flux[monitored_gas]*u.vector().get_local()[2*(all_molecules+1)-2-(int(reac_input['number_of_inerts'])-kjc)]))
				graph_data['conVtime_'+str(all_molecules-(int(reac_input['number_of_inerts'])-kjc))].append((new_val))
			

			mesh_size = reac_input['mesh_size']

			### If you are trying to fit parameters, then check if the point must be defined as a part of the objective function ###
			if reac_input['fit_parameters'].lower() == 'true':

				for k_fitting in range(0,len(legend_label[:-1])):
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
								jfunc_2 += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
								tot_objective += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
							else:
								pass

						except UnboundLocalError:
							if legend_label[k_fitting] != 'Inert':
								jfunc_2 = assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
								tot_objective += assemble(inner(u_n[k_fitting] - w3,u_n[k_fitting] - w3)*dP(1))
							else:
								pass

			if reac_input['mkm_analysis'].lower() == 'true' or reac_input['RRM_analysis'].lower() == 'true' or reac_input['petal_plots'].lower() == 'true':
				for jk in range(0,all_molecules):
					new_values = [] 
					cat_values = []

					#for cj in range(0,all_molecules):
					#	cat_values.append(0)

					#for kj in range(monitored_gas+1,len(necessary_values['reactants'])+1):
					#	new_values.append(0)
			
					for z in range(int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1,int((cat_location + 0.5*frac_length)*reac_input['mesh_size'])):
						
						try:
							###for z_num in range(0,1+necessary_values['molecules_in_gas_phase']):#.split(',')
							#print(z*(all_molecules+1)-(z_num))
							value_cat = u_n.vector().get_local()[z*(all_molecules)-(all_molecules)+jk]# = float(z_sur)
							cat_values.append(value_cat)
							###cat_values[z_num] += value_cat*dx_r
							#surf_data['conVtime_'+str(z_num+monitored_gas)] += value_new*dx_r
						except AttributeError:
							value_cat = u_n.vector().get_local()[z*(all_molecules)]# = float(z_sur)
							cat_values[0] += value_cat*dx_r
						#value_cat = u_n.vector().get_local()[z*(all_molecules+1)-all_molecules]
						#cat_values[necessary_values['molecules_in_gas_phase']] += value_cat*dx_r

						###try:
						###	swap_for_2 = necessary_values['molecules_in_gas_phase']
						###	for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
						###		#print(z*(all_molecules+1)-(swap_for_2+z_num))
						###		value_new = u_n.vector().get_local()[z*(all_molecules)-(swap_for_2+z_num)]# = float(z_sur)
						###		new_values[z_num] += value_new*dx_r
						###		#surf_data['conVtime_'+str(z_num+monitored_gas)] += value_new*dx_r
						###except AttributeError:
						###	value_new = u_n.vector().get_local()[z*(all_molecules)]# = float(z_sur)
						###	new_values[0] += value_new*dx_r
					cat_data['conVtime_'+str(jk)].append(cat_values)

			#cat_length = reac_input['length_reac']*reac_input['factor']
			#if ',' in str(reac_input['rangesurface_species']):
			#	print(reac_input['rangesurface_species'].split(','))
			#	for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
			#		cat_data['conVtime_'+str(z_num)].append(float(z_sur)/(cat_length))
			#else:
			#	cat_data['conVtime_'+str(0)].append(cat_values[0])

			#if ',' in str(reac_input['rangesurface_species']):
			#	for z_num,z_sur in enumerate(reac_input['ranges2urface_species'].split(',')):
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
			
	
			##if t > 0:
			if round(t,6) not in reactant_time:
				dt = solver_iteration(dt,reac_input['solver_method_options'],solver,dk,1.5,1.1)
				
			else:
				if reac_input['reactor_type'] == 'tap':
					if ',' in str(species_pulse_list):
						for k in range(0,int(reac_input['reactants_num'])):
							if reactant_time[k] == round(t,6):
							###u_n.vector()[int((all_molecules)*(reac_input['mesh_size']+1)-1)+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]
								u_n.vector()[int((all_molecules)*(reac_input['mesh_size']+1))+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]

					else:
						u_n.vector()[int((all_molecules)*(reac_input['mesh_size']+1)-1)-(all_molecules)] = float(species_pulse_list)*Inert_pulse_conc
					for k_step in range(0,int(reac_input['number_of_inerts'])):
						##u_n.vector()[int((all_molecules)*(reac_input['mesh_size']+1)-1-k_step)] = float(species_pulse_list[-1-k_step])*Inert_pulse_conc###??? Added the porosity contribution
						if reactant_time[-1-k_step] == round(t,6):
							u_n.vector()[int((all_molecules)*(reac_input['mesh_size']+1)-1-k_step)] = float(species_pulse_list[-1-k_step])*Inert_pulse_conc###??? Added the porosity contribution
					
				if k_pulse == 0:

					for z in range(int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1,int((cat_location + 0.5*frac_length)*reac_input['mesh_size'])):
						if ',' in str(reac_input['rangesurface_species']):
							###for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
							for z_num,z_sur in enumerate(reac_input['rangesurface_species']):	
								if int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1 <= 0:
									#u_n.vector()[(all_molecules)-(2+z_num)] = float(z_sur)
									u_n.vector()[(all_molecules)-(int(reac_input['number_of_inerts'])+z_num+1)] = float(z_sur)
								else:
									###u_n.vector()[z*(all_molecules)-(2+z_num)] = float(z_sur)
									u_n.vector()[z*(all_molecules)-(int(reac_input['number_of_inerts'])+z_num+1)] = float(z_sur)
						else:
							if int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1 <= 0:
								u_n.vector()[(all_molecules)-(2)] = float(species_pulse_list)
							else:
								u_n.vector()[z*(all_molecules)-(2)] = float(species_pulse_list)

				dt = solver_iteration(dt,reac_input['solver_method_options'],solver,dk,1.5,1.1)

			if reac_input['sensitivity_analysis'].lower() == 'true':
				time_size = time.time()
				u_final = u.split(deepcopy=False)
				should_it = int(round(t*reac_input['Time_steps']/reac_input['Time'],0))
				if should_it <= 50 and should_it%5 == 0:
					for k_step in range(0,monitored_gas):
						temp = call_sens_analysis(u_final[k_step],controls,dP(1))
						sensitivity_output[k_step].append(temp)
					simulation_time_list.append(time.time() - time_size)
					sens_time_list.append(t)


			if reac_input['ad_rrm'].lower() == 'true':
				time_size = time.time()
				u_final = u.split(deepcopy=False)
				should_it = int(round(t*reac_input['Time_steps']/reac_input['Time'],0))
				if should_it <= 50 and should_it%5 == 0:
					for k_step in range(0,monitored_gas):
						temp = call_ad_rrm_analysis(u_final[k_step],controls,dP(1))
						sensitivity_output[k_step].append(temp)
					simulation_time_list.append(time.time() - time_size)
					sys.exit()
					sens_time_list.append(t)

			progressBar(t, reac_input['Time'])
			u_n.assign(u)
			t += dt

		current_reac_set += 1
		if current_reac_set >  pulse_variation:
			current_reac_set = 1

		x_values = []
		it_times = []

		def print_fun(x):
			it_times.append(time.time())
			x_values.append(x.tolist())
			print(time.time())
			print(x)

		set_log_active(False)
		fitting_time = time.time()
		
		if reac_input['fit_parameters'].lower() == 'true':
			rf_2 = ReducedFunctional(jfunc_2, controls)
			low_bounds = []
			up_bounds = []
			try:
				for gt in range(0,len(controls)):
					low_bounds.append(0)
					up_bounds.append(np.inf)
				if reac_input['optimization_method'] == 1 or reac_input['optimization_method'] == '':
					u_opt_2 = minimize(rf_2, callback=print_fun, bounds = (low_bounds,up_bounds),tol=1e-13, options={"ftol":1e-13,"gtol":1e-13})
				elif reac_input['optimization_method'] == 2:
					u_opt_2 = minimize(rf_2, method = 'Newton-CG', callback=print_fun,tol=1e-13, options={"ftol":1e-13,"gtol":1e-13})
				elif reac_input['optimization_method'] == 3:
					u_opt_2 = minimize(rf_2, method = 'BFGS', callback=print_fun,tol=1e-13, options={"gtol":1e-13})
				elif reac_input['optimization_method'] == 4:
					u_opt_2 = minimize(rf_2, method = 'SLSQP', callback=print_fun, bounds = (low_bounds,up_bounds),tol=1e-13, options={"ftol":1e-13})
				elif reac_input['optimization_method'] == 5:
					u_opt_2 = minimize(rf_2, method = 'CG', callback=print_fun ,tol=1e-13, options={"gtol":1e-13})
				#elif reac_input['optimization_method'] == 6:
				#	u_opt_2 = minimize(rf_2, method = 'basinhopping', callback=print_fun, bounds = (low_bounds,up_bounds),tol=1e-13, options={"ftol":1e-13,"gtol":1e-13})
			
				else:
					print('Requested Optimization Method Does Not Exist')
					sys.exit()
			except RuntimeError:
				print('Optimization Failed or Entered Too Stiff Region')
				time.sleep(1.5)
				print('Will shortly begin generating the optimization gif.')
				time.sleep(5)
			#Extra:
			#u_opt = minimize(rf, method="COBYLA", bounds = ((0,0,0),(np.inf,np.inf,np.inf)), callback=print_fun,tol = 1e-10)#, options={"disp": True} # , method="COBYLA", bounds = ((0,0),(np.inf,np.inf), callback=print_fun,tol = 1e-5
			#u_opt = minimize(rf, method="TNC", callback=print_fun,tol = 1e-5, options={"disp": True})#,options={'gtol': 1e-05})#  , bounds=(0,1000)  #method="L-BFGS-B" ###, , constraints={'type':'ineq', 'fun': rf}
				

			with open('./'+reac_input['output_folder_name']+'_folder/fitting/your_file.txt', 'w') as f:
				f.write("Contents: "+str(it_times))
				f.write('\n')
				f.write("Constants: "+str(x_values))

				f.close

	############# Store the output data #######################################################
		if k_pulse == 0:
			dictionary_of_numpy_data = {}
			for j_species in range(0,monitored_gas+int(reac_input['number_of_inerts'])):
				dictionary_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionary_of_numpy_data[legend_label[j_species]] = np.vstack((dictionary_of_numpy_data[legend_label[j_species]],new_data))
			
			for k_species in range(monitored_gas+1,len(necessary_values['reactants'])+1):
				dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.vstack((dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]],new_data))
			
			#dictionary_of_gas_cat_data = {}
			
			#for z_species in range(0,all_molecules):
			#	dictionary_of_gas_cat_data[necessary_values['reactants'][z_species]] = np.empty(shape=[0, len(graph_data['timing'])])
			#	new_data = np.asarray(graph_data['timing'])
			#	dictionary_of_gas_cat_data[necessary_values['reactants'][z_species]] = np.vstack((dictionary_of_gas_cat_data[necessary_values['reactants'][z_species]],new_data))

#for z in range(int((cat_location - 0.5*frac_length)*reac_input['mesh_size'])-1,int((cat_location + 0.5*frac_length)*reac_input['mesh_size'])):
#	for z_num,z_sur in enumerate(reac_input['rangesurface_species'].split(',')):
#		value_new = u_n.vector().get_local()[z*(all_molecules+1)-(2+z_num)]# = float(z_sur)
#		surf_data['conVtime_'+str(z_num+monitored_gas)] += value_new*dx_r

	############# Visualize/graph the data #######################################################
		if reac_input['show_experimental_data'].lower() == 'true':
			data = pd.ExcelFile('../example_data/CO_oxidation_112.5_C.base_corrected.xlsx')
			df3 = pd.read_excel(data,'Averaged_last_50_pulse',header=None)
			df4 = pd.read_excel(data,'4AMU',header=None)
			df5 = pd.read_excel(data,'40AMU',header=None)
			df6 = pd.read_excel(data,'28AMU',header=None)
			df7 = pd.read_excel(data,'32AMU',header=None)
			df8 = pd.read_excel(data,'44AMU',header=None)

			ax2.plot(df3[0][:].convert_objects(convert_numeric=True),df3[2][:].convert_objects(convert_numeric=True),label='av CO',color = 'b')
			ax2.plot(df3[0][:].convert_objects(convert_numeric=True),df3[3][:].convert_objects(convert_numeric=True),label='av O2',color = 'y')
			ax2.plot(df3[0][:].convert_objects(convert_numeric=True),df3[1][:].convert_objects(convert_numeric=True),label='av CO2',color = 'g')
			ax2.plot(df3[0][:].convert_objects(convert_numeric=True),df3[4][:].convert_objects(convert_numeric=True),label='av Ar',color = 'r')
			ax2.plot(df3[0][:].convert_objects(convert_numeric=True),df3[5][:].convert_objects(convert_numeric=True),label='av He',color = 'k')

			#ax2.plot(df3[0][:1500].convert_objects(convert_numeric=True),df6[0][:1500].convert_objects(convert_numeric=True),label='CO 1',color='b')
			#ax2.plot(df3[0][:1500].convert_objects(convert_numeric=True),df7[62][:1500].convert_objects(convert_numeric=True),label='O2 1',color='y')
			#ax2.plot(df3[0][:1500].convert_objects(convert_numeric=True),df8[0][:1500].convert_objects(convert_numeric=True),label='CO2 1',color='g')
			#ax2.plot(df3[0][:1500].convert_objects(convert_numeric=True),df5[0][:1500].convert_objects(convert_numeric=True),label='Ar 1',color='r')		
			#ax2.plot(df3[0][:1000].convert_objects(convert_numeric=True),df6[2][:1000].convert_objects(convert_numeric=True),label='He 1',color='k')

		sig = 0.05
		mu = 0.0434
		
		alpha_1 = 0.0281
		alpha_2 = 0.0027
		
		beta_1 = 0.00166
		beta_2 = 0.00270
		
		theta_1 = 0
		theta_2 = 0
		theta_3 = 0
		
		w_1 = 2*3.14159*60
		w_2 = 2*3.14159*70
		w_3 = 2*3.14159*120

		for k,j in enumerate(graph_data):
			if j != 'timing':
				if reac_input['noise'].lower() == 'true':
					
					for z in range(0,int(reac_input['Time_steps'])):
						graph_data[j][z] += np.random.normal(0,1)*sig +beta_2*np.cos(w_2*(k*dt))
				if k_pulse > 0:
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
				else:
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)#/reac_input['Inert_pulse_size']*(reac_input['reactants_num']+1)
				new_data = np.asarray(graph_data[j])
				dictionary_of_numpy_data[legend_label[k]] = np.vstack((dictionary_of_numpy_data[legend_label[k]],new_data))
			else:
				pass
		
		if reac_input['sensitivity_analysis'].lower() == 'true':
			for k_sens_step in range(monitored_gas):
				sens_time = np.asarray(sens_time_list)
				#sens_time = np.asarray(graph_data['timing'][0:])
				sens_time = sens_time.T#np.transpose(sens_time)
				sensitivity_output_2 = np.asarray(sensitivity_output[k_sens_step])
				sensitivity_output_2 = np.append(sens_time[:,None],sensitivity_output_2,axis=1)	
				np.savetxt(reac_input['output_folder_name']+'_folder/sensitivity/'+legend_label[k_sens_step]+'/pulse_'+str(k_pulse+1)+'.csv',sensitivity_output_2,delimiter=",",header='t,'+','.join(legend_2))#
			df_sens_time = pd.DataFrame(simulation_time_list)
			df_sens_time.to_csv(reac_input['output_folder_name']+'_folder/sensitivity/time.csv',header=None)			

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
		#		dictionary_of_numpy_data[name_list[k-1]] = np.vstack((dictionary_of_numpy_data[name_list[k-1]][0],new_data))
		#	else:
		#		pass

	ax2.legend(title="Gas Species")
	plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
	#ax2.set_ylim(0,0)
	if reac_input['save_figure'].lower() == 'true':
		plt.savefig('./'+reac_input['output_folder_name']+'_folder/graphs/flux_data.png')
	
	if reac_input['Display_figure'].lower() == 'true':
		plt.show()
	
	for j_species in range(0,monitored_gas+int(reac_input['number_of_inerts'])):
		dictionary_of_numpy_data[legend_label[j_species]] = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
		np.savetxt('./'+reac_input['output_folder_name']+'_folder/flux_data/'+legend_label[j_species]+'.csv', dictionary_of_numpy_data[legend_label[j_species]], delimiter=",")
	###Current!!!
	if reac_input['mkm_analysis'].lower() == 'true' or reac_input['RRM_analysis'].lower() == 'true' or reac_input['petal_plots'].lower() == 'true':
		for j_species in range(0,all_molecules-int(reac_input['number_of_inerts'])):
			print(necessary_values['reactants'][j_species])
			np.savetxt('./'+reac_input['output_folder_name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.array(cat_data['conVtime_'+str(j_species)]), delimiter=",")


		#dictionary_of_gas_cat_data[legend_label[j_species]] = np.transpose(dictionary_of_gas_cat_data[legend_label[j_species]])
		#np.savetxt('./'+reac_input['output_file_name']+'_folder/thin_data/'+legend_label[j_species]+'.csv', dictionary_of_gas_cat_data[legend_label[j_species]], delimiter=",")

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

	if reac_input['fit_parameters'].lower() == 'true':
	###if reac_input['output_folder_name'] == '':
	#if reac_input['fit_parameters'].lower() == 'true':
		with open('./'+reac_input['output_folder_name']+'_folder/fitting/your_file-A.txt', 'r') as f:
			lines = f.readlines()
		f.close
		lines = [x.strip() for x in lines]
		times = lines[0]
		times = times.replace('Contents: ','')
		times = eval(times)

		constants = lines[1]
		constants = constants.replace('Constants: ','')
		constants = eval(constants)
		things = len(times)

		for k_num in range(0,things):
			#print('what?')
			alter = pd.read_csv('./input_file.csv',header=None)

			variables_to_change = ['Display_figure','fit_parameters','sensitivity_analysis','store_data','output_folder_name','Reaction_Information']
			
			for k in range(alter.shape[0]):
				if alter[0][k] in variables_to_change:
					if alter[0][k] == 'store_data':
						alter.iloc[k,1] = 'TRUE'
					elif alter[0][k] == 'output_folder_name':
						alter.iloc[k,1] = reac_input['output_folder_name']+'_folder/fitting/iter_'+str(k_num)
					elif alter[0][k] == 'Reaction_Information':
						value = 1
						current_rate = k+1
						kz = 0
						while kz < len(float_constants):
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
			#alter.to_csv('./input_file.csv',header=None,index=False)
			alter.to_csv('./input_file.csv',header=None,index=False)
			
			#sys.exit()
			
			
			call_sim()

		#generate_gif(legend_label[:len(legend_label)], reac_input['exp_data_folder']+'/flux_data', './'+reac_input['output_folder_name']+'_folder/fitting/', len(constants), constants, reactor_kinetics_input['reactions_test'], times)
		
		for k_num in range(0,things):
			shutil.rmtree('./'+reac_input['output_folder_name']+'_folder/fitting/iter_'+str(k_num)+'_folder') 


	return graph_data, legend_label, necessary_values['reactants']

def call_sim():
	
	reactor_kinetics_input,kinetic_parameters,kin_in = read_input()

	graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)

	###if reactor_kinetics_input['mkm_analysis'].lower() == 'true':
	###	MKM_graphs(kin_in,reactor_kinetics_input['reactions_test'],reactor_kinetics_input['output_folder_name']+'_folder/thin_data',reactor_kinetics_input['Display_figure'].lower())
	###if reactor_kinetics_input['petal_plots'].lower() == 'true':
	###		print('true_mkm')
	###if reactor_kinetics_input['RRM_analysis'].lower() == 'true':
	###	R_RRM_func(legend_label[0:len(legend_label)-1],os.getcwd(),reactor_kinetics_input['output_folder_name']+'_folder')
	###	if reactor_kinetics_input['petal_plots'].lower() == 'true':
	###		print('true_RRM')
	###if reactor_kinetics_input['mkm_analysis'].lower() == 'true' and reactor_kinetics_input['RRM_analysis'].lower() == 'true':
	###	jacobian_visual(kin_in,reactor_kinetics_input['reactions_test'],reactor_kinetics_input['output_folder_name']+'_folder/thin_data',reactor_kinetics_input['Display_figure'].lower(),legend_label[0:len(legend_label)-1],os.getcwd(),reactor_kinetics_input['output_file_name']+'_folder',legend_label[:-1],reactor_kinetics_input['number_of_pulses'])

call_sim()
