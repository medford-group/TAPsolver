from fenics import *
from fenics_adjoint import *
from funs_tap_sim import *
from variational_form_constructor import make_f_equation
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

	######## Reactor type & data storage ########################################################################


def tap_simulation_function(reactor_kinetics_input):
	reac_input = reactor_kinetics_input
	reactor_type = reac_input['reactor_type'] 	#currently t_pfr / tap / t_pfr_diff / batch {not ready just yet}
	output_file_name = reac_input['output_file_name'] # "./FILENAME.csv"
	theta = reac_input['theta'] 				# forward_euler = 0, backward_euler = 1, crank_nicolson = 1/2
	solver_method_options = reac_input['solver_method_options']	#'simple_adaptive' or None 
	Time = reac_input['Time']				# Total pulse time
	Time_steps = reac_input['Time_steps']		# Number of time steps  
	mesh_size = reac_input['mesh_size']	# Want to test and store all of these
	T = reac_input['T']
	sensitivity_analysis = reac_input['sensitivity_analysis'] # True = on / False = off
	Frequency_of_sensitivity = reac_input['Frequency_of_sensitivity']
	Display_figure = reac_input['Display_figure']
	save_figure = reac_input['save_figure']
	store_data = reac_input['store_data']
	
	############# Feed composition ############################################################################
	
	reactants_num = reac_input['reactants_num']
	Inert_pulse_size = reac_input['Inert_pulse_size']	###Size of inert pulse (# of molecules, not mol)#5e-8 #(5e-9)*
	reactant_ratio = reac_input['reactant_ratio']		###Ratio of reactant to the size of Inert_pulse_size
	number_of_pulses = reac_input['number_of_pulses']
	
	############# Reaction Equations ############################################################################
	
	### Pt oxidation
	# actual
	#reactions_test = ['O2 + 2* -> 2O*']
	# Toy model
	#reactions_test = ['A + * -> A*']
	
	### CO OXIDATION
	#Eley
	#reactions_test = ['A + * -> A*','B + A* -> * + C']
	#Eley with preadsorbed
	#reactions_test = ['CO + O* ->  CO2 + *']
	#Lang
	#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*']
	#Eley & Lang
	#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*','CO + O* -> CO2 + *']
	#El_La_aop
	#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*','CO + O* -> CO2 + *','O* + z* -> zO* + *']
	#Reece reaction network
	reactions_test = reac_input['reactions_test']
	
	
	### Ammonia decomposition 
	#reactions_test = ['NH3 + * <-> NH3*', 'H2 + 2* <-> 2H*', 'N2 + * <-> N2*', 'N2* + * <-> 2N*', 'NH3* + * <-> NH2* + H*', 'NH2* + * <-> NH* + H*', 'NH* + * <-> N* + H*']
	
	############# Reactor Parameters ############################################################################
	
	len_inert_1 = reac_input['len_inert_1']#2.54*(19/20)/2		###Length of the first inert zone (cm)
	len_cat = reac_input['len_cat']#2.54*(1/20)			###Length of the catalyst zone (cm) 
	len_inert_2 = reac_input['len_inert_2']#2.54*(19/20)/2		###Length of the second inert zone (cm)
	reac_radius = reac_input['reac_radius']#sqrt((1.8e-5)*(100**2)/3.14159)		###Radius of the reactor
	
	############# Transport properties ##########################################################################
	
	void_inert = reac_input['void_inert']		###Voidage of the catalyst
	void_cat = reac_input['void_cat']			###Voidage of the catalyst
	
	ref_rate_inert = reac_input['ref_rate_inert']		###Reference 
	ref_rate_cat = reac_input['ref_rate_cat']
	ref_mass = reac_input['ref_mass']
	#mass_list = np.array((16,16))
	mass_list = reac_input['mass_list']	
	#mass_list = np.array((17.03,2.015,28.01,40))
	velocity = reac_input['velocity']
	
	############# Catalyst properties ##########################################################################
	
	#Fraction of sites
	mcat = reac_input['mcat']   			###Catalyst mass in grams
	ssarea = reac_input['ssarea']			###Surface area of support in cm2/g
	bsarea = reac_input['bsarea']	###Total surface area of support
	sdensity = reac_input['sdensity']	###Surface Density (atoms/cm2)
	atloading = reac_input['atloading']		### Fractional metal loading of catalyst (assume same surface density as support)
	sarea = reac_input['sarea']	### Total surface area of metal
	open_sites = reac_input['open_sites']
	OA = reac_input['OA'] 		### OA site density (atoms/cm2)
	OB = reac_input['OB']	
	
	#rangesurface_species = [(0.25/3000)*6.022e23,(0.5/1000)*6.022e23]
	#rangesurface_species = [(5e-9)*6.022e23,(5e-9)*6.022e23]
	#rangesurface_species = [((1e7)*(5e-9)/(100**2))*6.022e23]
	rangesurface_species = reac_input['rangesurface_species']
	
	########## Input Rate Constants #############################################################################
	
	###### Time stepping/solver options
	
	#Ke0 = Constant(reac_input['Ke0'])
	
	#Ke0 = Constant(reac_input['Ke0'])
	#Kd0 = Constant(reac_input['Ke0'])
	
	#Ke1 = Constant(reac_input['Ke0'])		###Rate constants
	
	#Ke2 = Constant(reac_input['Ke0'])
	#Ke3 = Constant(reac_input['Ke0'])
	
	Ke0 = Constant(reac_input['Ke0'])
	Ke1 = Constant(reac_input['Ke1'])
	Kd0 = Constant(reac_input['Kd0'])
	Ke1 = Constant(reac_input['Ke1'])
	Kd1 = Constant(reac_input['Kd1'])
	Ke2 = Constant(reac_input['Ke2'])
	Kd2 = Constant(reac_input['Kd2'])
	Ke3 = Constant(reac_input['Ke3'])
	Kd3 = Constant(reac_input['Kd3'])
	Ke4 = Constant(reac_input['Ke4'])
	Kd4 = Constant(reac_input['Kd4'])
	Ke5 = Constant(reac_input['Ke5'])
	Kd5 = Constant(reac_input['Kd5'])
	Ke6 = Constant(reac_input['Ke6'])
	Kd6 = Constant(reac_input['Kd6'])
	
	#controls = [Control(Ke0),Control(Kd0),Control(Ke1),Control(Ke2),Control(Ke3)]
	#legend_2 = ['Ke0','Kd0','Ke1','Ke2','Ke3']
	
	
	############# Used to control the output/warnings of Fenics ############################################################################
	
	parameters["std_out_all_processes"] = False																							
	cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
	set_log_active(False)
	import warnings
	warnings.filterwarnings("ignore", category=DeprecationWarning)
	tol = 1E-14
	
	rangesurface_species = list(reversed(rangesurface_species))
	
	############# Establish transport properties in FEniCS ##################################################################################
	
	original_pulse_size = Inert_pulse_size
	dk = Constant(Time/Time_steps)
	eb = np.array((void_inert,void_cat,void_inert))
	
	ref_rate = np.append(ref_rate_inert,ref_rate_cat)
	ref_rate = np.append(ref_rate,ref_rate[0])  	###Add inert rate to end of ref_rate
	
	D = np.empty((len(mass_list),3))            ###Generate numpy matrix for diffusion coefficients
	def diff_func(ref_m,mol_mass,ref_r):     	###Define diffusion equation (having as a function might not be necessary)
		return ref_r*(mp.sqrt(40*T)/mp.sqrt(423.15*mol_mass))
	for k,j in enumerate(mass_list):			### Use for loop to generate diffusion coefficient
		for k_2,j_2 in enumerate(ref_rate):
			D[k,k_2] = diff_func(ref_mass,j,j_2)
	print(D)
	def vect_vol(x,y):
		return x*(reac_radius**(2))*np.pi*y
	
	r_param, dx_r, dx2_r, frac_length, cat_location = establish_grid_system(len_inert_1,len_cat,len_inert_2,mesh_size)
	
	############# Establish reactor parameters in FEniCS ###################################################################################
	
	ca = (reac_radius**2)*3.14159				###Cross sectional area of reactor
	point_volume = dx_r * ca 					###Volume at a specific point in the reactor
	specific_area = 1e8
	surface_area_at_point = point_volume*specific_area
	open_sites_per_point = sdensity*surface_area_at_point
	
	vect_vol_imp = np.vectorize(vect_vol)
	vr = np.sum(vect_vol_imp(r_param,eb))		### Volume of the reactor
	vc = r_param[1]*(reac_radius**(2))*np.pi*eb[1]		### Volume of the catalyst zone
	Q = sarea/vc 			### surface to gas scaling coefficient (1/cm)
	Inert_pulse_size = Inert_pulse_size/(dx_r*ca)
	
	necessary_values = make_f_equation(reactions_test,reactants_num,reactor_type,False)
	
	############# Establish the mesh and boundaries/domains in FEniCS #######################################################
	
	
	
	mesh = UnitIntervalMesh(int(mesh_size))
	P1 = FiniteElement('P',mesh.ufl_cell(),1)
	
	if reactions_test != ['INERT_ONLY']:
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
	
	#def pulse_point(x, on_boundary):
	#	return near(x[0],1/grid_points,tol)
	
	class integration_section(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (0.99,1.0))
	
	osub = integration_section()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	osub.mark(domains, 1)
	dP = Measure('vertex',domain = mesh, subdomain_data=domains)
	
	
	all_molecules = necessary_values['gas_num']+necessary_values['surf_num']+1
	u = Function(V)
	u_n = Function(V)
	u_temp = Function(V)
	u_1 = Function(V)
	u_2 = Function(V)
	u_3 = Function(V)
	
	graph_data,v_d,u_d,u_nd,sens_data = initialize_variable_dictionaries(necessary_values,all_molecules,V,u,u_n)
	
	dt = Time/Time_steps
	Length = np.sum(r_param)**2
	
	W = VectorFunctionSpace(mesh, 'P', 1)
	w = Function(W)
	w.vector()[:] = velocity
	
	
	
	############# Declare the equations used to describe the system in FEniCS #######################################################
	
	try:
		F = eval(necessary_values['F'])
	except NameError:
		error_output(reactions_test)
	#### Just an equations for reeces code (different from others since it requires the use of a scaling parameter ( Q ))
	#Q = Constant(5.5423e3)
	#F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)  + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)  + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0)  + dk*Constant(D[1][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1)  + dk*Constant(D[1][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + ((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + ((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + ((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + ((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + ((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + ((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1) + ((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + ((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1)- dk*(1.0* Q*Kd0*(u_d['u_3']**1.0)*v_d['v_1']*dx(1))+ dk*(1.0* Ke0*(u_d['u_1']**1.0)*v_d['v_1']*dx(1))+ dk*(1.0* Kd0*(u_d['u_3']**1.0)*v_d['v_3']*dx(1))- dk*(1.0* (Ke0/Q)*(u_d['u_1']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_4']*dx(1))- dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_6']*dx(1))- dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke3*(u_d['u_5']**1.0)*v_d['v_5']*dx(1))- dk*(1.0* Q*Ke3*(u_d['u_5']**1.0)*v_d['v_2']*dx(1))+ Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)+ dk * Constant(D[2][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0) + Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1) + dk * Constant(D[2][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)
	
	monitored_gas = necessary_values['molecules_in_gas_phase']
	
	bcs = define_boundary_conditions(reactor_type,reactions_test,necessary_values['molecules_in_gas_phase'],V,reactants_num,all_molecules,reactant_ratio,boundary_L,boundary_R)
	############# Initialize the graphs and set the graphical parameters #######################################################
	
	fig2,ax2,legend_label,header,colors = establish_output_graph(reactor_type,necessary_values['molecules_in_gas_phase'],necessary_values['reactants'])
	
	to_flux = flux_generation(reactor_type,monitored_gas,reactants_num,original_pulse_size,D,eb,dx_r,reac_radius,dx2_r)
	store_data_func(store_data,output_file_name)
	store_sens_analysis(sensitivity_analysis,output_file_name,monitored_gas,legend_label)
	############# Declare the solver to be used #######################################################
	
	J = derivative(F,u)
	problem = NonlinearVariationalProblem(F,u,bcs,J)
	solver = NonlinearVariationalSolver(problem)
	
	solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
	solver.parameters["newton_solver"]["absolute_tolerance"] = 1e3
	############# Define a method to iteratively choose time-step (will change based on meeting with tobin isaac) ##
	#%#%#%#%#% TASK #3	
	############# Solve the system #######################################################
	for k_pulse in range(0,number_of_pulses):
		tape.clear_tape()
		
		print("")
		print("Simulation Status: "+"Pulse #"+str(k_pulse+1))
		
		dk = Constant(Time/Time_steps)
		dt = Time/Time_steps
		
		start_time = time.time()
		
		sensitivity_output = {}
		for k_sens in range(monitored_gas):
			sensitivity_output[k_sens] = []
		
		graph_data['timing'] = []
		for k_gasses in range(0,necessary_values['molecules_in_gas_phase']):
			graph_data['conVtime_'+str(k_gasses)] = []
		graph_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
		
		t = 0
		while t <= Time:
			#####STORE THE RESULTS FOR EACH ITERATION
			graph_data['timing'].append(t)
			for k in range(0,monitored_gas):
				new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k+1]))#/(2*Inert_pulse_size)
				graph_data['conVtime_'+str(k)].append((new_val))
			new_val = ((to_flux[monitored_gas]*u.vector().get_local()[2*(all_molecules+1)-1]))#/(2*Inert_pulse_size)
			graph_data['conVtime_'+str(all_molecules-1)].append((new_val))
	
			#####Solve for remaining time steps
	
			if t > 0:
				dt = solver_iteration(dt,solver_method_options,solver,dk,1.5,1.1)
			else:
				if reactor_type == 'tap':
					for k in range(0,reactants_num):
						u_n.vector()[int((all_molecules+1)*(mesh_size+1)-1)+k-(all_molecules)] = reactant_ratio[k]*Inert_pulse_size#/(point_volume)
					u_n.vector()[int((all_molecules+1)*(mesh_size+1)-1)] = Inert_pulse_size#/(point_volume)
				
				if k_pulse == 0:
					for z in range(int((cat_location - 0.5*frac_length)*mesh_size)-1,int((cat_location + 0.5*frac_length)*mesh_size)):
						for z_num,z_sur in enumerate(rangesurface_species):
							if int((cat_location - 0.5*frac_length)*mesh_size)-1 <= 0:
								u_n.vector()[(all_molecules+1)-(2+z_num)] = z_sur
							else:
								u_n.vector()[z*(all_molecules+1)-(2+z_num)] = z_sur
	
				dt = solver_iteration(dt,solver_method_options,solver,dk,1.5,1.1)
				
			
			if sensitivity_analysis == True:
				u_final = u.split(deepcopy=False)
				for k_step in range(0,monitored_gas):
					temp = call_sens_analysis(u_final[k_step],controls,dP(1))
					sensitivity_output[k_step].append(temp)
			progressBar(t, Time)
			u_n.assign(u)
			t += dt
		print("")
		total_simulation_time = time.time() - start_time
		print("Complete in "+str(time.time() - start_time)+" seconds")
	
	############# Store the output data #######################################################
		if k_pulse == 0:
			dictionaray_of_numpy_data = {}
			for j_species in range(0,monitored_gas+1):
				print(legend_label[j_species])
				dictionaray_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionaray_of_numpy_data[legend_label[j_species]] = np.vstack((dictionaray_of_numpy_data[legend_label[j_species]],new_data))
	
	############# Visualize/graph the data #######################################################
	
		for k,j in enumerate(graph_data):
			if j != 'timing':
				if k_pulse > 0:
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
				else:
					ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)
				new_data = np.asarray(graph_data[j])
				dictionaray_of_numpy_data[legend_label[k]] = np.vstack((dictionaray_of_numpy_data[legend_label[k]],new_data))
			else:
				pass
	
		if sensitivity_analysis == True:
			for k_sens_step in range(monitored_gas):
				sens_time = np.asarray(graph_data['timing'][0:])
				sens_time = sens_time.T#np.transpose(sens_time)
				sensitivity_output_2 = np.asarray(sensitivity_output[k_sens_step])
				sensitivity_output_2 = np.append(sens_time[:,None],sensitivity_output_2,axis=1)	
				np.savetxt('./'+output_file_name+'_folder/sensitivity_'+output_file_name+'/'+legend_label[k_sens_step]+'/pulse_'+str(k_pulse+1)+'.csv',sensitivity_output_2,delimiter=",",header='t,'+','.join(legend_2))#
	
	ax2.legend(title="Gas Species")
	
	if save_figure == True:
		plt.savefig('./'+output_file_name+'_folder/'+output_file_name+'.png')
	
	if Display_figure == True:
		plt.show()
	
	for j_species in range(0,monitored_gas+1):
		dictionaray_of_numpy_data[legend_label[j_species]] = np.transpose(dictionaray_of_numpy_data[legend_label[j_species]])
		np.savetxt('./'+output_file_name+'_folder/'+legend_label[j_species]+'.csv', dictionaray_of_numpy_data[legend_label[j_species]], delimiter=",")
	
	for k,j in enumerate(graph_data):
		if j != 'timing':
			ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
		else:
			pass

	plt.clf()
	plt.close()

	return total_simulation_time, graph_data, legend_label
