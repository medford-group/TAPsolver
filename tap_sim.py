######## Python packages used / imported ####################################################################
#from pyadjoint import Block
from fenics import *
from fenics_adjoint import *
from variational_form_constructor import make_f_equation
from mpmath import nsum, exp, inf
import matplotlib.pyplot as plt
import fenics_adjoint.types.function_space
from fenics_adjoint.types.compat import assemble_adjoint_value,create_bc,function_from_vector,extract_subfunction,extract_mesh_from_form
#from fenics_adjoint.types.compat import function_space 
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import csv
import sys
import os
import ufl
from pyadjoint.tape import get_working_tape, stop_annotating, annotate_tape,no_annotations
from pyadjoint.enlisting import Enlist
#from pyadjoint.block import evaluate_adj
#from .tape import get_working_tape
#from timestepping import *
from progressbar import ProgressBar

######## Reactor type & data storage ########################################################################

reactor_type = 'tap' 	#currently t_pfr / tap / t_pfr_diff / batch {not ready just yet}
output_file_name = "pt_oxidation" # "./FILENAME.csv"
theta = 1 				# forward_euler = 0, backward_euler = 1, crank_nicolson = 1/2
solver_method_options = 'simple_adaptive'	#'simple_adaptive' or None 
Time = 1				# Total pulse time
Time_steps = 20		# Number of time steps 
mesh_spacing = 100 
mesh_size = 200
sensitivity_analysis = True # True = on / False = off
Frequency_of_sensitivity = 1

############# Feed composition ############################################################################

reactants_num = 1
Inert_pulse_size = 6.022e23	###Size of inert pulse (# of molecules, not mol)#5e-8 #(5e-9)*
reactant_ratio = [1]		###Ratio of reactant to the size of Inert_pulse_size
number_of_pulses = 5

############# Reaction Equations ############################################################################

### Pt oxidation
# actual
#reactions_test = ['O2 + 2* -> 2O*']
# Toy model
reactions_test = ['A + * -> A*']

### CO OXIDATION
#Eley
#reactions_test = ['A + * -> A*','B + A* -> * + C']
#Eley with preadsorbed
#reactions_test = ['CO + O* ->  CO2 + *']
#Lang
#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*']
#Eley & Lang
#reactions_test = ['CO + * <-> CO*','CO* + O* -> CO2 + 2*','CO + O* -> CO2 + *']
#Reece reaction network
#reactions_test = ['CO <-> CO*', 'CO* + OA* -> CO2*', 'CO* + OB* -> CO2*', 'CO2* -> CO2']

### Ammonia decomposition 
#reactions_test = ['NH3 + * <-> NH3*', 'H2 + 2* <-> 2H*', 'N2 + * <-> N2*', 'N2* + * <-> 2N*', 'NH3* + * <-> NH2* + H*', 'NH2* + * <-> NH* + H*', 'NH* + * <-> N* + H*']

############# Reactor Parameters ############################################################################

len_inert_1 = 2.1		###Length of the first inert zone (cm)
len_cat = 0.2			###Length of the catalyst zone (cm) 
len_inert_2 = 2.1		###Length of the second inert zone (cm)
reac_radius = 0.2		###Radius of the reactor

############# Transport properties ##########################################################################

void_inert = 0.4		###Voidage of the catalyst
void_cat = 0.4			###Voidage of the catalyst

ref_rate_inert = 43.5319		###Reference 
ref_rate_cat = 32.5319
ref_mass = 423.15
mass_list = np.array((16,40))
#mass_list = np.array((28,44,40))	
#mass_list = np.array((17.03,2.015,28.01,40))
velocity = 1

############# Catalyst properties ##########################################################################

#Fraction of sites
mcat = 0.025   			###Catalyst mass in grams
ssarea = 5E4			###Surface area of support in cm2/g
bsarea = ssarea*mcat	###Total surface area of support
sdensity = 1.4e15	###Surface Density (atoms/cm2)
atloading = 0.08		### Fractional metal loading of catalyst (assume same surface density as support)
sarea = bsarea*atloading	### Total surface area of metal
open_sites = 0.55*sdensity
OA = 0.15*sdensity 		### OA site density (atoms/cm2)
OB = 0.30*sdensity	

rangesurface_species = [(5e-9)*6.022e23]
#rangesurface_species = [(0.5/1000)*6.022e23]

########## Input Rate Constants #############################################################################

kb = 8.61733e-5
h = 4.1357e-15
T = 400+273.15

# Can redefine the rate constant equation if desired
def rate_con_gen(T,del_g):
	return (kb*T/h)*np.exp(-del_g/(kb*T))

#Ke0 = Constant(rate_con_gen(T,1.12-(-.26)))
#Kd0 = Constant(rate_con_gen(T,1.12))

#Ke1 = Constant(rate_con_gen(T,0.37))
#Kd1 = Constant(rate_con_gen(T,0.37-.3))

#Ke2 = Constant(rate_con_gen(T,0.94))
#Kd2 = Constant(rate_con_gen(T,0.94-0.2))

#Ke3 = Constant(rate_con_gen(T,0.64))
#Kd3 = Constant(rate_con_gen(T,0.64+1.1))

#Ke4 = Constant(rate_con_gen(T,1.36-.39))
#Kd4 = Constant(rate_con_gen(T,1.36))

#Ke5 = Constant(rate_con_gen(T,.56-.39))
#Kd5 = Constant(rate_con_gen(T,.56))

#Ke6 = Constant(rate_con_gen(T,))
#Kd6 = Constant(rate_con_gen(T,))

###### Time stepping/solver options

Ke0 = Constant(1000*(100**3)/(6.022e23))

#Ke0 = Constant(1.15*(100**3)/(6.022e23))
#Kd0 = Constant(10)

#Ke1 = Constant(400*(1000)/(6.022e23))		###Rate constants

#Ke2 = Constant((100**3)/(6.022e23))
#Ke0 = Constant(1.1108e07)
#Kd0 = Constant(1.6900e12)
#Ke1 = Constant(3.354e-7)
#Ke2 = Constant(3.287e-10)
#Ke3 = Constant(9.6077e9)

controls = [Control(Ke0)]
#legend_2 = ['Ke0','Kd0','Ke1','Ke2','Ke3']

#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%##%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%# FEniCS code beyond this point #%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%
#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%##%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#

############# Used to control the output/warnings of Fenics ############################################################################

parameters["std_out_all_processes"] = False																							
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
set_log_active(False)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
tol = 1E-14

if reactor_type == 't_pfr' or 't_pfr_diff' or 'tap':
	
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
	
	def vect_vol(x,y):
		return x*(reac_radius**(2))*np.pi*y
	
	############# Establish grid points in FEniCS #########################################################################################
	
	r_param = np.array((len_inert_1,len_cat,len_inert_2))
	grid_points = round(np.sum(r_param)*mesh_spacing)	###Number of grid points for gasses
	grid_points_2 = round(r_param[1]*mesh_spacing)		###Grid points for catalyst zone
	
	def last_grid_point(x):						###Last grid points for gaseous species
	    return x*mesh_spacing
	
	last_func = np.vectorize(last_grid_point)
	zone_sizes = last_func(r_param)
	grid_loc_1 = np.round(np.cumsum(zone_sizes))
	#!#!@#@#dx_r = np.sum(r_param)/grid_points 			###Distance (cm) per grid point
	dx_r = np.sum(r_param)/mesh_size
	dx2_r = dx_r*dx_r 							###(cm2) for calculation of alpha
	
	### Side work on developing general mesh implementation
	
	frac_length = r_param[1]/(np.sum(r_param))
	cat_location = (r_param[1]*0.5+r_param[0])/(np.sum(r_param))
	
	
	#######################################################
	
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
	
	############# Establish method to define and track reactants/intermediates/products ######################################################
	#
	graph_data = {}
	v_d = {}
	u_d = {}
	u_nd = {}
	sens_data = {}																															
	
	if reactions_test[0] != 'INERT_ONLY':
		inert_st = False
		necessary_values = make_f_equation(reactions_test,reactants_num,reactor_type,inert_st)
		for k in range(0,necessary_values['molecules_in_gas_phase']):
			graph_data['conVtime_'+str(k)] = []
			sens_data['conVtime_'+str(k)] = []
	
		graph_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
		sens_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
	else:
		inert_st = True
		necessary_values = make_f_equation(reactions_test,reactants_num,inert_st)
		necessary_values['molecules_in_gas_phase'] = 0
	
		graph_data['conVtime_'+str(0)] = []
		sens_data['conVtime_'+str(0)] = []
	
	############# Establish the mesh and boundaries/domains in FEniCS #######################################################
			
	#!#!@#@#mesh = UnitIntervalMesh(int(grid_points))
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
	
	#!#!@#@#
	#class thin_zone(SubDomain):
	#	def inside(self, x, on_boundary):
	#		return between(x[0], ( (grid_loc_1[0]/grid_points) , (grid_loc_1[1]/grid_points)))
	#thin_zone = thin_zone()
	
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
	
	def pulse_point(x, on_boundary):
		return near(x[0],1/grid_points,tol)
	
	class integration_section(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (0.96,1.0))
	
	osub = integration_section()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	osub.mark(domains, 1)
	dP = Measure('vertex',domain = mesh, subdomain_data=domains)
	
	############# Establish method to define and track reactants/intermediates/products #######################################################
	
	if reactions_test != ['INERT_ONLY']:
		all_molecules = necessary_values['gas_num']+necessary_values['surf_num']+1
		u = Function(V)
		u_n = Function(V)
		u_temp = Function(V)
		u_1 = Function(V)
		u_2 = Function(V)
		u_3 = Function(V)
		norm_1 = Function(V)
		norm_2 = Function(V)
		tempA = TestFunctions(V)
		tempB = split(u)
		tempC = split(u_n)
	else:
		all_molecules = 1
		u = Function(V)
		u_n = Function(V)
		v = TestFunction(V)
	
	if reactions_test != ['INERT_ONLY']:
		for kit in range(0,(all_molecules)+1):
			v_d['v_'+str(kit+1)] = tempA[kit]
			u_d['u_'+str(kit+1)] = tempB[kit]
			u_nd['u_n'+str(kit+1)] = tempC[kit]
	
	dt = Time/Time_steps
	Length = np.sum(r_param)**2
	
	W = VectorFunctionSpace(mesh, 'P', 1)
	w = Function(W)
	w.vector()[:] = velocity
	
	############# Declare the equations used to describe the system in FEniCS #######################################################
	
	try:
		pass
		F = eval(necessary_values['F'])
	except NameError:
		print("       ")
		print("       ")
		print("NAME ERROR")
		print("       ")
		print("       ")
		print("Likely Need to enter the rate constants.")
		print("There is currently no way to generate or ")
		print("gather the rate constants automatically")
		print("       ")
		print("Must follow the following format")
		print("       ")
		print("'Ke1' = forward rate constant for reaction # 1")
		print("'Kd1' = reverse rate constant for reaction # 1")
		print("       ")
		print("Rate constants for the following equations")
		print("must be included")
		print("       ")
		for j_nu,k_nu in enumerate(reactions_test):
			print("Reaction "+str(j_nu+1)+"     "+k_nu)
		print("       ")
		print("Could also have incorrect diffusion coefficient array defined")
		sys.exit()
	
	#### Just an equations for reeces code (different from others since it requires the use of a scaling parameter ( Q ))
	#Q = Constant(5.5423e3)
	#F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)  + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)  + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0)  + dk*Constant(D[1][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1)  + dk*Constant(D[1][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + ((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + ((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + ((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + ((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + ((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + ((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1) + ((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + ((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1)- dk*(1.0* Q*Kd0*(u_d['u_3']**1.0)*v_d['v_1']*dx(1))+ dk*(1.0* Ke0*(u_d['u_1']**1.0)*v_d['v_1']*dx(1))+ dk*(1.0* Kd0*(u_d['u_3']**1.0)*v_d['v_3']*dx(1))- dk*(1.0* (Ke0/Q)*(u_d['u_1']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_4']*dx(1))- dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_6']*dx(1))- dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke3*(u_d['u_5']**1.0)*v_d['v_5']*dx(1))- dk*(1.0* Q*Ke3*(u_d['u_5']**1.0)*v_d['v_2']*dx(1))+ Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)+ dk * Constant(D[2][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0) + Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1) + dk * Constant(D[2][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)
	
	
	############# Define the boundary conditions (depending on the reactor type) ###############
	
	monitored_gas = necessary_values['molecules_in_gas_phase']
	
	if reactor_type == 'tap':
		bcs = []
		if reactions_test != ['INERT_ONLY']:
			for k in range(0,necessary_values['molecules_in_gas_phase']):
				bcs.append(DirichletBC(V.sub(k),Constant(1),boundary_R))
			bcs.append(DirichletBC(V.sub(all_molecules),Constant(1),boundary_R))
		else:	
			bcs.append(DirichletBC(V,Constant(0),boundary_R))
	
	elif reactor_type == 't_pfr' or 't_pfr_diff':
		bcs = []
		for k in range(0,reactants_num):
			bcs.append(DirichletBC(V.sub(k),Constant(reactant_ratio[k]),boundary_L))
		for kin in range(reactants_num,monitored_gas):
			bcs.append(DirichletBC(V.sub(kin),Constant(0),boundary_L))
		bcs.append(DirichletBC(V.sub(monitored_gas),Constant(0),boundary_L))
	
	############# Initialize the graphs and set the graphical parameters #######################################################
	
	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('$t (s)$')
	if reactor_type == 'tap':
		ax2.set_ylabel('$Flux (molecules/s)$')
	elif reactor_type == 't_pfr' or 't_pfr_diff':
		ax2.set_ylabel('$Outlet Concentration (molecules/cm3)$')
	zef = np.linspace(dx_r, 1.-dx_r, grid_points+1)
	
	legend_label = []
	header = "time"
	if reactions_test[0] != 'INERT_ONLY':
		for k in range(0,necessary_values['molecules_in_gas_phase']):
			legend_label.append(necessary_values['reactants'][k])
			header = header+","+necessary_values['reactants'][k]
	legend_label.append("Inert")
	header = header+",Inert"
	colors = ['b','r','m','g']
	
	ax2.legend(title="Gas Species")
	#ax2.set_xlim(left=0,right=.3)
	
	to_flux = []
	for k in range(0,monitored_gas):
		to_flux.append(2*(1/original_pulse_size) *D[k][0] * dx_r*(reac_radius**2)*3.14159/(eb[0]*dx2_r))
	to_flux.append((2*(1/original_pulse_size) *D[monitored_gas][0] * dx_r*(reac_radius**2)*3.14159/(eb[0]*dx2_r)))
	
	if sensitivity_analysis == True:
		try:
			os.makedirs('./sensitivity_'+output_file_name)
		except OSError:
			pass

	############# Declare the solver to be used #######################################################
	
	J = derivative(F,u)
	problem = NonlinearVariationalProblem(F,u,bcs,J)
	solver = NonlinearVariationalSolver(problem)
	
	solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
	solver.parameters["newton_solver"]["absolute_tolerance"] = 1e3
	
	############# Define a method to iteratively choose time-step (will change based on meeting with tobin isaac) ##
	 
	def solver_iteration(time_step,method):
		if method == 'simple_adaptive':
			try:
				dk.assign(time_step/1.1)
				u_temp.assign(u)
				solver.solve()
				u_1.assign(u)
				u.assign(u_temp)
				dk.assign(time_step*1.1)
				solver.solve()
				u_3.assign(u)
				u.assign(u_temp)
				dk.assign(time_step)
				solver.solve()
				u_2.assign(u)
		
				#u.assign(u_temp)
				#print('')
				#print('norms')
				ref_norm = u_2.vector().norm("l2")
				norm_1 = ((u_1.vector() - u_2.vector()))
				norm1 = norm_1.norm("l2")/ref_norm
				#print(print(norm1))
				norm_2 = (u_2.vector() - u_3.vector())
				norm2 = norm_2.norm("l2")/ref_norm
				#print(print(norm2))
				#print('')
		
				if norm1 > 0.02:
					time_step = time_step/1.1
					dk.assign(time_step)
					u.assign(u_1)
					#print("")
					#print("decrease")
					#print("")
				elif norm2 < 0.02:
					time_step = time_step*1.1
					dk.assign(time_step)
					u.assign(u_3)
					#print("")
					#print("Increase")
				#print("")
				#time.sleep(0.4)
				return time_step
			
			except RuntimeError:
				time_step = time_step*0.5
				print(time_step)
				dk.assign(time_step)
				if time_step < 1e-5:
					print("Time step too low")
					print(time.time() - start_time)
					sys.exit()
				time_step=solver_iteration(time_step,solver_method_options)
				return time_step
		elif method == 'None':
			try:
				solver.solve()
				return time_step
	
			except RuntimeError:
				time_step = time_step*0.5
				print(time_step)
				dk.assign(time_step)
				if time_step < 1e-5:
					print("Time step too low")
					print(time.time() - start_time)
					sys.exit()
				time_step=solver_iteration(time_step,solver_method_options)
				return time_step
	
	
	############# Solve the system #######################################################
	
	for k_pulse in range(0,number_of_pulses):
		dk = Constant(Time/Time_steps)
		dt = Time/Time_steps
		start_time = time.time()
		sensitivity_output = []
		graph_data['timing'] = []
		for k_gasses in range(0,necessary_values['molecules_in_gas_phase']):
			graph_data['conVtime_'+str(k_gasses)] = []
		graph_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
		t = 0
		while t <= Time+0.01:
			#####STORE THE RESULTS FOR EACH ITERATION
			graph_data['timing'].append(t)																									
			if reactions_test != ['INERT_ONLY']:																							
				if reactor_type == 'tap':
					for k in range(0,monitored_gas):
						new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k+1]))#/(2*Inert_pulse_size)
						graph_data['conVtime_'+str(k)].append((new_val))
					new_val = ((to_flux[monitored_gas]*u.vector().get_local()[2*(all_molecules+1)-1]))#/(2*Inert_pulse_size)
					graph_data['conVtime_'+str(all_molecules-1)].append((new_val))
				elif reactor_type == 't_pfr' or 't_pfr_diff':
					for k in range(0,monitored_gas):
						new_val = (( u.vector().get_local()[(all_molecules)+k+1]))
						graph_data['conVtime_'+str(k)].append((new_val))
					new_val = ((u.vector().get_local()[2*(all_molecules+1)-1]))
					graph_data['conVtime_'+str(all_molecules-1)].append((new_val))
			else:
				graph_data['conVtime_0'].append((( (u.vector().get_local()[1]))))
			
			########################################
	
			#####Solve for remaining time steps
	
			if t > 0:			
				dt = solver_iteration(dt,solver_method_options)																											
				newish = u.split(deepcopy=False)
			
			#####Solve and initialize for first time step	
			
			else:
				if reactions_test != ['INERT_ONLY']:
					if reactor_type == 'tap':
						for k in range(0,reactants_num):#necessary_values['reactants_number']
							#!#!@#@#u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)+k-(all_molecules)] = reactant_ratio[k]*Inert_pulse_size#/(point_volume)
							u_n.vector()[int((all_molecules+1)*(mesh_size+1)-1)+k-(all_molecules)] = reactant_ratio[k]*Inert_pulse_size#/(point_volume)
						#!#!@#@#u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)] = Inert_pulse_size#/(point_volume)
						u_n.vector()[int((all_molecules+1)*(mesh_size+1)-1)] = Inert_pulse_size#/(point_volume)
					if k_pulse == 0:
						#!#!@#@#for z in range(int(grid_loc_1[0]-1),int(grid_loc_1[1])):
						if int((cat_location - 0.5*frac_length)*mesh_size)-1 <= 0:
							for z in range(int((cat_location - 0.5*frac_length)*mesh_size)-1,int((cat_location - 0.5*frac_length)*mesh_size)):
								for z_num,z_sur in enumerate(rangesurface_species):
									print(z)
									u_n.vector()[(all_molecules+1)-(2+z_num)] = z_sur

						else:	
							for z in range(int((cat_location - 0.5*frac_length)*mesh_size)-1,int((cat_location - 0.5*frac_length)*mesh_size)):
								for z_num,z_sur in enumerate(rangesurface_species):
									u_n.vector()[z*(all_molecules+1)-(2+z_num)] = z_sur
	
					newish = u_n.split(deepcopy=True)
					dt = solver_iteration(dt,solver_method_options)
					temp_again = u.split(deepcopy=True)
				else:
					u_n.vector()[int(grid_points)-1] = 10000
					solve(F== 0,u,bcs,solver_parameters=parms_solve)
			
			if sensitivity_analysis == True:
				if t < 0.3:
					for k_step in range(0,monitored_gas):
						u_final = u.split(deepcopy=False)
						print(legend_label[k_step]+" Sensitivity Analysis Being Performed")
						sens_func = assemble(inner(u_final[k_step],u_final[k_step])*dP(1))
						X = compute_gradient(sens_func,controls)#compute_gradient(sens_func,[control1,control2])
						new_data = np.asarray(graph_data['timing'])
						sensitivity_output.append(X)
	
			print("Current Time")
			print(t)
			u_n.assign(u)
			t += dt
		
		print(time.time() - start_time)
	
	############# Store the output data #######################################################
		if k_pulse == 0:
			dictionaray_of_numpy_data = {}
			for j_species in range(0,monitored_gas+1):
				dictionaray_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
				#output_data = np.empty(shape=[0, len(graph_data['timing'])])
				new_data = np.asarray(graph_data['timing'])
				dictionaray_of_numpy_data[legend_label[j_species]] = np.vstack((dictionaray_of_numpy_data[legend_label[j_species]],new_data))
	
	############# Visualize/graph the data #######################################################
	
		for k,j in enumerate(graph_data):
			if j != 'timing':
				ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
				new_data = np.asarray(graph_data[j])
				dictionaray_of_numpy_data[legend_label[k-1]] = np.vstack((dictionaray_of_numpy_data[legend_label[k-1]],new_data))
				#output_data = np.vstack((output_data,new_data))
			else:
				pass
		
		if sensitivity_analysis == True:
			np.savetxt('./sensitivity_'+output_file_name+'/pulse_'+str(k_pulse+1)+'.csv',sensitivity_output)
	
	plt.show()
	
	for j_species in range(0,monitored_gas+1):
		dictionaray_of_numpy_data[legend_label[j_species]] = np.transpose(dictionaray_of_numpy_data[legend_label[j_species]])
		np.savetxt('./'+legend_label[j_species]+'_'+output_file_name+'.csv', dictionaray_of_numpy_data[legend_label[j_species]], delimiter=",")
	#output_data = np.transpose(output_data)
	#np.savetxt(output_file_name, output_data, header=header, delimiter=",")
	
	############# Store the parameters used in the model #######################################################
	
	F = open("./"+output_file_name+".txt","w")
	F.write("Reactor Length: "+str(np.sum(r_param))+"\n")
	F.write("First Forward: "+str(1.1108e07)+"\n")
	F.write("Inert pulse size: "+str(Inert_pulse_size)+"\n")
	F.write("rangesurface_species: "+str(rangesurface_species)+"\n")
	F.write("reactions_test: "+str(reactions_test))
	F.close()
	
	
	for k,j in enumerate(graph_data):
		if j != 'timing':
			ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
		else:
			pass
	
	ax2.legend(title="Gas Species")
	
	fig3, ax3 = plt.subplots()
	plt.title("Sensitivity of Species "+legend_label[species_of_interst]+" to Rate Constants")
	ax3.set_xlabel('$t\ (s)$')
	ax3.set_ylabel('$Sensitivity$')
	zef_2 = np.linspace(0,2,sensitivity_output[species_of_interst].shape[0])
	
	for kit in range(0,sensitivity_output[species_of_interst].shape[1]):
		ax3.plot(zef_2,sensitivity_output[species_of_interst][:,kit],label=legend_2[kit], ls = '--', alpha=0.7)
	ax3.legend(title="parameters",loc=7)


	#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%##%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
	#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%# Batch Reactor option  %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
	#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%##%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#

else:

	graph_data = {}
	v_d = {}
	u_d = {}
	u_nd = {}
	sens_data = {}																															
	
	if reactions_test[0] != 'INERT_ONLY':
		inert_st = False
		necessary_values = make_f_equation(reactions_test,reactants_num,reactor_type,inert_st)
		for k in range(0,necessary_values['molecules_in_gas_phase']):
			graph_data['conVtime_'+str(k)] = []
			sens_data['conVtime_'+str(k)] = []
	
		graph_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
		sens_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
	else:
		inert_st = True
		necessary_values = make_f_equation(reactions_test,reactants_num,inert_st)
		necessary_values['molecules_in_gas_phase'] = 0
	
		graph_data['conVtime_'+str(0)] = []
		sens_data['conVtime_'+str(0)] = []
	
	############# Establish the mesh and boundaries/domains in FEniCS #######################################################
			
	#!#!@#@#mesh = UnitIntervalMesh(int(grid_points))
	mesh = UnitIntervalMesh(int(1))
	P1 = FiniteElement('P',mesh.ufl_cell(),1)
	
	if reactions_test != ['INERT_ONLY']:
		test_new = eval(necessary_values['element'])
		element = MixedElement(test_new)
		V = FunctionSpace(mesh,element)
		V_du = FunctionSpace(mesh,P1)
	else:
		V = FunctionSpace(mesh,P1)
		V_du = FunctionSpace(mesh,P1)

	
	#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%##%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
	#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%# Unnecessary code at the moment  %#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
	#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%##%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#%#
	

	#### Alternative Time Step
	
#cur_max = 0
#tot_max = 0

	#if t > 1:
	#if cur_max < (1/500)*tot_max:
	#	dt = 0.05
	#	dk.assign(dt)
	#if test1 < 2:
	#	dt = dt+0.0001
	#	dk.assign(dt)
	#print("terms")
	#print(new_error)
	#print((Old_error/new_error)**0.075)
	#print(((1e-8/new_error)**0.175))
	#print((((Old_error**2)/(new_error*Oldest_error))**0.01))
	#print("")
	#dt = ((Old_error/new_error)**0.075)*((1e-8/new_error)**0.00175)*dt*(((Old_error**2)/(new_error*Oldest_error))**0.01)
	#Oldest_error = Old_error
	#Old_error = new_error

#	max_list=[]
#	max_list.append(new_val)
#	max_list.append(new_val)
#	cur_max = max(max_list)
#	if cur_max > tot_max:
#		tot_max = cur_max
