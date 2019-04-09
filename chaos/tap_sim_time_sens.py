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

######################################################################################################

def evaluate(tape,step):
	if type(tape._blocks[step]) == fenics_adjoint.solving.SolveBlock:
		fwd_block_variable = tape._blocks[step].get_outputs()[0]
		#Takes u2 (the solution)
		u2 = fwd_block_variable.output
		#Extracts the function spece (or mesh I think)
		V = u2.function_space()
		#Creates an empty function based on the function space
		adj_var = Function(V)
		#Store the previous dJdu value (it is a general matrix and should be easy to do math on)
		dJdu = fwd_block_variable.adj_value

		replaced_coeffs = {}
		#Same value as u2 above
		tmp_u = tape._blocks[step].func
		#The weak form the of pde
		F_form = tape._blocks[step].lhs

		#Get dependenicies will list all boundaries and constants for each block
		for block_variable in tape._blocks[step].get_dependencies():
			coeff = block_variable.output

			#This if statement will only change the 
			if coeff in F_form.coefficients():
				replaced_coeffs[coeff] = block_variable.saved_output
		replaced_coeffs[tmp_u] = fwd_block_variable.saved_output
		#replace the coefficients for the ufl form
		F_form = backend.replace(F_form, replaced_coeffs)

		dFdu = backend.derivative(F_form, fwd_block_variable.saved_output, backend.TrialFunction(u.function_space()))
		dFdu_form = backend.adjoint(dFdu)
		dFdu = assemble_adjoint_value(dFdu_form, **tape._blocks[step].assemble_kwargs)
		
		dJdu_copy = dJdu.copy()
		test_old = dFdu
		bcs = []
		for bc in tape._blocks[step].bcs:
			if isinstance(bc, backend.DirichletBC):
				bc = create_bc(bc, homogenize=True)
		bcs.append(bc)
		bc.apply(dFdu, dJdu)
		#print(type(dFdu))
		#print((adj_var.vector()))
		#print(dJdu.get_local)
		#print(type(dJdu))
		#print(dJdu.get_local())
		#linear solve?
		backend.solve(dFdu, adj_var.vector(), dJdu)
		#print((adj_var.vector()))
		#sys.exit()
		#Only used in the boundary condition part
		adj_var_bdy = function_from_vector(V, dJdu_copy - assemble_adjoint_value(backend.action(dFdu_form, adj_var)))
		print("####################################################")
		for block_variable in tape._blocks[step].get_dependencies():
			c = block_variable.output
			if c != tape._blocks[step].func or tape._blocks[step].linear:
				c_rep = replaced_coeffs.get(c, c)

				if isinstance(c, backend.Function):
					dFdm = -backend.derivative(F_form, c_rep, backend.TrialFunction(c.function_space()))
					dFdm = backend.adjoint(dFdm)
					dFdm = dFdm*adj_var
					dFdm = assemble_adjoint_value(dFdm, **tape._blocks[step].assemble_kwargs)
					#print(type(dFdm))
					#sys.exit()
					block_variable.add_adj_output(dFdm)
				elif isinstance(c, backend.Constant):
					mesh = extract_mesh_from_form(F_form)
					dFdm = -backend.derivative(F_form, c_rep, backend.TrialFunction(c._ad_function_space(mesh)))
					dFdm = backend.adjoint(dFdm)
					dFdm = dFdm*adj_var
					dFdm = assemble_adjoint_value(dFdm, **tape._blocks[step].assemble_kwargs)
					print("constant")
					print((dFdm.get_local()))
					#time.sleep(.5)
					block_variable.add_adj_output(dFdm)
				elif isinstance(c, backend.DirichletBC):
					tmp_bc = create_bc(c, value=extract_subfunction(adj_var_bdy, c.function_space()))
					block_variable.add_adj_output([tmp_bc])
				elif isinstance(c, backend.Expression):
					mesh = F_form.ufl_domain().ufl_cargo()
					c_fs = c._ad_function_space(mesh)
					dFdm = -backend.derivative(F_form, c_rep, backend.TrialFunction(c_fs))
					dFdm = backend.adjoint(dFdm)
					dFdm = dFdm * adj_var
					dFdm = assemble_adjoint_value(dFdm, **tape._blocks[step].assemble_kwargs)
					block_variable.add_adj_output([[dFdm, c_fs]])
		#time.sleep(1)

	else:
		tape._blocks[step].evaluate_adj()
	
	return tape

def compute_gradient(J,m,block_idx=0,options=None,tape=None):
	pbar = ProgressBar()
	X = np.empty(shape=[0, len(m)])
	options = {} if options is None else options
	tape = get_working_tape() if tape is None else tape
	tape.reset_variables()
	J.adj_value = 1.0
	
	## MY METHOD FOR CALCULATING
	#for k in range(len(tape.get_blocks())-1,0,-1):#pbar(range(len(tape.get_blocks())-1,0,-1)):
	#	with stop_annotating():
	#		tape = evaluate(tape,k)
	#	#print("TEST")
	print((m[0].values()[0]))
	m = Enlist(m)
	print((m[0].get_derivative(options=options).values()))
	sys.exit()
	#if k < len(tape.get_blocks())-(len(bcs)) and k%(len(bcs)+2) == 0:
	grads = [i.get_derivative(options=options) for i in m]
	print(type(m[0]))
	value_list = []
	for nu in grads:
		value_list.append(nu.values()[0])
	temp = np.array((value_list))
	X = np.vstack((X,temp))
	#pbar.flush()
	return X


#########################################################################################################3
### INPUTS
#O2 irr
#reactions_test = ['O2 + 2* -> 2O*']
#Eley
#reactions_test = ['A + * -> A*','B + A* -> * + C']
#Eley with preadsorbed
#reactions_test = ['CO + O* ->  CO2 + *']
#Lang
#reactions_test = ['CO + * -> CO*','CO* + O* -> CO2 + 2*']

#Reece reaction network
reactions_test = ['CO <-> CO*', 'CO* + OA* -> CO2*', 'CO* + OB* -> CO2*', 'CO2* -> CO2']

#Original reaction 

reactants_num = 1
Inert_pulse_size = 4.43e15	###Size of inert pulse (# of molecules, not mol)#5e-8
reactant_ratio = [1]		###Ratio of reactant to the size of Inert_pulse_size

#### TIME 
Time = .06		###Total pulse time
Time_steps = 10	###Number of time steps 

#### Setting void & Length 
T = 400+273.15			###Reactor Temperature
reac_radius = 0.2		###Radius of the reactor
void_inert = 0.5970		###Voidage of the catalyst
void_cat = 0.6970			###Voidage of the catalyst

len_inert_1 = 2.228		###Length of the first inert zone (cm)
len_cat = 0.206			###Length of the catalyst zone (cm) 
len_inert_2 = 2.4		###Length of the second inert zone (cm)

ref_rate_inert = 43.5319		###Reference 
ref_rate_cat = 32.5319
ref_mass = 423.15
#mass_list = np.array((28,40))
mass_list = np.array((28,44,40))	

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
rangesurface_species = [0,OB,0,OA,0]

def rate_con_gen(A,Ea):
	return A*np.exp((-Ea/(R_or_kb*T)))

Ke0 = Constant(1.1108e7)#Constant(100*(100**3)/(6.022e23))		###Rate constants
Kd0 = Constant(1.6900e12)
Ke1 = Constant(3.354e-7)
Ke2 = Constant(3.287e-10)
Ke3 = Constant(9.6077e9)

mesh_spacing = 100/2 


######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
#Options for reactions

#reactions_test =  ['CO + * <-> CO*', 'CO* + OA* <-> CO2*', 'CO* + OB* <-> CO2*', 'CO2* <-> CO2 + *']								
#reactions_test =  ['CO + * <-> CO*',    'CO* + OA* -> CO2*',    'CO* + OB* -> CO2 + *',    'CO2* -> CO2 + *']
#reactions_test = ['NH3 + * -> NH3*',  'N2 + * -> N2*',   'CO + * -> CO*']
#eactions_test = ['NH3 -> NH3*',  'N2 -> N2*',   'CO -> CO*']
#reactions_test = ['NH3 <-> NH3*',  'N2 <-> N2*',  'H2* <-> 2H*',  'N2* <-> 2N*',  'N* <-> NH* + *',  'NH* <-> NH2* + *',  'NH2* + H* <-> NH3*']
#reactions_test = ['NH3 + * <-> NH3*',  'N2 + * <-> N2*',  'H2 + 2* <-> 2H*',  'N2* + * <-> 2N*',  'N* + H* <-> NH* + *',  'NH* + H* <-> NH2* + *',  'NH2* + H* <-> NH3* + *']
#reactions_test = ['A <-> A*', 'B -> B*', 'A* + B* -> C*', 'A* + B* -> D']
#reactions_test =  ['CO <-> CO*', 'CO* + OA* -> CO2*', 'CO* + OB* -> CO2*', 'CO2* -> CO2']
#reactions_test = ['N2 + * <-> N2*']
#reactions_test =  ['CO + * <-> CO*']
#reactions_test =  ['CO <-> CO*', 'CO* + OA* -> CO2*', 'CO* + OB* -> CO2*', 'CO2* -> CO2']
#reactions_test = ['INERT_ONLY']

######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################
######################################################################################################################################################################################

#dolfin_adjoint.adj_checkpointing(strategy='multistage', steps=11,snaps_on_disk=2, snaps_in_ram=2, verbose=True)

###################################################	### Used to control the output/warnings of Fenics
parameters["std_out_all_processes"] = False																							
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
set_log_active(False)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

###################################################	### Different reactor parameters to set

dk = Constant(Time/Time_steps)
eb = np.array((void_inert,void_cat,void_inert))

r_param = np.array((len_inert_1,len_cat,len_inert_2))

grid_points = round(np.sum(r_param)*mesh_spacing)	###Number of grid points for gasses
grid_points_2 = round(r_param[1]*mesh_spacing)		###Grid points for catalyst zone

def last_grid_point(x):						###Last grid points for gaseous species
	return x*mesh_spacing

last_func = np.vectorize(last_grid_point)
zone_sizes = last_func(r_param)
grid_loc_1 = np.round(np.cumsum(zone_sizes))
dx_r = np.sum(r_param)/grid_points 			###Distance (cm) per grid point
dx2_r = dx_r*dx_r 							###(cm2) for calculation of alpha
ca = (reac_radius**2)*3.14159				###Cross sectional area of reactor
point_volume = dx_r * ca 					###Volume at a specific point in the reactor
specific_area = 1e8
surface_area_at_point = point_volume*specific_area
open_sites_per_point = sdensity*surface_area_at_point

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

print(D)

vect_vol_imp = np.vectorize(vect_vol)
vr = np.sum(vect_vol_imp(r_param,eb))		### Volume of the reactor
vc = r_param[1]*(reac_radius**(2))*np.pi*eb[1]		### Volume of the catalyst zone
Q = sarea/vc 			### surface to gas scaling coefficient (1/cm)
Inert_pulse_size = Inert_pulse_size/(dx_r*(reac_radius**2)*3.14159)
#print(dx_r)
#print(reac_radius)
###################################################	### Initiating FEniCS components
graph_data = {}
v_d = {}
u_d = {}
u_nd = {}
sens_data = {}

tol = 1E-14																															

graph_data['timing'] = []

if reactions_test[0] != 'INERT_ONLY':
	inert_st = False
	necessary_values = make_f_equation(reactions_test,reactants_num,inert_st)
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
		
mesh = UnitIntervalMesh(int(grid_points))
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
		return between(x[0], ( (grid_loc_1[0]/grid_points) , (grid_loc_1[1]/grid_points)))
thin_zone = thin_zone()

domains = CellFunction("size_t", mesh)
domains.set_all(0)
thin_zone.mark(domains,1)

dx = Measure("dx")[domains]

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0],0,tol)

def boundary_R(x, on_boundary):
	return on_boundary and near(x[0],1,tol)

def pulse_point(x, on_boundary):
	return near(x[0],1/grid_points,tol)

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
Q = Constant(5.5423e3)

###########################################################	### Call on function definition script and test for parameter input errors
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


F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)  + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)  + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0)  + dk*Constant(D[1][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1)  + dk*Constant(D[1][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + ((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + ((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + ((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + ((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + ((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + ((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1) + ((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + ((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1)- dk*(1.0* Q*Kd0*(u_d['u_3']**1.0)*v_d['v_1']*dx(1))+ dk*(1.0* Ke0*(u_d['u_1']**1.0)*v_d['v_1']*dx(1))+ dk*(1.0* Kd0*(u_d['u_3']**1.0)*v_d['v_3']*dx(1))- dk*(1.0* (Ke0/Q)*(u_d['u_1']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_4']*dx(1))- dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_3']*dx(1))+ dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_6']*dx(1))- dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_5']*dx(1))+ dk*(1.0* Ke3*(u_d['u_5']**1.0)*v_d['v_5']*dx(1))- dk*(1.0* Q*Ke3*(u_d['u_5']**1.0)*v_d['v_2']*dx(1))+ Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)+ dk * Constant(D[2][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0) + Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1) + dk * Constant(D[2][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)

######CAN ERASE########

#Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)  + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + 
#Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)  + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + 
#Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0)  + dk*Constant(D[1][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + 
#Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1)  + dk*Constant(D[1][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + 

#((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + 
#((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + 
#((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + 
#((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + 
#((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + 
#((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + 
#((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + 
#((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1) + 
#((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + 
#((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1)- 

#dk*(1.0* Q*Kd0*(u_d['u_3']**1.0)*v_d['v_1']*dx(1))+ 
#dk*(1.0* Ke0*(u_d['u_1']**1.0)*v_d['v_1']*dx(1))+ 

#dk*(1.0* Kd0*(u_d['u_3']**1.0)*v_d['v_3']*dx(1))- 
#dk*(1.0* (Ke0/Q)*(u_d['u_1']**1.0)*v_d['v_3']*dx(1))+ 

#dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_3']*dx(1))+
#dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_4']*dx(1))-
#dk*(1.0* Ke1*(u_d['u_3']**1.0)*(u_d['u_4']**1.0)*v_d['v_5']*dx(1))+

#dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_3']*dx(1))+ 
#dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_6']*dx(1))- 
#dk*(1.0* Ke2*(u_d['u_3']**1.0)*(u_d['u_6']**1.0)*v_d['v_5']*dx(1))+ 

#dk*(1.0* Ke3*(u_d['u_5']**1.0)*v_d['v_5']*dx(1))- 
#dk*(1.0* Q*Ke3*(u_d['u_5']**1.0)*v_d['v_2']*dx(1))+ 

#Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)+ dk * Constant(D[2][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0) +
#Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1) + dk * Constant(D[2][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)

#######################


#######################################################	### Define the zero concentration BC for each gas
bcs = []
if reactions_test != ['INERT_ONLY']:
	for k in range(0,necessary_values['molecules_in_gas_phase']):
		bcs.append(DirichletBC(V.sub(k),Constant(1),boundary_R))
	bcs.append(DirichletBC(V.sub(all_molecules),Constant(1),boundary_R))
else:	
	bcs.append(DirichletBC(V,Constant(0),boundary_R))

###################################################	### Define parameters / graph used in next sections
fig2, ax2 = plt.subplots()
ax2.set_xlabel('$t (s)$')
ax2.set_ylabel('$Flux (molecules/s)$')
zef = np.linspace(dx_r, 1.-dx_r, grid_points+1)

to_flux = []
monitored_gas = necessary_values['molecules_in_gas_phase']

for k in range(0,monitored_gas):
	to_flux.append(2 *D[k][0] * dx_r*(reac_radius**2)*3.14159/(eb[0]*dx2_r))
to_flux.append((2* D[monitored_gas][0] * dx_r*(reac_radius**2)*3.14159/(eb[0]*dx2_r)))
###################################################	### Setup Solver
J = derivative(F,u)
problem = NonlinearVariationalProblem(F,u,bcs,J)
solver = NonlinearVariationalSolver(problem)

solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e0

###################################################	### Used as a make shift adaptive time step 
def solver_iteration(time_step):
	try:
		##dk.assign(time_step/1.1)
		##u_temp.assign(u)
		solver.solve()
		##u_1.assign(u)
		##u.assign(u_temp)
		##dk.assign(time_step*1.1)
		##solver.solve()
		##u_3.assign(u)
		##u.assign(u_temp)
		##dk.assign(time_step)
		##solver.solve()
		##u_2.assign(u)

		#u.assign(u_temp)
		#print('')
		#print('norms')
		##ref_norm = u_2.vector().norm("l2")
		##norm_1 = ((u_1.vector() - u_2.vector()))
		##norm1 = norm_1.norm("l2")/ref_norm
		#print(print(norm1))
		##norm_2 = (u_2.vector() - u_3.vector())
		##norm2 = norm_2.norm("l2")/ref_norm
		#print(print(norm2))
		#print('')

		##if norm1 > 0.02:
		##	time_step = time_step/1.5
		##	dk.assign(time_step)
		##	u.assign(u_1)
			#print("")
			#print("decrease")
			#print("")
		##elif norm2 < 0.02:
		##	time_step = time_step*1.1
		##	dk.assign(time_step)
		##	u.assign(u_3)
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
		time_step=solver_iteration(time_step)
		return time_step
		
###################################################	### Solve the nonlinear variational problem 

class integration_section(SubDomain):
	def inside(self, x, on_boundary):
		return between(x[0], (0.96,1.0))

osub = integration_section()
domains = MeshFunction("size_t", mesh,0)
domains.set_all(0)
osub.mark(domains, 1)
dP = Measure('vertex',domain = mesh, subdomain_data=domains)


controls = [Control(Ke0),Control(Kd0),Control(Ke1),Control(Ke2),Control(Ke3)]

legend_2 = ['Ke0','Kd0','Ke1','Ke2','Ke3']

sensitivity_output = []

start_time = time.time()
t = 0
while t <= Time+0.01:	
	graph_data['timing'].append(t)																									### Store each time step

	if reactions_test != ['INERT_ONLY']:																							### Store the outlet flux for each observed gas
		for k in range(0,monitored_gas):
			new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k+1]))#/(2*Inert_pulse_size)
			graph_data['conVtime_'+str(k)].append((new_val))
		new_val = ((to_flux[monitored_gas]*u.vector().get_local()[2*(all_molecules+1)-1]))#/(2*Inert_pulse_size)
		graph_data['conVtime_'+str(all_molecules-1)].append((new_val))
		
	else:
		graph_data['conVtime_0'].append((( (u.vector().get_local()[1]))))

	if t > 0:			
		dt = solver_iteration(dt)																											
		newish = u.split(deepcopy=False)
		u_n.vector()[int(grid_points)-1] = 10000
		u_final = u.split(deepcopy=False)
		#print(legend_label[k_step]+" Sensitivity Analysis Being Performed")
		sens_func = assemble(inner(u_final[0],u_final[0])*dP(1))
		X = compute_gradient(sens_func,controls)#compute_gradient(sens_func,[control1,control2])
		print("output X")
		print(X)

		#ax2.plot(zef,newish[2].vector().get_local(), ls = '--', alpha=0.7)
		
	else:
		if reactions_test != ['INERT_ONLY']:
			for z in range(int(grid_loc_1[0]-1),int(grid_loc_1[1])):
				for z_num,z_sur in enumerate(rangesurface_species):
					u_n.vector()[z*(all_molecules+1)-(2+z_num)] = z_sur
			for k in range(0,reactants_num):#necessary_values['reactants_number']
				u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)+k-(all_molecules)] = reactant_ratio[k]*Inert_pulse_size#/(point_volume)
			u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)] = Inert_pulse_size#/(point_volume)
			newish = u_n.split(deepcopy=True)
			dt = solver_iteration(dt)
			temp_again = u.split(deepcopy=True)
		else:
			u_n.vector()[int(grid_points)-1] = 10000
			solve(F== 0,u,bcs,solver_parameters=parms_solve)
			u_final = u.split(deepcopy=False)
			#print(legend_label[k_step]+" Sensitivity Analysis Being Performed")
			sens_func = assemble(inner(u_final[0],u_final[0])*dP(1))
			X = compute_gradient(sens_func,controls)#compute_gradient(sens_func,[control1,control2])
	#sensitivity_output.append(X)


	print("Current Time")
	print(t)
	u_n.assign(u)
	t += dt

print(time.time() - start_time)

###########################################	### Generate the Legend & Color options for output graph
legend_label = []
if reactions_test[0] != 'INERT_ONLY':
	for k in range(0,necessary_values['molecules_in_gas_phase']):
		legend_label.append(necessary_values['reactants'][k])
legend_label.append("Inert")
colors = ['b','r','m','g','o']

###########################################	### Graph the output Curves
for k,j in enumerate(graph_data):
	if j != 'timing':
		ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
	else:
		pass


#print(legend_label)
#ax2.plot(graph_data['timing'],graph_data['conVtime_1'],color=colors[0],label=legend_label[1], ls = '--', alpha=0.7)
#ax2.set_xlim(left=0,right=.3)
ax2.legend(title="Gas Species")

##########################################	### Sensitivity analysis


####################################################################



for k_step in range(0,monitored_gas):
	print(legend_label[k_step]+" Sensitivity Analysis Being Performed")
	sens_func = assemble(inner(u_final[k_step],u_final[k_step])*dP(1))
	X = compute_gradient(sens_func,controls)#compute_gradient(sens_func,[control1,control2])
	#sensitivity_output.append(X)
	print(X)
sys.exit()
###############################

################################################################################################################################	### Generate the Legend & Color options for output graph


################################################################################################################################	### Graph the output Curves
for k,j in enumerate(graph_data):
	if j != 'timing':
		ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
	else:
		pass

ax2.legend(title="Gas Species")

#print(len(graph_data['timing']))
#print(sensitivity_output[0].shape[0])
#sys.exit()
###############################

species_of_interst = 0

fig3, ax3 = plt.subplots()
plt.title("Sensitivity of Species "+legend_label[species_of_interst]+" to Rate Constants")
ax3.set_xlabel('$t\ (s)$')
ax3.set_ylabel('$Sensitivity$')
zef_2 = np.linspace(0,2,sensitivity_output[species_of_interst].shape[0])

for kit in range(0,sensitivity_output[species_of_interst].shape[1]):
	ax3.plot(zef_2,sensitivity_output[species_of_interst][:,kit],label=legend_2[kit], ls = '--', alpha=0.7)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax3.legend(title="parameters",loc=7)

species_of_interst = 1

fig4, ax4 = plt.subplots()
plt.title("Sensitivity of Species "+legend_label[species_of_interst]+" to Rate Constants")
ax4.set_xlabel('$t\ (s)$')
ax4.set_ylabel('$Sensitivity$')
zef_3 = np.linspace(0,2,sensitivity_output[species_of_interst].shape[0])

for kit in range(0,sensitivity_output[species_of_interst].shape[1]):
	ax4.plot(zef_3,sensitivity_output[species_of_interst][:,kit],label=legend_2[kit], ls = '--', alpha=0.7)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax4.legend(title="parameters",loc=7)
sys.exit()
plt.show()

sys.exit()

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
