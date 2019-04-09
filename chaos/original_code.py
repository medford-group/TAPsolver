from fenics import *
from dolfin_adjoint import *
from reece_f_constructor import make_f_equation
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

#dolfin_adjoint.adj_checkpointing(strategy='multistage', steps=11,snaps_on_disk=2, snaps_in_ram=2, verbose=True)

################################################################################################################################	### Used to control the output/warnings of Fenics
parameters["std_out_all_processes"] = False																							
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
set_log_active(False)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

################################################################################################################################	### Choose a set of reactions to investigate
#reactions_test =  ['CO + * <-> CO*', 'CO* + OA* <-> CO2*', 'CO* + OB* <-> CO2*', 'CO2* <-> CO2 + *']								
#reactions_test =  ['CO + * <-> CO*',    'CO* + OA* -> CO2*',    'CO* + OB* -> CO2 + *',    'CO2* -> CO2 + *']
#reactions_test = ['NH3 + * -> NH3*',  'N2 + * -> N2*',   'CO + * -> CO*']
#eactions_test = ['NH3 -> NH3*',  'N2 -> N2*',   'CO -> CO*']
#reactions_test = ['NH3 <-> NH3*',  'N2 <-> N2*',  'H2* <-> 2H*',  'N2* <-> 2N*',  'N* <-> NH* + *',  'NH* <-> NH2* + *',  'NH2* + H* <-> NH3*']
#reactions_test = ['NH3 + * <-> NH3*',  'N2 + * <-> N2*',  'H2 + 2* <-> 2H*',  'N2* + * <-> 2N*',  'N* + H* <-> NH* + *',  'NH* + H* <-> NH2* + *',  'NH2* + H* <-> NH3* + *']

#reactions_test = ['A <-> A*', 'B -> B*', 'A* + B* -> C*', 'A* + B* -> D']
reactions_test =  ['CO <-> CO*', 'CO* + OA* -> CO2*', 'CO* + OB* -> CO2*', 'CO2* -> CO2']
#reactions_test = ['N2 + * <-> N2*']
#reactions_test =  ['CO + * <-> CO*']
#reactions_test = ['INERT_ONLY']
reactants_num = 1																													### Need to set the number of gas reactants being fed into the reactor

def solver_iteration(time_step):
	try:
		test1,test2 = solver.solve()
		return time_step,test1,test2 
	except RuntimeError:
		time_step = time_step*0.5
		dk.assign(time_step)
		#if time_step < 1e-5:
		#	print("Time step too low")
		#	sys.exit()
		time_step,test1,test2=solver_iteration(time_step)
		return time_step,test1,test2
	#
################################################################################################################################	### Different reactor parameters to set
Time = 2
Time_steps = 1000
dk = Constant(Time/Time_steps)
#dk = Time/Time_steps
mass_list = np.array((28,44,40))																									### Mass of each molecular species. Must list inert last
r_param = np.array((2.228,0.206,2.400))																								### Reactor Zone Dimensions

eb = np.array((0.5970,0.6970,0.5970))																								### Reactor Voidage
T = 400+273.15																														### Temperature (K)
h = 4.135667e-15
kb = 8.6173303e-5
R = 8.314
R_or_kb = R
mcat = 0.025   																														### catalyst mass in grams
ssarea = 5E4																														### Surface area of support in cm2/g
bsarea = ssarea*mcat																												### Total surface area of support
sdensity = 1.40E15																													### Surface Density (atoms/cm2)
atloading = 0.08																													### Fractional metal loading of catalyst (assume same surface density as support)
sarea = bsarea*atloading																											### Total surface area of metal
open_sites = 0.55*sdensity
OA = 0.15*sdensity 																													### OA site density (atoms/cm2)
OB = 0.30*sdensity																													### OB site density (atoms/cm2)

psizei = 4.43E15																													### Pulse size Argon (molecules)
psizer = psizei																														### Pulse size CO (1:1 Ar to CO mixture)

ref_mass = 423.15	
ref_rate = np.array((43.5319,32.5319))
ref_rate = np.append(ref_rate,ref_rate[0])                          																### Add inert rate to end of ref_rate
D = np.empty((len(mass_list),3))                          																			### Generate numpy matrix for diffusion coefficients

def diff_func(ref_m,mol_mass,ref_r):     																							### Define diffusion equation (having as a function might not be necessary)
	return ref_r*(mp.sqrt(40*T)/mp.sqrt(423.15*mol_mass))

for k,j in enumerate(mass_list):																									### Use for loop to generate diffusion coefficient
	for k_2,j_2 in enumerate(ref_rate):
		D[k,k_2] = diff_func(ref_mass,j,j_2)

grid_points = round(np.sum(r_param)*100/2)																							### Number of grid points for gasses
grid_points_2 = round(r_param[1]*100/2)																								### Grid points for catalyst zone

def last_grid_point(x):																												### Last grid points for gaseous species
    return x*100/2

last_func = np.vectorize(last_grid_point)
zone_sizes = last_func(r_param)
grid_loc_1 = np.round(np.cumsum(zone_sizes))
surf_points = grid_loc_1[1] - grid_loc_1[0]																							### Number of surface species points

dx_r = np.sum(r_param)/grid_points 																									### Distance (cm) per grid point
dx2_r = dx_r*dx_r 																													### (cm2) for calculation of alpha
reac_radius = 0.2
ca = (reac_radius**2)*3.14159
point_volume = dx_r * ca

def vect_vol(x,y):
	return x*(0.2**(2))*np.pi*y

vect_vol_imp = np.vectorize(vect_vol)

vr = np.sum(vect_vol_imp(r_param,eb))																								### Volume of the reactor
vc = r_param[1]*(0.2**(2))*np.pi*eb[1]																								### Volume of the catalyst zone

Q = sarea/vc 																														### surface to gas scaling coefficient (1/cm)
#print(D)
#sys.exit()
################################################################################################################################	### Initiating FEniCS components
graph_data = {}
v_d = {}
u_d = {}
u_nd = {}
sens_data = {}

tol = 1E-14																															### Set the tolerance in FEniCS (used in defining boundaries)

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

################################################################################################################################	### Define the rate constants
def rate_con_gen(A,Ea):																												### Arrhenius Function
	return A*np.exp((-Ea/(R_or_kb*T)))

A_o = kb*T/h

#A = [A_o,A_o,   A_o,A_o,   A_o,A_o,   A_o,A_o,   A_o,A_o,   A_o,A_o,   A_o,A_o]																		### Pre-exponential factors
A = [1.0000,1.000E15,1.000E-5,1.000E-7,1.000E10]
#Ea = [0.86,1.12,  1.57,0.94,   0.86,0.29,   0.64,1.74,   0.66,0.88,   0.56,0.17,    1.36,0.97]
Ea = [0,35.723e3,19e3,32.0000e3,0.22400e3]
#Ea = [1.14-0.11,1.14,1 ,19*1000,32*1000,0.222*1000]																					### Activation Energy

Ke0 = Constant(80*101325/mp.sqrt(2*np.pi*(1.381E-23)*(0.028/6.022E23)*T)/sdensity/1E4*atloading*A[0])
Kd0 = Constant(rate_con_gen(A[1],Ea[1]))#Constant(rate_con_gen(A[1],Ea[1]))
Ke1 = Constant(rate_con_gen(A[2],Ea[2]))
#Kd1 = Constant(rate_con_gen(A[3],Ea[3]))
Ke2 = Constant(rate_con_gen(A[3],Ea[3]))
#Kd2 = Constant(rate_con_gen(A[5],Ea[5]))
Ke3 = Constant(rate_con_gen(A[4],Ea[4]))
#Kd3 = Constant(rate_con_gen(A[7],Ea[7]))
#Ke4 = Constant(rate_con_gen(A[8],Ea[8]))
#Kd4 = Constant(rate_con_gen(A[9],Ea[9]))
#Ke5 = Constant(rate_con_gen(A[10],Ea[10]))
#Kd5 = Constant(rate_con_gen(A[11],Ea[11]))
#Ke6 = Constant(rate_con_gen(A[12],Ea[12]))
#Kd6 = Constant(rate_con_gen(A[13],Ea[13]))

################################################################################################################################	### Call on function definition script and test for error
try:
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
	sys.exit()
#print(Q)
#print(necessary_values['reactants'])
#print(open_sites)

F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)  + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)  + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + Constant(eb[0])*((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0)  + dk*Constant(D[1][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(0) + Constant(eb[1])*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1)  + dk*Constant(D[1][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_2']), grad(v_d['v_2']))*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + dk*((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1) + ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0) + dk*((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + ((u_d['u_5'] - u_nd['u_n5']))*v_d['v_5']*dx(0) + dk*((u_d['u_5'] - u_nd['u_n5']) )*v_d['v_5']*dx(1) + ((u_d['u_6'] - u_nd['u_n6']))*v_d['v_6']*dx(0) + dk*((u_d['u_6'] - u_nd['u_n6']) )*v_d['v_6']*dx(1) + ((u_d['u_7'] - u_nd['u_n7']))*v_d['v_7']*dx(0) + dk*((u_d['u_7'] - u_nd['u_n7']) )*v_d['v_7']*dx(1)- dk*Constant(Q)*(1.0* Kd0*u_d['u_3']*v_d['v_1']*dx(1))+ dk*(1.0* Ke0*u_d['u_1']*v_d['v_1']*dx(1))+ dk*(1.0* Kd0*u_d['u_3']*v_d['v_3']*dx(1))- dk*(Constant(1/Q))*(1.0* Ke0*u_d['u_1']*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*u_d['u_3']*u_d['u_4']*v_d['v_3']*dx(1))+ dk*(1.0* Ke1*u_d['u_3']*u_d['u_4']*v_d['v_4']*dx(1))- dk*(1.0* Ke1*u_d['u_3']*u_d['u_4']*v_d['v_5']*dx(1))+ dk*(1.0* Ke2*u_d['u_3']*u_d['u_6']*v_d['v_3']*dx(1))+ dk*(1.0* Ke2*u_d['u_3']*u_d['u_6']*v_d['v_6']*dx(1))- dk*(1.0* Ke2*u_d['u_3']*u_d['u_6']*v_d['v_5']*dx(1))+ dk*(1.0* Ke3*u_d['u_5']*v_d['v_5']*dx(1))- dk*Constant(Q)*(1.0* Ke3*u_d['u_5']*v_d['v_2']*dx(1))+ Constant(eb[0]) * ((u_d['u_8'] - u_nd['u_n8']))*v_d['v_8']*dx(0)+ dk * Constant(D[2][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(0) + Constant(eb[1]) * ((u_d['u_8'] - u_nd['u_n8']) )*v_d['v_8']*dx(1) + dk * Constant(D[2][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_8']), grad(v_d['v_8']))*dx(1)
#F = Constant(eb[0])*((u_d['u_1'] - u_nd['u_n1']))*v_d['v_1']*dx(0)  + dk*Constant(D[0][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(0) + Constant(eb[1])*((u_d['u_1'] - u_nd['u_n1']) )*v_d['v_1']*dx(1)  + dk*Constant(D[0][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_1']), grad(v_d['v_1']))*dx(1) + ((u_d['u_2'] - u_nd['u_n2']))*v_d['v_2']*dx(0) + dk*((u_d['u_2'] - u_nd['u_n2']) )*v_d['v_2']*dx(1) + ((u_d['u_3'] - u_nd['u_n3']))*v_d['v_3']*dx(0) + dk*((u_d['u_3'] - u_nd['u_n3']) )*v_d['v_3']*dx(1)+ dk*(1.0* Ke0*u_d['u_1']*u_d['u_3']*v_d['v_1']*dx(1))+ dk*(1.0* Ke0*u_d['u_1']*u_d['u_3']*v_d['v_3']*dx(1))- dk*Constant(1/Q)*(1.0* Ke0*u_d['u_1']*u_d['u_3']*v_d['v_2']*dx(1))+ Constant(eb[0]) * ((u_d['u_4'] - u_nd['u_n4']))*v_d['v_4']*dx(0)+ dk * Constant(D[1][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_4']), grad(v_d['v_4']))*dx(0) + Constant(eb[1]) * ((u_d['u_4'] - u_nd['u_n4']) )*v_d['v_4']*dx(1) + dk * Constant(D[1][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_4']), grad(v_d['v_4']))*dx(1)
################################################################################################################################	### Define the zero concentration BC for each gas
bcs = []
if reactions_test != ['INERT_ONLY']:
	for k in range(0,necessary_values['molecules_in_gas_phase']):
		bcs.append(DirichletBC(V.sub(k),Constant(0),boundary_R))
	bcs.append(DirichletBC(V.sub(all_molecules),Constant(0),boundary_R))
else:	
	bcs.append(DirichletBC(V,Constant(0),boundary_R))

################################################################################################################################	### Define parameters / graph used in next sections
fig2, ax2 = plt.subplots()
ax2.set_xlabel('$t (s)$')
ax2.set_ylabel('$Flux (molecules/s)$')
zef = np.linspace(dx_r, 1.-dx_r, grid_points+1)

to_flux = []
monitored_gas = necessary_values['molecules_in_gas_phase']

for k in range(0,monitored_gas):
	to_flux.append(2 *D[k][0] * dx_r*(0.2**2)*3.14159/(eb[0]*dx2_r))
to_flux.append((2* D[monitored_gas][0] * dx_r*(0.2**2)*3.14159/(eb[0]*dx2_r)))

################################################################################################################################	### Solve the system
t = 0

J = derivative(F,u)
problem = NonlinearVariationalProblem(F,u,bcs,J)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8

cur_max = 0
tot_max = 0

start_time = time.time()

while t <= Time+0.01:	
	graph_data['timing'].append(t)																									### Store each time step
	max_list=[]
	if reactions_test != ['INERT_ONLY']:																							### Store the outlet flux for each observed gas
		for k in range(0,monitored_gas):
			new_val = (to_flux[k]*( u.vector().get_local()[(all_molecules)+k+1]))
			max_list.append(new_val)
			graph_data['conVtime_'+str(k)].append((new_val))
		graph_data['conVtime_'+str(all_molecules-1)].append(((to_flux[monitored_gas]*u.vector().get_local()[2*(all_molecules+1)-1])))
	
	else:
		graph_data['conVtime_0'].append((( (u.vector().get_local()[1]))))
	#print(F)
	#time.sleep(2)
	cur_max = max(max_list)
	if cur_max > tot_max:
		tot_max = cur_max

	if t > 0:																														###
		#print(type(dt))
		#print("TEST")
		#time.sleep(3)
		dt,test1,test2 = solver_iteration(dt)
		#test1,test2 = solver.solve()
		#time.sleep(.2)
		#u0.assign(u)
		#data = solve(F == 0,u,bcs,solver_parameters=parms_solve)
		newish = u.split(deepcopy=True)
		#ax2.plot(zef,newish[2].vector().get_local(), ls = '--', alpha=0.7)
		#sys.exit()
	else:
		if reactions_test != ['INERT_ONLY']:
			for z in range(int(grid_loc_1[0]-1),int(grid_loc_1[1])):
				#u_n.vector()[z*(all_molecules+1)-2] = open_sites
				u_n.vector()[z*(all_molecules+1)-(8-5)] = OB
				u_n.vector()[z*(all_molecules+1)-(8-3)] = OA
			
			
			for k in range(0,reactants_num):#necessary_values['reactants_number']
				u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)+k-(all_molecules)] = psizei/(point_volume)
			u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)] = psizei/(point_volume)
			newish = u_n.split(deepcopy=True)
			dt,test1,test2 = solver_iteration(dt)
			#solve(F == 0,u,bcs,solver_parameters=parms_solve)
			temp_again = u.split(deepcopy=True)
			#test1 = 1
			#test2 = 1
		else:
			u_n.vector()[int(grid_points)-1] = 10000
			solve(F== 0,u,bcs,solver_parameters=parms_solve)
	
	print("TIME STEP")
	print(t)
	u_n.assign(u)
	if cur_max < (1/500)*tot_max:
		dt = 0.05
		dk.assign(dt)
	if test1 < 2:
		dt = dt+0.0001
		dk.assign(dt)
	t += dt
print(time.time() - start_time)
################################################################################################################################	### Generate the Legend & Color options for output graph
legend_label = []
if reactions_test[0] != 'INERT_ONLY':
	for k in range(0,necessary_values['molecules_in_gas_phase']):
		legend_label.append(necessary_values['reactants'][k])
legend_label.append("Inert")
colors = ['b','r','m','g','o']

################################################################################################################################	### Graph the output Curves
for k,j in enumerate(graph_data):
	if j != 'timing':
		#print(legend_label[k-1])
		ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
		#print(j)
	else:
		pass
#ax2.set_xlim(left=0.4,right=0.6)
ax2.legend(title="Gas Species")
plt.show()
sys.exit()