from fenics import *
from fenics_adjoint import *
from pfr_variational import make_f_equation
from mpmath import nsum, exp, inf
#import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import csv
import sys
import os

parameters["std_out_all_processes"] = False																							
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
#set_log_active(False)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

reactions_test = ['O2 <-> 2O','CO + O -> CO2']

#reactions_test = ['A <-> B']

#reactions_test = ['A <-> B','B <-> C','A + B <-> D', 'C + D <-> E']

reactants_num = 1
Inert_pulse_size = (2.5e-8)*(6.022e23)																											### Size of inert pulse (# of molecules, not mol)
reactant_ratio = [1,1]

#### TIME 
Time = 0.6
Time_steps = 50

T = 400+273.15
reac_radius = 0.2

len_inert_1 = 2
len_cat = 0.1
len_inert_2 = 2

Ke0 = Constant(5000)
Kd0 = Constant(200)

Ke1 = Constant(20)
Kd1 = Constant(1000)

Ke2= Constant(230)
Kd2 = Constant(400*(100**3)/(6.022e23))

Ke3 = Constant(1000)
Kd3 = Constant(2000)


#########################

#########################

dk = Constant(Time/Time_steps)

r_param = np.array((len_inert_1,len_cat,len_inert_2))

grid_points = round(np.sum(r_param)*100/2)																							### Number of grid points for gasses
grid_points_2 = round(r_param[1]*100/2)																								### Grid points for catalyst zone

def last_grid_point(x):																												### Last grid points for gaseous species
    return x*100/2

last_func = np.vectorize(last_grid_point)
zone_sizes = last_func(r_param)
grid_loc_1 = np.round(np.cumsum(zone_sizes))
dx_r = np.sum(r_param)/grid_points 																									### Distance (cm) per grid point
dx2_r = dx_r*dx_r 																													### (cm2) for calculation of alpha
ca = (reac_radius**2)*3.14159																										### cross sectional area of reactor
point_volume = dx_r * ca 																											### volume at a specific point in the reactor

graph_data = {}
v_d = {}
u_d = {}
u_nd = {}
sens_data = {}

tol = 1E-14																															

graph_data['timing'] = []

inert_st = False
necessary_values = make_f_equation(reactions_test,reactants_num,inert_st)
for k in range(0,necessary_values['molecules_in_gas_phase']):
	graph_data['conVtime_'+str(k)] = []
	sens_data['conVtime_'+str(k)] = []
graph_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []
sens_data['conVtime_'+str(necessary_values['gas_num']+necessary_values['surf_num'])] = []

mesh = UnitIntervalMesh(int(grid_points))
P1 = FiniteElement('P',mesh.ufl_cell(),1)

test_new = eval(necessary_values['element'])
element = MixedElement(test_new)
V = FunctionSpace(mesh,element)
V_du = FunctionSpace(mesh,P1)

class thin_zone(SubDomain):
	def inside(self, x, on_boundary):
		return between(x[0], ( (grid_loc_1[0]/grid_points) , (grid_loc_1[1]/grid_points)))
thin_zone = thin_zone()

domains = MeshFunction("size_t", mesh,0)#CellFunction("size_t", mesh)
domains.set_all(0)
thin_zone.mark(domains,1)

dx = Measure("dx")[domains]

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0],0,tol)

def boundary_R(x, on_boundary):
	return on_boundary and near(x[0],1,tol)

boundaries = MeshFunction("size_t",mesh,0)
boundary_L = CompiledSubDomain('near(x[0], 0)')#'on_boundary && near(x[0], 0, tol)',tol=1E-14
#boundary_out = CompiledSubDomain('near(x[0], 0.98)')
#boundary_out.mark(boundaries,2)
boundary_L.mark(boundaries,1)
ds = Measure('ds')[boundaries]

class integration_section(SubDomain):
	def inside(self, x, on_boundary):
		return between(x[0], (0.96,1.0))

all_molecules = necessary_values['gas_num']+necessary_values['surf_num']+1
u = Function(V)
u_n = Function(V)
tempA = TestFunctions(V)
tempB = split(u)
tempC = split(u_n)

for kit in range(0,(all_molecules)+1):
	v_d['v_'+str(kit+1)] = tempA[kit]
	u_d['u_'+str(kit+1)] = tempB[kit]
	u_nd['u_n'+str(kit+1)] = tempC[kit]

#w = Constant(1)
W = VectorFunctionSpace(mesh, 'P', 1)
w = Function(W)

dt = Time/Time_steps

Length = np.sum(r_param)**2
#print(w.vector().get_local())
w.vector()[:] = 1
#print(w.vector().get_local())
#sys.exit()

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

bcs = []
#for k in range(0,reactants_num):
bcs.append(DirichletBC(V.sub(0),Constant(1),boundary_L))
bcs.append(DirichletBC(V.sub(1),Constant(0),boundary_L))
#bcs.append(DirichletBC(V.sub(2),Constant(0),boundary_L))
#bcs.append(DirichletBC(V.sub(3),Constant(0),boundary_L))
#bcs.append(DirichletBC(V.sub(4),Constant(0),boundary_L))
#bcs.append(DirichletBC(V.sub(5),Constant(0),boundary_L))
#bcs.append(DirichletBC(V.sub(6),Constant(1),boundary_L))

#	bcs.append(DirichletBC(V.sub(k),Constant(0),boundary_R))
bcs.append(DirichletBC(V.sub(all_molecules),Constant(1),boundary_L))
#bcs.append(DirichletBC(V.sub(all_molecules),Constant(0),boundary_R))

################################################################################################################################	### Define parameters / graph used in next sections
#fig2, ax2 = plt.subplots()
#ax2.set_xlabel('$t (s)$')
#ax2.set_ylabel('$Flux (molecules/s)$')
zef = np.linspace(dx_r, 1.-dx_r, grid_points+1)

to_flux = []
monitored_gas = necessary_values['molecules_in_gas_phase']

for k in range(0,monitored_gas):
	to_flux.append(1)
to_flux.append((1))

################################################################################################################################	### Solve the system
t = 0

J = derivative(F,u)
problem = NonlinearVariationalProblem(F,u,bcs,J)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8


#A = np.array((0,0.5,0.5,1))	
#B = np.array([[0,0,0,0],[0.5,0,0,0],[0,1/2,0,0],[0,0,1,0]])
#C = np.array((1/6,1/3,1/3,1/6))
#scheme = BDF1(F, u,t=dk,bcs=bcs)
#scheme = ButcherMultiStageScheme(F, u,time=dk,a=B,b=A,c=C,order=1,bcs=bcs)
#solver = RKSolver(scheme,1)

def solver_iteration(time_step):
	try:
		test1,test2 = solver.solve()
		return time_step,test1,test2 
	except RuntimeError:
		time_step = time_step*0.5
		dk.assign(time_step)
		if time_step < 1e-5:
			print("Time step too low")
			print(time.time() - start_time)
			sys.exit()
		time_step,test1,test2=solver_iteration(time_step)
		return time_step,test1,test2

cur_max = 0
tot_max = 0
start_time = time.time()
#g1 = Constant(0.001)
#L = g1*v_d['v_4']*ds(1)



osub = integration_section()
domains = MeshFunction("size_t", mesh,0)
domains.set_all(0)
osub.mark(domains, 1)
dP = Measure('vertex',domain = mesh, subdomain_data=domains)

sens_data['conVtime_1'].append(0)
sens_data['conVtime_2'].append(0)
start_time = time.time()
while t <= Time+0.01:	
	graph_data['timing'].append(t)																									### Store each time step
	#test_new = project(u,V)
	#sens_func = assemble(inner(u,u)*dP(1))
	#test_new_split = test_new.split(deepcopy=True)
	max_list=[]
	if reactions_test != ['INERT_ONLY']:																							### Store the outlet flux for each observed gas
		for k in range(0,monitored_gas):
			new_val = (( u.vector().get_local()[(all_molecules)+k+1]))
			max_list.append(new_val)
			graph_data['conVtime_'+str(k)].append((new_val))
		new_val = ((u.vector().get_local()[2*(all_molecules+1)-1]))
		graph_data['conVtime_'+str(all_molecules-1)].append((new_val))
		max_list.append(new_val)
	else:
		graph_data['conVtime_0'].append((( (u.vector().get_local()[1]))))

	cur_max = max(max_list)
	if cur_max > tot_max:
		tot_max = cur_max

	if t > 0:																														
		#dt,test1,test2 = solver_iteration(dt)
		#sens_func = assemble(test_new[k]*dP(1))
		
		#ens_func = assemble(test_new[k][0]*dP(1))
		#start_time = time.time()
		#dJdm1,dJdm2 = compute_gradient(sens_func,[control1,control2],forget=True)#compute_gradient(sens_func,[control1,control2])
		#print("")
		#print('sensitivity time')
		#print(time.time() - start_time)
		#print("dJdm: values for control 1 & 2")
		#print(dJdm1.values(),dJdm2.values())
		solve(F==0,u,bcs)
		newish = u.split(deepcopy=False)
		#ax2.plot(zef,newish[3].vector().get_local(), ls = '--', alpha=0.7)
	else:
		if reactions_test != ['INERT_ONLY']:
			#for z in range(int(grid_loc_1[0]-1),int(grid_loc_1[1])):
			#	u_n.vector()[z*(all_molecules+1)-2] = 100
				#u_n.vector()[z*(all_molecules+1)-(8-5)] = OB
				#u_n.vector()[z*(all_molecules+1)-(8-3)] = OA
			
			#print(type(sens_func))

			
			#for k in range(0,reactants_num):#necessary_values['reactants_number']
			#	u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)+k-(all_molecules)] = reactant_ratio[k]*Inert_pulse_size/(point_volume)
			#u_n.vector()[int((all_molecules+1)*(grid_points+1)-1)] = 10
			#newish = u_n.split(deepcopy=True)
			solve(F==0,u,bcs)
			#dt,test1,test2 = solver_iteration(dt)
			temp_again = u.split(deepcopy=False)
		else:
			u_n.vector()[int(grid_points)-1] = 10000
			solve(F== 0,u,bcs,solver_parameters=parms_solve)

	print("TIME STEP")
	print(t)
	u_n.assign(u)
	#if t > 1:
	#if cur_max < (1/500)*tot_max:
	#	dt = 0.05
	#	dk.assign(dt)
	#if test1 < 2:
	#	dt = dt+0.0001
	#	dk.assign(dt)
	t += dt
print('')
print('Total simulation time')
print(time.time() - start_time)
u_final = u.split(deepcopy=False)
#print(len(u_final))
#sys.exit()
#Ke0,Kd0,Ke1,Kd1,Ke2,Kd2,Ke3,Kd3

control1 = Control(Ke0)
control2 = Control(Kd0)
control3 = Control(Ke1)
control4 = Control(Kd1)
control5 = Control(Ke2)
control6 = Control(Kd2)
control7 = Control(Ke3)
control8 = Control(Kd3)
#sens_func = assemble(inner(u,u)*dP(1))
sens_func = assemble(inner(u_final[0],u_final[0])*dP(1))
#ens_func = assemble(test_new[k][0]*dP(1))
start_time = time.time()
print('sensitivity time')
dJdm1,dJdm2 = compute_gradient(sens_func,[control1,control2])#compute_gradient(sens_func,[control1,control2])
print(time.time() - start_time)
print("dJdm for each control")
print(dJdm1.values(),dJdm2.values())
sys.exit()
sens_func = assemble(inner(u_final[2],u_final[2])*dP(1))
#ens_func = assemble(test_new[k][0]*dP(1))
start_time = time.time()
print("")
print('sensitivity time')
dJdm1,dJdm2,dJdm3,dJdm4,dJdm5,dJdm6,dJdm7,dJdm8 = compute_gradient(sens_func,[control1,control2,control3,control4,control5,control6,control7,control8])#compute_gradient(sens_func,[control1,control2])
print(time.time() - start_time)
print('')
print("dJdm for each control")
print(dJdm1.values(),dJdm2.values(),dJdm3.values(),dJdm4.values(),dJdm5.values(),dJdm6.values(),dJdm7.values(),dJdm8.values())

#print(time.time() - start_time)


################################################################################################################################	### Generate the Legend & Color options for output graph
legend_label = []
if reactions_test[0] != 'INERT_ONLY':
	for k in range(0,necessary_values['molecules_in_gas_phase']):
		legend_label.append(necessary_values['reactants'][k])
legend_label.append("Inert")
colors = ['b','r','m','g','b','r']

################################################################################################################################	### Graph the output Curves
for k,j in enumerate(graph_data):
	if j != 'timing':
		ax2.plot(graph_data['timing'],graph_data[j],color=colors[k-1],label=legend_label[k-1], ls = '--', alpha=0.7)
	else:
		pass

ax2.legend(title="Gas Species")

plt.show()
sys.exit()