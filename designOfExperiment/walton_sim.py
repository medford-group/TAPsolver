#Import all packages
from fenics import *
from fenics_adjoint import *
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
import warnings

#Control FEniCS output / information

def progressBar(value, endvalue, bar_length=20):
	percent = float(value) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, round(percent * 100,3)))
	sys.stdout.flush()

def run_simulation(equation_def,monitored_gases,sim_dict):

	if (sim_dict['store_data'].lower() == 'true' or sim_dict['store_graph'].lower() == 'true'):
		path = './'+sim_dict['folder_name']+'_folder/'
		try:  
			os.mkdir(path)
		except OSError:  
			print ("Creation of the directory %s failed" % path)
			pass
		else:  
			print ("Successfully created the directory %s " % path)
		#if reac_input['save_figure'].lower() == 'true':

		user_data = pd.read_csv('./input_file.csv',header=None)
		user_data.to_csv(path+'input_file.csv',header=None)

		if sim_dict['store_data'].lower() == 'true':
			path = './'+sim_dict['folder_name']+'_folder/data/'
			try:  
				os.mkdir(path)
			except OSError:  
				print ("Creation of the directory %s failed" % path)
				pass
			else:
				print ("Successfully created the directory %s " % path)

		if sim_dict['store_graph'].lower() == 'true':
			path = './'+sim_dict['folder_name']+'_folder/graphs/'
			try:  
				os.mkdir(path)
			except OSError:  
				print ("Creation of the directory %s failed" % path)
				pass
			else:  
				print ("Successfully created the directory %s " % path)

	K_list = equation_def.iloc[1,1:].tolist()
	M_list = equation_def.iloc[2,1:].tolist()
	ka_list = equation_def.iloc[3,1:].tolist()
	V_list = equation_def.iloc[4,1:].tolist()
	a_list = equation_def.iloc[5,1:].tolist()
	b_list = equation_def.iloc[6,1:].tolist()
	bc_list = equation_def.iloc[7,1:].tolist()

	colors = ['b','r','g','m','k','y','c']
	ls_list = ['-','--','-.',':']

	#Define input constants
	mesh_size = int(sim_dict['mesh_size'])
	length = sim_dict['length']
	K_num = 1
	M_num = 1
	ka_num = 1
	V_num = 1
	a_num = 2
	b_num = 3
	
	K_con = Constant(K_num)
	ka_con = Constant(ka_num)
	V_con = Constant(V_num)
	a_con = Constant(a_num)
	b_con = Constant(b_num)
	
	
	parameters["std_out_all_processes"] = False																							
	cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
	#set_log_active(False)
	warnings.filterwarnings("ignore", category=DeprecationWarning)
	tol = 1E-14
	
	graph_data = {}
	v_d = {}
	u_d = {}
	u_nd = {}

	graph_data = {}
	for k in range(0,len(monitored_gases)*2):
		graph_data['conVtime_'+str(k)] = []

	mesh = UnitIntervalMesh(mesh_size)
	P1 = FiniteElement('P',mesh.ufl_cell(),1)

	element = '['

	for k in range(0,len(monitored_gases)*2):
		if (k+1) < len(monitored_gases)*2:
			element = element + 'P1,'
		else:
			element = element + 'P1]'

	test_new = eval(element)
	element = MixedElement(test_new)
	V = FunctionSpace(mesh,element)
	v_du = FunctionSpace(mesh,P1)
	
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
		
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)

	bcs=[]
	for k in range(0,len(monitored_gases)):
		bcs.append(DirichletBC(V.sub(k),Constant(float(bc_list[k])),boundary_L))

	all_molecules = len(monitored_gases)*2
	u = Function(V)
	u_n = Function(V)
	tempA = TestFunctions(V)
	tempB = split(u)
	tempC = split(u_n)

	for kit in range(0,len(monitored_gases)*2):
		v_d['v_'+str(kit)] = tempA[kit]
		u_d['u_'+str(kit)] = tempB[kit]
		u_nd['u_n'+str(kit)] = tempC[kit]

	w_list = []

	for k in range(0,len(monitored_gases)):
		W = VectorFunctionSpace(mesh, 'P', 1)
		w = Function(W)
		w.vector()[:] = float(M_list[k])
		w_list.append(w)

	sim_time = sim_dict['sim_time']
	time_steps = sim_dict['time_steps']
	
	dt = sim_time/time_steps
	dk = Constant(dt)

	F_str = ''

	for k in range(0,len(monitored_gases)):
		F_str += ' (u_d["u_'+str(k)+'"] - u_nd["u_n'+str(k)+'"])*v_d["v_'+str(k)+'"]*dx +'
		F_str += '(dk/Constant(length))*dot(w_list['+str(k)+'], grad(u_d["u_'+str(k)+'"]))*v_d["v_'+str(k)+'"]*dx' 
		F_str += ' + ((u_d["u_'+str(k+len(monitored_gases))+'"] - u_nd["u_n'+str(k+len(monitored_gases))+'"]))*v_d["v_'+str(k+len(monitored_gases))+'"]*dx' 
		F_str += '+ dk*Constant(ka_list['+str(k)+'])*((Constant(V_list['+str(k)+'])*Constant(b_list['+str(k)+'])*u_d["u_'+str(k)+'"])/(1 + Constant(a_list['+str(k)+'])*u_d["u_'+str(k)+'"]) - u_d["u_'+str(k+len(monitored_gases))+'"])*v_d["v_'+str(k)+'"]*dx'
		F_str += '- dk*Constant(ka_list['+str(k)+'])*((Constant(V_list['+str(k)+'])*Constant(b_list['+str(k)+'])*u_d["u_'+str(k)+'"])/(1 + Constant(a_list['+str(k)+'])*u_d["u_'+str(k)+'"]) - u_d["u_'+str(k+len(monitored_gases))+'"])*v_d["v_'+str(k+len(monitored_gases))+'"]*dx'
		
		if k != len(monitored_gases)-1:
			F_str += '+'
	print(ka_list)
	#ka_list[0] = Control(ka_list[0])
	
	F = eval(F_str)

	J = derivative(F,u)
	problem = NonlinearVariationalProblem(F,u,bcs,J)
	solver = NonlinearVariationalSolver(problem)
	solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
	#x1=0
	#x2=sim_time+1
	#list2 = [x/sim_time for x in range(x1, x2)]
	#list2 = list(reversed(list2))
	
	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('Time (s)')
	ax2.set_ylabel('Outlet Concentration')
	time_list = []
	t = 0
	#print()
	for n in range(int(time_steps)):
		progressBar(t, sim_time)
		solver.solve()
		time_list.append(t)
		for k in range(0,len(monitored_gases)*2):
			new_val = (( u_n.vector().get_local()[k]))
			#print(new_val)
			graph_data['conVtime_'+str(k)].append((new_val))
		t += dt
		#time.sleep(.1)
		#solve(F == 0, u)
		u_n.assign(u)

	for k in range(len(monitored_gases)):
		ax2.plot(time_list,graph_data['conVtime_'+str(k)],color=colors[k],label=monitored_gases[k],ls=ls_list[k])

	if sim_dict['store_data'].lower() == 'true':
		np.savetxt('./'+sim_dict['folder_name']+'_folder/data/time.csv', time_list, delimiter=",")
		for j_species in range(0,len(monitored_gases)):
			#dictionary_of_numpy_data[legend_label[j_species]] = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
			np.savetxt('./'+sim_dict['folder_name']+'_folder/data/'+monitored_gases[j_species]+'.csv', graph_data['conVtime_'+str(j_species)], delimiter=",")

	if sim_dict['store_graph'].lower() == 'true':
		plt.savefig('./'+sim_dict['folder_name']+'_folder/graphs/concentration_profile.png')
	
	ax2.legend(title="Species")
	
	if sim_dict['display_graph'].lower() == 'true':
		plt.show()

	sys.exit()

user_data = pd.read_csv('./input_file.csv',header=None)

rows_1, cols_1 = np.where(user_data == 'equation_def')
rows_2, cols_2 = np.where(user_data == 'sim_options')

equation_def = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
sim_options = user_data.iloc[1+rows_2[0]:,:]
monitored_gases = equation_def.iloc[0,1:].tolist()

sim_dict = {}

for k in range(0,len(sim_options.index)):
	try:
		sim_dict[sim_options.iloc[k,0]] = float(sim_options.iloc[k,1]) 
	except ValueError:
		sim_dict[sim_options.iloc[k,0]] = sim_options.iloc[k,1]

run_simulation(equation_def,monitored_gases,sim_dict)