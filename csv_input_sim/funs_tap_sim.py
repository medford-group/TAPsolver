from fenics import *
from fenics_adjoint import *
from pyadjoint.enlisting import Enlist
import matplotlib.pyplot as plt
import numpy as np
import os

def call_sens_analysis(u_value,control_list,domain):
	flux =  u_value*u_value*domain
	sens_func = assemble(flux)
	X = compute_gradient(sens_func,control_list)###################
	m = Enlist(control_list)
	grads = [i.get_derivative(options=None) for i in m]
	value_list = []
	Z = np.array(len(control_list))
	for nu in grads:
		value_list.append(nu.values().item(0))
	temp = np.array((value_list))
	return temp
	#sensitivity_output[k_step].append(temp)


def call_solver(dk,u_temp,u,u_new,solver_def,keep_sol = True):
	u_temp.assign(u)
	solver_def.solve()
	u_new.assign(u)
	if keep_sol == False:
		u.assign(u_temp)
	return u_new


def norm_comp(u1,u2,u3,d_t,i_t):
	ref_norm = u2.vector().norm("l2")
	norm1 = u1.vector().norm('l2')/ref_norm
	norm2 = u3.vector().norm('l2')/ref_norm

	if norm1 > 0.02:
		time_step = time_step/d_t
		dk.assign(time_step)
		u.assign(u1)
	elif norm2 < 0.01:
		time_step = time_step*i_t
		dk.assign(time_step)
		u.assign(u3)
	return time_step

def solver_iteration(time_step,method,solver,dk,dec_tim,inc_tim):
	try:
		if method == 'simple_adaptive':
		
			u_temp.assign(u)
			uout_1 = call_solver(dk.assign(time_step/dec_tim),u_temp,u,u_1,solver)
			uout_3 = call_solver(dk.assign(time_step*inc_tim),u_temp,u,u_3,solver)
			uout_2 = call_solver(dk.assign(time_step),u_temp,u,u_2,solver,keep_sol=True)
			time_step = norm_comp(uout_1,uout_2,uout_3,dec_tim,inc_tim)
			return time_step
		
		elif method == 'None':
			solver.solve()
			return time_step
	except RuntimeError:
		time_step = time_step*0.5
		print(time_step)
		dk.assign(time_step)
		if time_step < 1e-10:
			print("Time step too low")
			print(time.time() - start_time)
			sys.exit()
		time_step=solver_iteration(time_step,method,solver,dk,1.5,1.1)
		return time_step


def progressBar(value, endvalue, bar_length=20):
	percent = float(value) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, round(percent * 100,3)))
	sys.stdout.flush()


def error_output(elem_reacs):
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
	for j_nu,k_nu in enumerate(elem_reacs):
		print("Reaction "+str(j_nu+1)+"     "+k_nu)
	print("       ")
	print("Could also have incorrect diffusion coefficient array defined")
	sys.exit()

	
def store_sens_analysis(yes_no,output_file,gasses,legend_ref):
	if yes_no == True:
		try:
			os.makedirs('./'+output_file+'_folder/sensitivity_'+output_file)
		except OSError:
			pass
		for k_sens_folder in range(gasses):	
			try:
				os.makedirs('./'+output_file+'_folder/sensitivity_'+output_file+'/'+legend_ref[k_sens_folder])
			except OSError:
				pass


def store_data_func(yes_no,output_file):
	if yes_no == True:
		try:
			os.makedirs('./'+output_file+'_folder')
		except OSError:
			pass


def flux_generation(reactor,gasses,reactants,pulse_size,Diff,voidage,dx,radius,dx2_r):
	to_flux = []
	if reactor == 'tap':
		for k in range(0,gasses):
			#print(dx)
			#sys.exit()
			#to_flux.append(1)
			to_flux.append(2 (dx*(radius**2)*3.14159) * (Diff[k][0] /(voidage[0]*dx2_r)))#(1/((1)*pulse_size)) *
			#to_flux.append(2*Diff[k][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r))#(1/((1)*pulse_size)) *
			#to_flux.append(2*(1/((1+reactants)*pulse_size)) *Diff[k][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r))
		#to_flux.append(1)
		to_flux.append((2*Diff[gasses][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r)))#*(1/((1)*pulse_size)) *
		#to_flux.append((2*(1/((1+reactants)*pulse_size)) *Diff[gasses][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r)))
	elif reactor_type == 't_pfr' or 't_pfr_diff':
		for k in range(0,gasses):
			to_flux.append(1)
		to_flux.append(1)
	#print(to_flux)
	#sys.exit()
	return to_flux


def define_boundary_conditions(reactor,elem_list,nec_values,V_sec,reacs_num,all_mol,reac_ratio,L_bound,R_bound):
	if reactor == 'tap':
		bcs = []
		if elem_list != ['INERT_ONLY']:
			for k in range(0,nec_values):
				bcs.append(DirichletBC(V_sec.sub(k),Constant(0),R_bound))
			bcs.append(DirichletBC(V_sec.sub(all_mol),Constant(0),R_bound))
		else:	
			bcs.append(DirichletBC(V_sec,Constant(0),R_bound))
	
	elif reactor == 't_pfr' or 't_pfr_diff':
		bcs = []
		for k in range(0,int(reacs_num)):
			bcs.append(DirichletBC(V_sec.sub(k),Constant(reac_ratio[k]),L_bound))
		for kin in range(reacs_num,all_mol):
			bcs.append(DirichletBC(V_sec.sub(kin),Constant(0),L_bound))
		bcs.append(DirichletBC(V_sec.sub(all_mol),Constant(0),L_bound))
	return bcs

def initialize_variable_dictionaries(nec,moles,V_nu,u_nu,un_nu):
	graph_data = {}
	v_d = {}
	u_d = {}
	u_nd = {}
	surf_data = {}
	sens_data = {}	
	cat_data = {}

	species_count = nec['gas_num']+nec['surf_num']

	for k in range(0,nec['molecules_in_gas_phase']):
		graph_data['conVtime_'+str(k)] = []
		cat_data['conVtime_'+str(k)] = []
		sens_data['conVtime_'+str(k)] = []
	graph_data['conVtime_'+str(species_count)] = []
	#cat_data['conVtime_'+str(nec['molecules_in_gas_phase'])] = []
	sens_data['conVtime_'+str(species_count)] = []

	for kj in range(nec['molecules_in_gas_phase'],len(nec['reactants'])-1):
		surf_data['conVtime_'+str(kj)] = []

	tempA = TestFunctions(V_nu)
	tempB = split(u_nu)
	tempC = split(un_nu)

	for kit in range(0,(moles)+1):
		v_d['v_'+str(kit+1)] = tempA[kit]
		u_d['u_'+str(kit+1)] = tempB[kit]
		u_nd['u_n'+str(kit+1)] = tempC[kit]

	return graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data

def establish_grid_system(in_1,cat,in_2,mesh_size):

	r_param = np.array((in_1,cat,in_2))
	#grid_points = round(np.sum(r_param)*mesh_space)	###Number of grid points for gasses

	def last_grid_point(x):						###Last grid points for gaseous species
	    return x*mesh_size

	last_func = np.vectorize(last_grid_point)
	zone_sizes = last_func(r_param)
	grid_loc_1 = np.round(np.cumsum(zone_sizes))
	dx_r = np.sum(r_param)/mesh_size
	dx2_r = dx_r*dx_r 							###(cm2) for calculation of alpha

	### Side work on developing general mesh implementation

	frac_length = r_param[1]/(np.sum(r_param))
	cat_location = (r_param[1]*0.5+r_param[0])/(np.sum(r_param))

	return r_param,dx_r,dx2_r,frac_length,cat_location

def establish_output_graph(reactor,gas_phase,reacs):
	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('$t (s)$')
	if reactor == 'tap':
		ax2.set_ylabel('$Dimensionless Flux (1/s)$')
	elif reactor == 't_pfr' or 't_pfr_diff':
		ax2.set_ylabel('$Outlet Concentration (molecules/cm3)$')
	#zef = np.linspace(dx_r, 1.-dx_r, grid_points+1)
	
	legend_label = []
	header = "time"
	for k in range(0,gas_phase):
		legend_label.append(reacs[k])
		header = header+","+reacs[k]
	legend_label.append("Inert")
	header = header+",Inert"
	colors = ['b','r','g','m','k','y','c']

	return fig2,ax2,legend_label,header,colors