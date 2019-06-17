from fenics import *
from fenics_adjoint import *
from pyadjoint.enlisting import Enlist
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import imageio
import pandas as pd
import numpy as np
import math
import os


def read_input():
	
	"""
	Returns a dictionary of simulation parameters from 
	"""

	user_data = pd.read_csv('./input_file.csv',header=None)
	
	rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
	rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
	rows_4, cols_4 = np.where(user_data == 'Reaction_Information')

	reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = user_data.iloc[1+rows_2[0]:rows_3[0]-1,:]
	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]
	reaction_info = user_data.iloc[1+rows_4[0]:,:]

	reactor_kinetics_input = {}
	
	for k in range(0,len(reactor_info.index)):
		try:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]

	for k in range(0,len(data_storage.index)):
		try:
			reactor_kinetics_input[data_storage.iloc[k,0]] = float(data_storage.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[data_storage.iloc[k,0]] = data_storage.iloc[k,1]

	for k in range(0,len(feed_surf_info.index)):
		try:
			reactor_kinetics_input[feed_surf_info.iloc[k,0]] = float(feed_surf_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[feed_surf_info.iloc[k,0]] = feed_surf_info.iloc[k,1]

	reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
	reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']

	reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()

	kinetic_parameters = {}
	
	for j in range(0,len(reaction_info.index)):
		kinetic_parameters['kf'+str(j)] = float(reaction_info.iloc[j,1])
		if str(reaction_info.iloc[j,2]) != 'nan':
			kinetic_parameters['kb'+str(j)] = float(reaction_info.iloc[j,2])
		else:
			pass
	kin_in = kinetic_parameters.copy()

	return reactor_kinetics_input,kinetic_parameters,kin_in

def call_sens_analysis(u_value,control_list,domain):

	""" Perform the sensitivty analysis"""

	flux =  u_value*u_value*domain
	sens_func = assemble(flux)
	X = compute_gradient(sens_func,control_list)###################
	#X = hessian(sens_func,control_list)###################
	m = Enlist(control_list)
	grads = [i.get_derivative(options=None) for i in m]
	value_list = []
	Z = np.array(len(control_list))

	for nu in grads:
		value_list.append(nu.values().item(0))
	
	temp = np.array((value_list))
	
	return temp
	

def call_ad_rrm_analysis(u_value,control_list,domain):

	"""Similar to call_sens_analysis, calculates the derivatives for thin zone"""

	flux =  u_value*u_value*domain
	sens_func = assemble(flux)
	X = compute_hessian(sens_func,control_list)###################
	m = Enlist(control_list)
	grads = [i.get_derivative(options=None) for i in m]
	value_list = []
	Z = np.array(len(control_list))
	
	for nu in grads:
		value_list.append(nu.values().item(0))
	
	temp = np.array((value_list))
	
	return temp
	#sensitivity_output[k_step].append(temp)

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
		print('Time Step Failure')
		sys.exit()

	#except RuntimeError:
	#	time_step = time_step*0.5
	#	print(time_step)
	#	dk.assign(time_step)
	#	if time_step < 1e-6:
	#		print("Time step too low")
	#		print(time.time() - start_time)
	#		sys.exit()
	#	time_step=solver_iteration(time_step,method,solver,dk,1.5,1.1)
	#	return time_step

def flux_generation(reactor,gasses,reactants,pulse_size,Diff,voidage,dx,radius,dx2_r):
	
	"""Scale the output values of flux with the constant found here (messy now due to different trials)"""

	to_flux = []

	if reactor == 'tap':

		for k in range(0,gasses):
			to_flux.append( (Diff[k][0]*voidage[0] /(dx)) * (radius**2)*3.14159/(pulse_size) ) 
			#to_flux.append(2 *(dx*(radius**2)*3.14159) * (Diff[k][0] /(dx*voidage[0])))#(1/((1)*pulse_size)) *###??? changed from 1 to the new form
			#to_flux.append(2*Diff[k][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r))#(1/((1)*pulse_size)) *
			#to_flux.append(2*(1/((1+reactants)*pulse_size)) *Diff[k][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r))
		#to_flux.append( (Diff[gasses][0]*voidage[0] /(dx)) * (radius**2)*3.14159/(pulse_size) )
		#to_flux.append((2*Diff[gasses][0] * (dx*(radius**2)*3.14159)/(dx*voidage[0])))#*(1/((1)*pulse_size)) *
		#to_flux.append((2*(1/((1+reactants)*pulse_size)) *Diff[gasses][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r)))
	elif reactor_type == 't_pfr' or 't_pfr_diff':
		for k in range(0,gasses):
			to_flux.append(1)
		
		to_flux.append(1)
	
	return to_flux

def define_boundary_conditions(reactor,elem_list,nec_values,V_sec,reacs_num,all_mol,reac_ratio,L_bound,R_bound,number_of_inerts):

	"""Define the appropriate boundary conditions for all monitored species"""

	if reactor == 'tap':
		bcs = []

		if elem_list != ['INERT_ONLY']:

			for k in range(0,nec_values):
				bcs.append(DirichletBC(V_sec.sub(k),Constant(0),R_bound))
			
			for k in range(0,int(number_of_inerts)):
				bcs.append(DirichletBC(V_sec.sub(all_mol-(1+k)),Constant(0),R_bound))
			
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
	
	"""For all monitored parameters (gasses/surface species), establish necessary test and trial functions, as well as dictionaries for value storing"""

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
	sens_data['conVtime_'+str(species_count)] = []
	
	for kj in range(nec['molecules_in_gas_phase'],len(nec['reactants'])):
		surf_data['conVtime_'+str(kj)] = []
	
	tempA = TestFunctions(V_nu)
	tempB = split(u_nu)
	tempC = split(un_nu)
	
	for kit in range(0,int(moles)):
		v_d['v_'+str(kit+1)] = tempA[kit]
		u_d['u_'+str(kit+1)] = tempB[kit]
		u_nd['u_n'+str(kit+1)] = tempC[kit]
	
	return graph_data,v_d,u_d,u_nd,sens_data,surf_data,cat_data


def establish_grid_system(in_1,cat,in_2,mesh_size):

	"""Generate the FEniCS Mesh"""

	r_param = np.array((in_1,cat,in_2))
	
	def last_grid_point(x):
	    return x*mesh_size

	last_func = np.vectorize(last_grid_point)
	zone_sizes = last_func(r_param)
	grid_loc_1 = np.round(np.cumsum(zone_sizes))
	dx_r = np.sum(r_param)/mesh_size
	dx2_r = dx_r*dx_r

	frac_length = r_param[1]/(np.sum(r_param))
	cat_location = (r_param[1]*0.5+r_param[0])/(np.sum(r_param))

	return r_param,dx_r,dx2_r,frac_length,cat_location


def establish_output_graph(reactor,gas_phase,reacs,inerts):

	"""Generate the outlet flux graph (initialize, anyway)"""

	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('$t\ (s)$')
	if reactor == 'tap':
		ax2.set_ylabel('$Outlet\ Flux\ (1/s)$')
	elif reactor == 't_pfr' or 't_pfr_diff':
		ax2.set_ylabel('$Outlet Concentration (molecules/cm3)$')
	
	legend_label = []
	header = "time"
	for k in range(0,gas_phase):
		legend_label.append(reacs[k])
		header = header+","+reacs[k]
	for j in range(0,inerts):
		legend_label.append("Inert-"+str(1+j))
	header = header+",Inert"
	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']

	return fig2,ax2,legend_label,header,colors


def exp_data_fitting_all(species_list,time,folder):

	"""Define the objective function in terms of every simulated point"""

	print("This method is not recommended and will lead to an endless calculation. Will continue after six seconds.")
	time.sleep(6)

	def interp(y2,y1,x2,x1,b,xn):
		print('used')
		return ((y2 - y1)/(x2 - x1))*xn + b

	user_data = {}
	species_list = species_list[:len(species_list)]
	print(species_list)
	time_steps = int(time)
	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+species_list[k]+'.csv',header=None)
	curve_fitting = {}
	exp_data = user_data
	
	for k_new in species_list:
		time_step = []
		times = []
		values = []
		for j in range(0,time_steps):
			time_step.append(j)
			times.append(round(user_data[k_new].iloc[j,0],6))
			values.append(user_data[k_new].iloc[j,1])

		data = {}
		data['time_step'] = time_step
		data['times'] = times
		data['values'] = values

		curve_fitting[k_new] = data

	return curve_fitting


def exp_data_fitting(species_list,sim_steps,folder,time,points):

	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = time/sim_steps
	
	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1

	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'

	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	
	curve_fitting = {}
	exp_data = user_data

	for k_new in species_list:

		def find_experimental_point(n,exp_step):
			""" Find an appropriate intensity point for the fitting process """
			approx_exp_n = n*(syn_time_step)/exp_step
			#print(n*(syn_time_step))
			
			if approx_exp_n != n:
				high = math.ceil(approx_exp_n)
				low = int(approx_exp_n)
				#print(interp(user_data[k_new][1][high],user_data[k_new][1][low],user_data[k_new][0][high],user_data[k_new][0][low],n*(syn_time_step)))
				return interp(user_data[k_new][1][high],user_data[k_new][1][low],user_data[k_new][0][high],user_data[k_new][0][low],n*(syn_time_step))

			else:
				return user_data[k_new][1][n]

		def exp_point_to_syn_point(n_exp,exp_step):
			"""Align an experimental data point with the associated (or nearest synthetic point)"""
			approx_syn_n = n_exp*exp_step/(syn_time_step)
			
			if int(approx_syn_n) > 0:
				return find_experimental_point(int(approx_syn_n),exp_step)
			else:
				return find_experimental_point(math.ceil(approx_syn_n),exp_step) 

		time_step = []
		times = []
		values = []
		exp_time_step = user_data[k_new][0][1]
		near_start = round(user_data[k_new].iloc[30,0],6)/(time/sim_steps)
	
		peak_loc = user_data[k_new].iloc[user_data[k_new][1].idxmax()]
		
		near_peak = peak_loc[0]/(time/sim_steps)
		peak2 = user_data[k_new].loc[user_data[k_new][0] == peak_loc[0]].index

		test3 = int(round((peak2[0]+1)/2,0))
		mid_loc = user_data[k_new].iloc[test3,:]
		near_mid = mid_loc[0]/(time/sim_steps)

		time_step.append(int(near_mid))
		times.append(int(near_mid)*(syn_time_step))
		values.append(find_experimental_point(int(near_mid),exp_time_step))

		if points > 1:
			time_step.append(int(near_peak))
			times.append(int(near_peak)*(syn_time_step))
			values.append(find_experimental_point(int(near_peak),exp_time_step))

		if points > 2:
			value_test = 0.85*peak_loc[1]
			sort_3 = user_data[k_new].iloc[(user_data[k_new][1]-value_test).abs().argsort()[:]]

			for k in range(0,sort_3.index[0]):
				if sort_3.iloc[k,0] > peak_loc[0]:
					thr_point = sort_3.iloc[k]
					break
				else:
					pass

			near_3 = thr_point[0]/(time/sim_steps)

			time_step.append(int(near_3))
			times.append(int(near_3)*(syn_time_step))
			values.append(find_experimental_point(int(near_peak),exp_time_step))
		
		if points > 3:
			value_test = 0.75*peak_loc[1]
			sort_4 = user_data[k_new].iloc[(user_data[k_new][1]-value_test).abs().argsort()[:]]

			for k in range(0,sort_4.index[0]):
				if sort_4.iloc[k,0] > peak_loc[0]:
					four_point = sort_4.iloc[k]
					break
				else:
					pass

			near_4 = four_point[0]/(time/sim_steps)

			time_step.append(int(near_4))
			times.append(int(near_4)*(syn_time_step))
			values.append(find_experimental_point(int(near_peak),exp_time_step))
			

		data = {}
		data['time_step'] = time_step
		data['times'] = times
		data['values'] = values

		curve_fitting[k_new] = data

	return curve_fitting


def generate_gif(molecules,exp_loc,fit_loc,all_steps,constants,reactions,time_data):
	
	"""
	Return a gif showing changes made during the optimization process
	"""

	def add_subplot_axes(ax,rect,axisbg='w'):
		
		"""
		Generates the subplot to help visualize some other quantitiy
		"""

		fig = plt.gcf()
		box = ax.get_position()
		width = box.width
		height = box.height
		inax_position  = ax.transAxes.transform(rect[0:2])
		transFigure = fig.transFigure.inverted()
		infig_position = transFigure.transform(inax_position)    
		
		x = infig_position[0]
		y = infig_position[1]
		width *= rect[2]
		height *= rect[3]
		subax = fig.add_axes([x,y,width,height])
		
		x_labelsize = subax.get_xticklabels()[0].get_size()
		y_labelsize = subax.get_yticklabels()[0].get_size()
		x_labelsize *= rect[2]**0.5
		y_labelsize *= rect[3]**0.5
		
		subax.xaxis.set_tick_params(labelsize=x_labelsize)
		subax.yaxis.set_tick_params(labelsize=y_labelsize)
		
		return subax
	
	x_data = list(range(0, all_steps))
	y_data = [0]
	
	for k in range(0,all_steps-1):
		y_data.append( (time_data[k+1] - time_data[k]) / 60 )

	def tap_plot(step):
		fig, ax = plt.subplots(figsize=(10,5))
	
		ax.grid()
		
		exp_data = {}
		sim_data = {}
		for k_names in molecules:

			exp_data[k_names] = pd.read_csv(exp_loc+'/'+k_names+'.csv',header=None)
			sim_data[k_names] = pd.read_csv(fit_loc+'/iter_'+str(step)+'_folder/flux_data/'+k_names+'.csv',header=None)
		
		for k_names in molecules:
			ax.plot(exp_data[k_names][0], exp_data[k_names][1])
		
		for k_names in molecules:
			ax.plot(sim_data[k_names][0], sim_data[k_names][1],ls='--')
	
		ax.set_xlabel('Time (s)', fontsize=16)
		ax.set_ylabel('Flux (nmol/s)', fontsize=16)
		textstr = 'Rate Constants: '
		
		for k_len in range(0,len(constants[0])):
			textstr = '\n'.join((textstr,'k'+str(1+k_len)+': '+'{:0.3e}'.format(constants[step][k_len])))

		props = dict(facecolor='white')
		ax.text(0.8, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props,ha='center')
		textstr = 'Elementary Reactions:'
		
		for k_reacs in range(0,len(reactions)):
			textstr = '\n'.join((textstr,reactions[k_reacs]))

		ax.text(0.5, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props,ha='center')
		ax.text(0.15, 0.95,'Iteration: '+str(step),transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props)
		peak_peak = 0

		for k_names in molecules[:-1]:
			peak_loc = exp_data[k_names].iloc[exp_data[k_names][1].idxmax()]
			plt.plot(peak_loc[0], peak_loc[1], 'ro')
			if peak_loc[1] > peak_peak:
				peak_peak = peak_loc[1]

			peak2 = exp_data[k_names].loc[exp_data[k_names][0] == exp_data[0]].index
			test3 = int(round((peak2[0]+1)/2,0))
			mid_loc = exp_data[k_names].iloc[test3,:]
			plt.plot(mid_loc[0], mid_loc[1], 'ro')

			#peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
			#peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
			#plt.plot(peak_loc[0], peak_loc[1], 'ro')

		ax.set_ylim(0,peak_peak*1.1)
		
		subpos = [0.4,0.13,0.5,0.4]
		subax1 = add_subplot_axes(ax,subpos)
		subax1.plot(x_data[:step],y_data[:step])
		subax1.set_title('Time per Iteration')
		subax1.set_ylabel('Time (minutes)')
		subax1.set_xlabel('Iteration #')
		subax1.set_xlim(0,all_steps)
		subax1.set_ylim(0,max(y_data)*1.1)
		
		fig.canvas.draw()
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
		return image
	
	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	imageio.mimsave(fit_loc+'/output.gif', [tap_plot(i) for i in range(all_steps)], fps=4)


"""Functions used to keep output organized"""

def generate_folder(path_name):
	try:  
		os.mkdir(path_name)
	except OSError:  
		###print ("Creation of the directory %s failed" % path_name)
		pass
	else:  
		###print ("Successfully created the directory %s " % path_name)
		pass
	
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


def progressBar(value, endvalue, bar_length=20):
	
	""" Generate the progress bar"""

	percent = float(value) / endvalue
	arrow = '-' * int(round(percent * bar_length)-1) + '>'
	spaces = ' ' * (bar_length - len(arrow))
	sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, round(percent * 100,3)))
	sys.stdout.flush()

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


def error_output(elem_reacs):
	
	"""Return this error in the event that the user doesn't define rate constants correct"""

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
	print("'kf1' = forward rate constant for reaction # 1")
	print("'kb1' = reverse rate constant for reaction # 1")
	print("       ")
	print("Rate constants for the following equations")
	print("must be included")
	print("       ")
	for j_nu,k_nu in enumerate(elem_reacs):
		print("Reaction "+str(j_nu+1)+"     "+k_nu)
	print("       ")
	print("Could also have incorrect diffusion coefficient array defined")
	sys.exit()