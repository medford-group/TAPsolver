from fenics import *
from fenics_adjoint import *
from pyadjoint.enlisting import Enlist
#import matplotlib
#matplotlib.use('agg')
import time
import matplotlib.pyplot as plt
import imageio
import pandas as pd
import numpy as np
import math
import os
import fenics_adjoint

def readInput(input_file,inputForm = 'old'):
	
	"""
	Convert the input file into dictionaries for TAPsolver to use
	"""
	if inputForm == 'old':
	
		user_data = pd.read_csv(input_file,header=None)
		
		rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
		rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	#	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
		rows_4, cols_4 = np.where(user_data == 'Reaction_Information')
	
		reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
		feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
	#	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]
		reaction_info = user_data.iloc[1+rows_4[0]:,:]
	
		reactor_kinetics_input = {}
		
		for k in range(0,len(reactor_info.index)):
			try:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
	
	#	for k in range(0,len(data_storage.index)):
	#		try:
	#			reactor_kinetics_input[data_storage.iloc[k,0]] = float(data_storage.iloc[k,1]) 
	#		except ValueError:
	#			reactor_kinetics_input[data_storage.iloc[k,0]] = data_storage.iloc[k,1]
	
		for k in range(0,len(feed_surf_info.index)):
			try:
				reactor_kinetics_input[feed_surf_info.iloc[k,0]] = float(feed_surf_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[feed_surf_info.iloc[k,0]] = feed_surf_info.iloc[k,1]
	
	
		### Change this set of code!
		#reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
		#reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
		#reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	
		reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']*(reactor_kinetics_input['Catalyst Location'] - 0.5*(reactor_kinetics_input['Catalyst Fraction']))#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
		reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length'] - (reactor_kinetics_input['len_cat'] + reactor_kinetics_input['len_inert_1'])#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	
		reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()

	else:
		
		user_data = pd.read_csv(input_file,header=None)
		
		rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
		rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	#	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
		rows_4, cols_4 = np.where(user_data == 'Reaction_Information')
		
		reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 

		feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
	#	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]
		reaction_info = user_data.iloc[1+rows_4[0]:,:]
		
		reactor_kinetics_input = {}
		
		number_of_gasses = 0
		number_of_surface = 0
		for z in feed_surf_info.iloc[0,1:]:
			if type(z) == (str):
				number_of_gasses += 1 
		for z in feed_surf_info.iloc[5,1:]:
			if type(z) == (str):
				number_of_surface += 1 

		gas_species = feed_surf_info.iloc[0,1:number_of_gasses+1]
		surface_species = feed_surf_info.iloc[6,1:number_of_surface+1]
		
		reactor_kinetics_input['Number of Inerts'] = 0
		reactor_kinetics_input['Pulse Size'] = ''
		reactor_kinetics_input['Pulse Time'] = ''
		reactor_kinetics_input['Mass List'] = ''
		reactor_kinetics_input['Number of Reactants'] = 0

		for jnum,j in enumerate(gas_species):
			
			if j.find('Inert') == 0:
				reactor_kinetics_input['Number of Inerts'] += 1
			if jnum == len(gas_species)-1:

				reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])
				reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])
				reactor_kinetics_input['Mass List'] = reactor_kinetics_input['Mass List']+str(feed_surf_info.iloc[3,1+jnum])
				reactor_kinetics_input['Number of Reactants'] += 1
			else:
				if (float(feed_surf_info.iloc[1,1+jnum]) != 0.0) and (j.find('Inert') == False):
					reactor_kinetics_input['Number of Reactants'] += 1
				reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])+','
				reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])+','
				reactor_kinetics_input['Mass List'] = reactor_kinetics_input['Mass List']+str(feed_surf_info.iloc[3,1+jnum])+','
		
		reactor_kinetics_input['Initial Surface Composition'] = ''

		for jnum,j in enumerate(surface_species):

			if jnum == len(surface_species)-1:
				reactor_kinetics_input['Initial Surface Composition'] = reactor_kinetics_input['Initial Surface Composition']+str(feed_surf_info.iloc[6,1+jnum])
			else:
				reactor_kinetics_input['Initial Surface Composition'] = reactor_kinetics_input['Initial Surface Composition']+str(feed_surf_info.iloc[6,1+jnum])+','

		for k in range(0,len(reactor_info.index)):
			try:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
		
		reactor_kinetics_input['Reactor Length'] =  float(reactor_info.iloc[0,1]) + float(reactor_info.iloc[0,2]) + float(reactor_info.iloc[0,3])		
		
		reactor_kinetics_input['Catalyst Fraction'] = float(reactor_info.iloc[0,2])/reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['Catalyst Location'] = (float(reactor_info.iloc[0,1])+(float(reactor_info.iloc[0,2])/2))/reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['Void Fraction Inert'] = float(reactor_info.iloc[1,1])
		reactor_kinetics_input['Void Fraction Catalyst'] = float(reactor_info.iloc[1,2])

		reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']*(reactor_kinetics_input['Catalyst Location'] - 0.5*(reactor_kinetics_input['Catalyst Fraction']))#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
		reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length'] - (reactor_kinetics_input['len_cat'] + reactor_kinetics_input['len_inert_1'])#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	
		reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()


	
	kinetic_parameters = {}
	Ao = {}
	Ea = {}
	Ga = {}
	dG = {}

	fittingParametersList = []

	gForward = []
	arrForward = []
	arrBackward = []
	
	for j in range(0,len(reaction_info.index)):

		if reaction_info.iloc[j,1].find("#") > 0:
			Anew, Eanew = reaction_info.iloc[j,1].split("#")
			if Anew.find("!") < 0:
				Ga['Ga'+str(j)] = float(Anew)
				kinetic_parameters['Ga'+str(j)] = float(Anew)
				fittingParametersList.append('Ga'+str(j))
			
			else:
				Ga['Ga'+str(j)] = float(Anew[:-1])

			if Eanew.find("!") < 0:
				dG['dG'+str(j)] = float(Eanew)
				fittingParametersList.append('dG'+str(j))
			
			else:
				dG['dG'+str(j)] = float(Eanew[:-1])

			gForward.append(j)

		elif reaction_info.iloc[j,1].find("$") > 0:
			Anew, Eanew = reaction_info.iloc[j,1].split("$")
			if Anew.find("!") < 0:
				Ao['Aof'+str(j)] = float(Anew)
				fittingParametersList.append('Aof'+str(j))
			
			else:
				Ao['Aof'+str(j)] = float(Anew[:-1])

			if Eanew.find("!") < 0:
				Ea['Eaf'+str(j)] = float(Eanew)
				fittingParametersList.append('Eaf'+str(j))
			
			else:
				Ea['Eaf'+str(j)] = float(Eanew[:-1])

			arrForward.append(j)

		else:
			if reaction_info.iloc[j,1].find("!") < 0:
				kinetic_parameters['kf'+str(j)] = float(reaction_info.iloc[j,1])
				fittingParametersList.append('kf'+str(j))
			
			else:
				new_value = float(reaction_info.iloc[j,1][:-1])
				kinetic_parameters['kf'+str(j)] = new_value#float(reaction_info.iloc[j,1])

		if str(reaction_info.iloc[j,2]) != 'nan':
			if str(reaction_info.iloc[j,2]).find("$") > 0:
				Anew, Eanew = str(reaction_info.iloc[j,2]).split("$")
				if Anew.find("!") < 0:
					Ao['Aob'+str(j)] = float(Anew)
					fittingParametersList.append('Aob'+str(j))
					
				else:
					Ao['Aob'+str(j)] = float(Anew[:-1])						

				if Eanew.find("!") < 0:
					Ea['Eab'+str(j)] = float(Eanew)
					fittingParametersList.append('Eab'+str(j))
					
				else:
					Ea['Eab'+str(j)] = float(Eanew[:-1])

				arrBackward.append(j)

			else:

				if str(reaction_info.iloc[j,2]).find("!") < 0:
					kinetic_parameters['kb'+str(j)] = float(reaction_info.iloc[j,2])
					fittingParametersList.append('kb'+str(j))

				else:
					new_value = float(reaction_info.iloc[j,2][:-1])
					kinetic_parameters['kb'+str(j)] = new_value
		else:
			pass

	kin_in = kinetic_parameters.copy()
	Ao_in = Ao.copy()
	Ea_in = Ea.copy()
	Ga_in = Ga.copy()
	dG_in = dG.copy()

	return reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,fittingParametersList,arrForward,arrBackward


def solverIteration(time_step,method,solver,dk,dec_tim,inc_tim):
	
	"""
	Old time step iteration method. Could be useful in future implementations when runge-kutta stepping is optional.
	"""

	try:
		if method == 'simple_adaptive':
	
			u_temp.assign(u)
			uout_1 = call_solver(dk.assign(time_step/dec_tim),u_temp,u,u_1,solver)
			uout_3 = call_solver(dk.assign(time_step*inc_tim),u_temp,u,u_3,solver)
			uout_2 = call_solver(dk.assign(time_step),u_temp,u,u_2,solver,keep_sol=True)
			time_step = norm_comp(uout_1,uout_2,uout_3,dec_tim,inc_tim)

			return time_step
		
		elif method == 'None':
			solver.solve(annotate = False) # ### can pass annotate = False if I don't want it to record the solution
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

def fluxGeneration(reactor,gasses,reactants,pulse_size,Diff,voidage,dx,radius,dx2_r,outScale):
	
	"""Scale the output values of flux with the constant found here (messy now due to different trials)"""

	to_flux = []

	if reactor == 'tap':

		for k in range(0,gasses):
			#to_flux.append( (Diff[k][0]/(dx*Diff[4][0])) ) 
			if outScale.lower() == 'true':

				#!#!#!#!#!#!#!#to_flux.append( (Diff[k][0]*voidage[0] /(dx)) * (radius**2)*3.14159 / pulse_size ) #0.53*
				to_flux.append( 2*(Diff[k][0] /(dx)) * (radius**2)*3.14159 / pulse_size ) #0.53*   ## 
			else:
				#!#!#!#!#!#!#!#to_flux.append((Diff[k][0]*voidage[0] /(dx)) * (radius**2)*3.14159)
				to_flux.append(2*(Diff[k][0] /(dx)) * (radius**2)*3.14159 )
			#to_flux.append(2 *(dx*(radius**2)*3.14159) * (Diff[k][0] /(dx*voidage[0])))#(1/((1)*pulse_size)) *###??? changed from 1 to the new form
			#to_flux.append(2*Diff[k][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r))#(1/((1)*pulse_size)) *
			#to_flux.append(2*(1/((1+reactants)*pulse_size)) *Diff[k][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r))
		#to_flux.append( (Diff[gasses][0]*voidage[0] /(dx)) * (radius**2)*3.14159/(pulse_size) )
		#to_flux.append((2*Diff[gasses][0] * (dx*(radius**2)*3.14159)/(dx*voidage[0])))#*(1/((1)*pulse_size)) *
		#to_flux.append((2*(1/((1+reactants)*pulse_size)) *Diff[gasses][0] * dx*(radius**2)*3.14159/(voidage[0]*dx2_r)))
	elif reactor == 't_pfr' or 't_pfr_diff':
		print('to flux generation')
		for k in range(0,gasses):
			to_flux.append(1)
		
		to_flux.append(1)
	
	return to_flux

def defineBCs(reactor,elem_list,nec_values,V_sec,reacs_num,all_mol,reac_ratio,L_bound,R_bound,number_of_inerts):

	"""Define the appropriate boundary conditions for all monitored species"""

	# Zero Flux boundary condition at entrance of reactor (dcdx = 0 @ L = 0)
	if reactor == 'tap':
		bcs = []

		if elem_list != ['INERT_ONLY']:

			for k in range(0,nec_values):
				bcs.append(DirichletBC(V_sec.sub(k),Constant(0),R_bound))
			
			for k in range(0,int(number_of_inerts)):
				bcs.append(DirichletBC(V_sec.sub(all_mol-(1+k)),Constant(0),R_bound))
			
		else:	
			bcs.append(DirichletBC(V_sec,Constant(0),R_bound))
	
	# Vacuum BC at outlet of reactor (C = 0 @ L = L_reactor)
	elif reactor == 't_pfr' or 't_pfr_diff':
		bcs = []
		newValues = reac_ratio.split(',')
		for k in range(0,len(newValues)):
			bcs.append(DirichletBC(V_sec.sub(k),Constant(int(newValues[k])),L_bound))
		
		#for kin in range(int(reacs_num),len(newValues)-int(number_of_inerts)):
		#	bcs.append(DirichletBC(V_sec.sub(kin),Constant(0),L_bound))
		
		#for kfin in range(len(newValues),int(number_of_inerts)):
		#	bcs.append(DirichletBC(V_sec.sub(all_mol),Constant(0),L_bound))
	
	return bcs


def initializeVariableDictionaries(nec,moles,V_nu,u_nu,un_nu):
	
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


def establishMesh(in_1,cat,in_2,mesh_size):

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

def generateGraphAgain(reactor,gas_phase,reacs,inerts,scaleGraph):

	"""Generate the outlet flux graph (initialize, anyway)"""

	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('$Time\ (s)$', fontsize = 14)
	if reactor == 'tap':
		if scaleGraph.lower() == 'true':
			ax2.set_ylabel('$Outlet\ Flow\ (1/s)$', fontsize = 14)
		else:
			ax2.set_ylabel('$Outlet\ Flow\ (nmol/s)$', fontsize = 14)
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

def establishOutletGraph(reactor,gas_phase,reacs,inerts,scaleGraph):

	"""Generate the outlet flux graph (initialize, anyway)"""

	fig2, ax2 = plt.subplots()
	ax2.set_xlabel('$Time\ (s)$', fontsize = 14)
	if reactor == 'tap':
		if scaleGraph.lower() == 'true':
			ax2.set_ylabel('$Outlet\ Flow\ (1/s)$', fontsize = 14)
		else:
			ax2.set_ylabel('$Outlet\ Flow\ (nmol/s)$', fontsize = 14)
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

def knudsenTest(species_list,sim_steps,folder,time,points,intensity,fraction):

	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'

	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	print()
	print('Knudsen Regime Fingerprint Test (should be ~ 0.31):')
	curve_fitting = {}
	exp_data = user_data
	for k_num_new, k_new in enumerate(species_list):
		peak_loc = user_data[k_new].iloc[user_data[k_new][1].idxmax()]
		near_peak = peak_loc[0]/(time/sim_steps)
		peak2 = user_data[k_new].loc[user_data[k_new][0] == peak_loc[0]].index

		print(k_new+': '+str(peak_loc[0]*peak_loc[1]/(intensity*float(fraction[k_num_new]))))
		print()

####Fit every point
def curveFitting(species_list,sim_steps,folder,timeTot,points,objSpecies):
	frequency = 5
	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = timeTot/sim_steps

	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1
		
	user_data = {}

	if species_list[0].find('Inert') == 0:
				
		for klabel,k in enumerate(objSpecies):

			if objSpecies[klabel] == '1':

				user_data[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
		
	else:
		species_list = species_list[:len(species_list)]#'./experimental_data/'

		for k in range(0,len(species_list)):
			if objSpecies[k] == '1':
				user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	
	curve_fitting = {}
	exp_data = user_data
	
	for k_newNum, k_new in enumerate(species_list):

		if objSpecies[k_newNum] == '1':

			def find_experimental_point(n,exp_step):
				""" Find an appropriate intensity point for the fitting process """

				approx_exp_n = n*(syn_time_step)/exp_step
				if approx_exp_n != n:
					high = math.ceil(approx_exp_n)
					low = int(approx_exp_n)
					return interp(user_data[k_new][1][high-1],user_data[k_new][1][low-1],user_data[k_new][0][high-1],user_data[k_new][0][low-1],n*(syn_time_step))

				else:
					user_data[k_new][1][n]
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

			exp_time_step = user_data[k_new][0][user_data[k_new].shape[0]-2] - user_data[k_new][0][user_data[k_new].shape[0]-3]
			near_start = round(user_data[k_new].iloc[30,0],6)/(timeTot/sim_steps)
			
			for k in range(15,int(sim_steps),frequency):
				time_step.append(k)
				times.append(k*(syn_time_step))
				#print(exp_time_step)

				if round(exp_time_step,5) == syn_time_step:
					values.append(user_data[k_new][1][k])
				else:
					values.append(find_experimental_point(k,exp_time_step))#k*(syn_time_step)
				
				#print(find_experimental_point(k*(syn_time_step),exp_time_step))

			data = {}
			data['time_step'] = time_step
			data['times'] = times
			data['values'] = values

			curve_fitting[k_new] = data

	return curve_fitting


def pointFitting(species_list,sim_steps,folder,time,points,objSpecies):

	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = time/sim_steps
	
	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1

	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'

	for k in range(0,len(species_list)):
		if objSpecies[k] == '1':
			user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	
	curve_fitting = {}
	exp_data = user_data

	for k_newNum, k_new in enumerate(species_list):
		if objSpecies[k_newNum] == '1':
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
			exp_time_step = user_data[k_new][0][1] - user_data[k_new][0][0]

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
				value_test = 0.9*peak_loc[1]
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
				values.append(find_experimental_point(int(near_3),exp_time_step))
			
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
				values.append(find_experimental_point(int(near_4),exp_time_step))
				

			data = {}
			data['time_step'] = time_step
			data['times'] = times
			data['values'] = values

			curve_fitting[k_new] = data
			
	return curve_fitting


def generateGif(molecules,exp_loc,fit_loc,all_steps,constants,reactions,time_data):
	
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
			ax.scatter(exp_data[k_names][0], exp_data[k_names][1]/(exp_data[k_names][1].max()),label="Exp. "+k_names,alpha=0.3)
		
		for k_names in molecules:
			ax.plot(sim_data[k_names][0], sim_data[k_names][1]/(exp_data[k_names][1].max()),label="Syn. "+k_names,ls='--')
	
		ax.set_xlabel('Time (s)', fontsize=16)
		ax.set_ylabel('Normalized Flow (1/s)', fontsize=16)
		#textstr = 'Rate Constants: '
		
		#for k_len in range(0,len(constants[0])):
		#	textstr = '\n'.join((textstr,'k'+str(1+k_len)+': '+'{:0.3e}'.format(constants[step][k_len])))

		props = dict(facecolor='white')
		#ax.text(0.8, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props,ha='center')
		#textstr = 'Elementary Reactions:'
		
		#for k_reacs in range(0,len(reactions)):
		#	textstr = '\n'.join((textstr,reactions[k_reacs]))

		#ax.text(0.5, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props,ha='center')
		ax.text(0.02, 0.95,'Iteration: '+str(step),transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props)
		peak_peak = 0

		#for k_names in molecules[:-1]:
		#	peak_loc = exp_data[k_names].iloc[exp_data[k_names][1].idxmax()]
		#	plt.plot(peak_loc[0], peak_loc[1], 'ro')
		#	if peak_loc[1] > peak_peak:
		#		peak_peak = peak_loc[1]

		#	peak2 = exp_data[k_names].loc[exp_data[k_names][0] == exp_data[0]].index
		#	test3 = int(round((peak2[0]+1)/2,0))
		#	mid_loc = exp_data[k_names].iloc[test3,:]
		#	plt.plot(mid_loc[0], mid_loc[1], 'ro')

		#	#peak_loc = exp_data_2.iloc[exp_data_2[1].idxmax()]
		#	#peak2 = exp_data_2.loc[exp_data_2[0] == peak_loc[0]].index
		#	#plt.plot(peak_loc[0], peak_loc[1], 'ro')

		ax.legend(title='Gas Species',loc='center left')
		ax.set_xscale('log')
		ax.set_xlim(0.002,1)
		ax.set_ylim(0,1.5)

		#subpos = [0.4,0.13,0.5,0.4]
		#subax1 = add_subplot_axes(ax,subpos)
		#subax1.plot(x_data[:step],y_data[:step])
		#subax1.set_title('Time per Iteration')
		#subax1.set_ylabel('Time (minutes)')
		#subax1.set_xlabel('Iteration #')
		#subax1.set_xlim(0,all_steps)
		#subax1.set_ylim(0,max(y_data)*1.1)
		
		fig.canvas.draw()
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
		return image
	
	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	components = list(range(0,all_steps))
	for zoo in range(0,10):
		components.append(all_steps-1) 

	imageio.mimsave(fit_loc+'/output.gif', [tap_plot(i) for i in components], fps=4)


"""Functions used to keep output organized"""

def generateFolder(path_name):
	try:  
		os.mkdir(path_name)
	except OSError:  
		###print ("Creation of the directory %s failed" % path_name)
		pass
	else:  
		###print ("Successfully created the directory %s " % path_name)
		pass
	
def storeSens(yes_no,output_file,gasses,legend_ref):	
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

def storeDataFunc(yes_no,output_file):
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

def callSolver(dk,u_temp,u,u_new,solver_def,keep_sol = True):
	u_temp.assign(u)
	solver_def.solve()
	u_new.assign(u)
	if keep_sol == False:
		u.assign(u_temp)
	return u_new

def normComp(u1,u2,u3,d_t,i_t):
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

def processTime(start_time):
	if (time.time() - start_time) < 120:
		return 'Completed in: '+str(round((time.time() - start_time),3))+' seconds'
	elif (time.time() - start_time)/60 < 120:
		return 'Completed in: '+str(round((time.time() - start_time)/60,3))+' minutes'
	else:
		return 'Completed in: '+str(round((time.time() - start_time)/3600,3))+' hours'

def evalCB(j, m):
	print('eval')
	print(j)
	print(m)


def errorOutput(elem_reacs):
	
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

def molecularProperties(gasSpecies,propValue):

	molecularValues = {}

	# Oxygen Scrambling
	molecularValues['O218'] = {'mass':36}
	molecularValues['O216'] = {'mass':32}
	molecularValues['O2'] = {'mass':32}
	molecularValues['O18O16'] = {'mass':34}
	
	# Natural Gas
	molecularValues['CO'] = {'mass':28.01}
	molecularValues['CO2'] = {'mass':44.01}
	molecularValues['CH2'] = {'mass':14.01}
	molecularValues['CH3'] = {'mass':15.01}
	molecularValues['CH4'] = {'mass':16.01}
	molecularValues['C2H4'] = {'mass':32}
	molecularValues['C2H6'] = {'mass':34}
	molecularValues['C3H4'] = {'mass':42.03}
	molecularValues['C3H6'] = {'mass':42.03}
	molecularValues['C3H8'] = {'mass':44.03}
	molecularValues['C4H6'] = {'mass':54.04}
	molecularValues['C4H8'] = {'mass':56.04}
	molecularValues['C4H10'] = {'mass':58.04}
	molecularValues['C5H10'] = {'mass':70.05}
	molecularValues['C5H12'] = {'mass':72.05}
	molecularValues['H2S'] = {'mass':34}
	molecularValues['H2O'] = {'mass':18}
	molecularValues['SO2'] = {'mass':64}

	# Atmospheric / Exhaust / Fertilizer
	molecularValues['N2'] = {'mass':28}
	molecularValues['NO'] = {'mass':30}
	molecularValues['NO2'] = {'mass':46}
	molecularValues['NO3'] = {'mass':62}
	molecularValues['N2O5'] = {'mass':108}
	molecularValues['NH2'] = {'mass':16}
	molecularValues['NH3'] = {'mass':17}
	molecularValues['HNO3'] = {'mass':63}
	molecularValues['O3'] = {'mass':48}

	# Other diatomic / common species
	molecularValues['Cl2'] = {'mass':71}
	molecularValues['Br2'] = {'mass':159.8}
	molecularValues['I2'] = {'mass':252}
	molecularValues['H2'] = {'mass':2.01}
	molecularValues['He'] = {'mass':4}
	molecularValues['Ar'] = {'mass':40}
	molecularValues['HCl'] = {'mass':36.5}

	return molecularValues[gasSpecies][propValue]