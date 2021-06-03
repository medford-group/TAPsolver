# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

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

def readInput(sim_file,inputForm = 'old'):
	
	"""
	Convert the input file into dictionaries for TAPsolver to use
	"""
	if inputForm == 'old':
	
		user_data = pd.read_csv(sim_file,header=None)
		
		rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
		rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(user_data == 'Reaction_Information')
	
		reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
		feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
		reaction_info = user_data.iloc[1+rows_4[0]:,:]
	
		reactor_kinetics_input = {}
		
		for k in range(0,len(reactor_info.index)):
			try:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
	
	
		for k in range(0,len(feed_surf_info.index)):
			try:
				reactor_kinetics_input[feed_surf_info.iloc[k,0]] = float(feed_surf_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[feed_surf_info.iloc[k,0]] = feed_surf_info.iloc[k,1]
	
	
		reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']*(reactor_kinetics_input['Catalyst Location'] - 0.5*(reactor_kinetics_input['Catalyst Fraction']))#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
		reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length'] - (reactor_kinetics_input['len_cat'] + reactor_kinetics_input['len_inert_1'])#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	
		reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()

	else:
		
		user_data = pd.read_csv(sim_file,header=None)
		
		rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
		rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(user_data == 'Reaction_Information')

		linkedKinetics = False
		if user_data[0].str.contains('Linked Kinetics').any():
			linkedKinetics = True
			rows_5, cols_5 = np.where(user_data == 'Linked Kinetics')
		thermoConstraints = False
		if user_data[0].str.contains('Thermodynamic Constraints').any() and linkedKinetics == True:
			thermoConstraints = True
			rows_6, cols_6 = np.where(user_data == 'Thermodynamic Constraints')    		
		elif user_data[0].str.contains('Thermodynamic Constraints').any():
			thermoConstraints = True
			rows_5, cols_5 = np.where(user_data == 'Thermodynamic Constraints')    		
		reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:]
		feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
	#	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]
		if thermoConstraints == False and linkedKinetics == False:
			reaction_info = user_data.iloc[1+rows_4[0]:,:]
		elif linkedKinetics == True and thermoConstraints == False:
			reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			linked_kinetics = user_data.iloc[1+rows_5[0]:,:]
		elif linkedKinetics == False and thermoConstraints == True:		
			reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			thermo_constraints = user_data.iloc[1+rows_5[0]:,:]
		elif linkedKinetics == True and thermoConstraints == True:
			reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			linked_kinetics = user_data.iloc[(1+rows_5[0]):(rows_6[0]-1),:]
			thermo_constraints = user_data.iloc[1+rows_6[0]:,:]		


		#thermoConstraints = False
		#if user_data[0].str.contains('Thermodynamic Constraints').any():
		#	thermoConstraints = True
		#	rows_5, cols_5 = np.where(user_data == 'Thermodynamic Constraints')    		
		#
		#reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
		#feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
		#
		#if thermoConstraints == False:
		#	reaction_info = user_data.iloc[1+rows_4[0]:,:]
		#else:
		#	reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
		#	thermo_constraints = user_data.iloc[1+rows_5[0]:,:]

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
				
			else:
				reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])+','
				reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])+','
				reactor_kinetics_input['Mass List'] = reactor_kinetics_input['Mass List']+str(feed_surf_info.iloc[3,1+jnum])+','
			reactor_kinetics_input['Number of Reactants'] += 1
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
	link = {}

	fittingParametersList = []

	gForward = []
	arrForward = []
	arrBackward = []
	linkForward = []
	linkBackward = []

	for j in range(0,len(reaction_info.index)):
		if reaction_info.iloc[j,1].find("#") > 0:
			Anew, Eanew = reaction_info.iloc[j,1].split("#")
			if Anew.find("!") < 0:
				#!?
				if Anew.find("{") == 0:
					link['Ga'+str(j)] = float(Anew)
					fittingParametersList.append('Ga'+str(j))
				else:
					Ga['Ga'+str(j)] = float(Anew)
					kinetic_parameters['Ga'+str(j)] = float(Anew)
					fittingParametersList.append('Ga'+str(j))
				
			else:
				#!?
				if Anew.find("{") == 0:
					fittingParametersList.append(Anew)
				else:
					Ga['Ga'+str(j)] = float(Anew[:-1])
			if Eanew.find("!") < 0:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					dG['dG'+str(j)] = float(Eanew)
					fittingParametersList.append('dG'+str(j))
			
			else:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					dG['dG'+str(j)] = float(Eanew[:-1])
			gForward.append(j)
		elif reaction_info.iloc[j,1].find("$") > 0:
			Anew, Eanew = reaction_info.iloc[j,1].split("$")
			if Anew.find("!") < 0:
				#!?
				if Anew.find("{") == 0:
					fittingParametersList.append(Anew)
				else:	
					Ao['Aof'+str(j)] = float(Anew)
					fittingParametersList.append('Aof'+str(j))
			
			else:
				#!?
				if Anew.find("{") == 0:
					fittingParametersList.append(Anew)
				else:
					Ao['Aof'+str(j)] = float(Anew[:-1])
			if Eanew.find("!") < 0:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					Ea['Eaf'+str(j)] = float(Eanew)
					fittingParametersList.append('Eaf'+str(j))
			
			else:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					Ea['Eaf'+str(j)] = float(Eanew[:-1])
			arrForward.append(j)
		else:
			if reaction_info.iloc[j,1].find("!") < 0:
				if reaction_info.iloc[j,1].find("{") == 0:
					link['kf'+str(j)] = reaction_info.iloc[j,1][1:-1]
					if reaction_info.iloc[j,1] not in fittingParametersList:
						fittingParametersList.append(reaction_info.iloc[j,1])
					linkForward.append(j)
				else:
					kinetic_parameters['kf'+str(j)] = float(reaction_info.iloc[j,1])
					fittingParametersList.append('kf'+str(j))
			
			else:
				if reaction_info.iloc[j,1].find("{") == 0:
					link['kf'+str(j)] = reaction_info.iloc[j,1][1:-1]
					if reaction_info.iloc[j,1] not in fittingParametersList:
						fittingParametersList.append(reaction_info.iloc[j,1])
					linkForward.append(j)
				else:
					new_value = float(reaction_info.iloc[j,1])
					kinetic_parameters['kf'+str(j)] = new_value#float(reaction_info.iloc[j,1])
		if str(reaction_info.iloc[j,2]) != 'nan':
			if str(reaction_info.iloc[j,2]).find("$") > 0:
				Anew, Eanew = str(reaction_info.iloc[j,2]).split("$")
				if Anew.find("!") < 0:
					#!?
					if Anew.find("{") == 0:
						fittingParametersList.append(Anew)
					else:
						Ao['Aob'+str(j)] = float(Anew)
						fittingParametersList.append('Aob'+str(j))
					
				else:
					#!?
					if Anew.find("{") == 0:
						fittingParametersList.append(Anew)
					else:
						Ao['Aob'+str(j)] = float(Anew[:-1])						
				if Eanew.find("!") < 0:
					#!?
					if Eanew.find("{") == 0:
						fittingParametersList.append(Eanew)
					else:
						Ea['Eab'+str(j)] = float(Eanew)
						fittingParametersList.append('Eab'+str(j))
					
				else:
					#!?
					if Eanew.find("{") == 0:
						fittingParametersList.append(Eanew)
					else:
						Ea['Eab'+str(j)] = float(Eanew[:-1])
				arrBackward.append(j)
			else:
				if str(reaction_info.iloc[j,2]).find("!") < 0:
				
					if reaction_info.iloc[j,2].find("{") == 0:
						link['kb'+str(j)] = reaction_info.iloc[j,2][1:-1]
						if reaction_info.iloc[j,2] not in fittingParametersList:
							fittingParametersList.append(reaction_info.iloc[j,2])
						linkBackward.append(j)
					else:
						kinetic_parameters['kb'+str(j)] = float(reaction_info.iloc[j,2])
						fittingParametersList.append('kb'+str(j))
				else:
					if reaction_info.iloc[j,2].find("{") == 0:
						link['kb'+str(j)] = reaction_info.iloc[j,2][1:-2]
						if reaction_info.iloc[j,2] not in fittingParametersList:
							fittingParametersList.append(reaction_info.iloc[j,2])
						linkBackward.append(j)
					else:
						new_value = float(reaction_info.iloc[j,2][:-1])
						kinetic_parameters['kb'+str(j)] = new_value
		else:
			pass

#	for j in range(0,len(reaction_info.index)):
#		if reaction_info.iloc[j,1].find("#") > 0:
#			Anew, Eanew = reaction_info.iloc[j,1].split("#")
#			if Anew.find("!") < 0:
#				Ga['Ga'+str(j)] = float(Anew)
#				kinetic_parameters['Ga'+str(j)] = float(Anew)
#				fittingParametersList.append('Ga'+str(j))
#			
#			else:
#				Ga['Ga'+str(j)] = float(Anew[:-1])
#
#			if Eanew.find("!") < 0:
#				dG['dG'+str(j)] = float(Eanew)
#				fittingParametersList.append('dG'+str(j))
#			
#			else:
#				dG['dG'+str(j)] = float(Eanew[:-1])
#
#			gForward.append(j)
#
#		elif reaction_info.iloc[j,1].find("$") > 0:
#			Anew, Eanew = reaction_info.iloc[j,1].split("$")
#			if Anew.find("!") < 0:
#				Ao['Aof'+str(j)] = float(Anew)
#				fittingParametersList.append('Aof'+str(j))
#			
#			else:
#				Ao['Aof'+str(j)] = float(Anew[:-1])
#
#			if Eanew.find("!") < 0:
#				Ea['Eaf'+str(j)] = float(Eanew)
#				fittingParametersList.append('Eaf'+str(j))
#			
#			else:
#				Ea['Eaf'+str(j)] = float(Eanew[:-1])
#
#			arrForward.append(j)
#
#		else:
#			if reaction_info.iloc[j,1].find("!") < 0:
#				kinetic_parameters['kf'+str(j)] = float(reaction_info.iloc[j,1])
#				fittingParametersList.append('kf'+str(j))
#			
#			else:
#				new_value = float(reaction_info.iloc[j,1][:-1])
#				kinetic_parameters['kf'+str(j)] = new_value#float(reaction_info.iloc[j,1])
#
#		if str(reaction_info.iloc[j,2]) != 'nan':
#			if str(reaction_info.iloc[j,2]).find("$") > 0:
#				Anew, Eanew = str(reaction_info.iloc[j,2]).split("$")
#				if Anew.find("!") < 0:
#					Ao['Aob'+str(j)] = float(Anew)
#					fittingParametersList.append('Aob'+str(j))
#					
#				else:
#					Ao['Aob'+str(j)] = float(Anew[:-1])						
#
#				if Eanew.find("!") < 0:
#					Ea['Eab'+str(j)] = float(Eanew)
#					fittingParametersList.append('Eab'+str(j))
#					
#				else:
#					Ea['Eab'+str(j)] = float(Eanew[:-1])
#
#				arrBackward.append(j)
#
#			else:
#
#				if str(reaction_info.iloc[j,2]).find("!") < 0:
#					kinetic_parameters['kb'+str(j)] = float(reaction_info.iloc[j,2])
#					fittingParametersList.append('kb'+str(j))
#
#				else:
#					new_value = float(reaction_info.iloc[j,2][:-1])
#					kinetic_parameters['kb'+str(j)] = new_value
#		else:
#			pass
#
	kin_in = kinetic_parameters.copy()
	Ao_in = Ao.copy()
	Ea_in = Ea.copy()
	Ga_in = Ga.copy()
	dG_in = dG.copy()

	linkedParameters = {}
	if linkedKinetics == True:
		for j in range(0,len(linked_kinetics.index)):
			linkedParameters[linked_kinetics.iloc[j,0]] = float(linked_kinetics.iloc[j,1])
	reactor_kinetics_input['linked parameters'] = linkedParameters
	reactor_kinetics_input['linked names'] = link
	reactor_kinetics_input['link forward'] = linkForward
	reactor_kinetics_input['link backward'] = linkBackward

	thermo_equations = []
	thermo_values = []
	if thermoConstraints == True:
		for j in range(0,len(thermo_constraints.index)):
			thermo_equations.append(thermo_constraints.iloc[j,0])
			try:
				if type(float(thermo_constraints.iloc[j,1])) == float and math.isnan(float(thermo_constraints.iloc[j,1])) == False:
					thermo_values.append(thermo_constraints.iloc[j,1])
				else:
					thermo_values.append('none')
			except:
				print(type(thermo_constraints.iloc[j,1]))
				pass

	reactor_kinetics_input['thermo equations'] = thermo_equations
	reactor_kinetics_input['thermo values'] = thermo_values

	return reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,fittingParametersList,arrForward,arrBackward

def readBatchInput(sim_file):

	user_data = pd.read_csv(sim_file,header=None)
		
	rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
	rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	rows_4, cols_4 = np.where(user_data == 'Reaction_Information')

	thermoConstraints = False
	if user_data[0].str.contains('Thermodynamic Constraints').any():
		thermoConstraints = True
		rows_5, cols_5 = np.where(user_data == 'Thermodynamic Constraints')    		

	reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
#	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]

	if thermoConstraints == False:
		reaction_info = user_data.iloc[1+rows_4[0]:,:]
	else:
		reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
		thermo_constraints = user_data.iloc[1+rows_5[0]:,:]

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
	surface_species = feed_surf_info.iloc[5,1:number_of_surface+1]
		
	reactor_kinetics_input['Pulse Size'] = ''
	reactor_kinetics_input['Pulse Time'] = ''
	reactor_kinetics_input['Number of Reactants'] = 0

	for jnum,j in enumerate(gas_species):
			
		if j.find('Inert') == 0:
			reactor_kinetics_input['Number of Inerts'] += 1
		if jnum == len(gas_species)-1:

			reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])
			reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])
				
		else:
			reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])+','
			reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])+','
		reactor_kinetics_input['Number of Reactants'] += 1
	reactor_kinetics_input['Initial Surface Composition'] = ''

	for jnum,j in enumerate(surface_species):

		if jnum == len(surface_species)-1:
			reactor_kinetics_input['Initial Surface Composition'] = reactor_kinetics_input['Initial Surface Composition']+str(feed_surf_info.iloc[5,1+jnum])
		else:
			reactor_kinetics_input['Initial Surface Composition'] = reactor_kinetics_input['Initial Surface Composition']+str(feed_surf_info.iloc[5,1+jnum])+','

	for k in range(0,len(reactor_info.index)):
		try:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1])
		except ValueError:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
			
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
				kinetic_parameters['kf'+str(j)] = new_value

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

	linkedParameters = {}
	if linkedKinetics == True:
		for j in range(0,len(linked_kinetics.index)):
			linkedParameters[linked_kinetics.iloc[j,0]] = float(linked_kinetics.iloc[j,1])
	reactor_kinetics_input['linked parameters'] = linkedParameters
	reactor_kinetics_input['linked names'] = link
	reactor_kinetics_input['link forward'] = linkForward
	reactor_kinetics_input['link backward'] = linkBackward

	thermo_equations = []
	thermo_values = []
	if thermoConstraints == True:
		for j in range(0,len(thermo_constraints.index)):
			thermo_equations.append(thermo_constraints.iloc[j,0])
			try:
				if type(float(thermo_constraints.iloc[j,1])) == float and math.isnan(float(thermo_constraints.iloc[j,1])) == False:
					thermo_values.append(thermo_constraints.iloc[j,1])
				else:
					thermo_values.append('none')
			except:
				print(type(thermo_constraints.iloc[j,1]))
				pass

	reactor_kinetics_input['thermo equations'] = thermo_equations
	reactor_kinetics_input['thermo values'] = thermo_values
	return reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,fittingParametersList,arrForward,arrBackward



def solverIteration(time_step,method,solver,dk,dec_tim,inc_tim):

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

def fluxGeneration(reactor,gasses,reactants,pulse_size,Diff,voidage,dx,radius,dx2_r,outScale):
	
	"""Scale the output values of flux with the constant found here (messy now due to different trials)"""

	to_flux = []

	if reactor == 'tap':

		for k in range(0,gasses):
			if outScale.lower() == 'true':
				to_flux.append( 2*(Diff[k][0] /(dx)) * (radius**2)*3.14159 / pulse_size ) #0.53*   ## 
			else:
				to_flux.append(2*(Diff[k][0] /(dx)) * (radius**2)*3.14159 )
	elif reactor == 't_pfr' or 't_pfr_diff':
		print('to flux generation')

		for k in range(0,gasses):
			to_flux.append(1)
		
		to_flux.append(1)
	
	return to_flux

def defineBCs(reactor,elem_list,nec_values,V_sec,reacs_num,all_mol,reac_ratio,L_bound,R_bound,number_of_inerts):

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
	
	# Vacuum BC at outlet of reactor (C = 0 @ L = L_reactor)
	elif reactor == 't_pfr' or 't_pfr_diff':
		bcs = []
		newValues = reac_ratio.split(',')
		for k in range(0,len(newValues)):
			bcs.append(DirichletBC(V_sec.sub(k),Constant(int(newValues[k])),L_bound))
	
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

	elif reactor == 'batch':
		print('batch')
		ax2.set_ylabel('$Concentration$')

	legend_label = []
	header = "time"
	for k in range(0,gas_phase):
		legend_label.append(reacs[k])
		header = header+","+reacs[k]
	if reactor != 'batch':
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
	frequency = 1
	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = timeTot/sim_steps

	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1
		
	user_data = {}

	if species_list[0].find('Inert') == 0:
				
		for klabel,k in enumerate(objSpecies):

			if objSpecies[klabel] == '1':
				try:
					fitStartValue = False
					user_data[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
				except:
					fitStartValue = True
					user_data[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
	else:
		species_list = species_list[:len(species_list)]#'./experimental_data/'
		for k in range(0,len(species_list)):
			if objSpecies[k] == '1':
				try:
					fitStartValue = False
					user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
				except:
					fitStartValue = True
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
			near_start = round(user_data[k_new].iloc[2,0],6)/(timeTot/sim_steps)
			if fitStartValue == True:
				fitStartTime = 0
				frequency = 1
			else: 
				fitStartTime = 2
				frequency = 1

			for k in range(fitStartTime,int(sim_steps),frequency):
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

####Fit every point
def stdEstablishment(species_list,sim_steps,folder,timeTot,points,objSpecies,stdValue=0.0): # uniform or varied
	frequency = 3
	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = timeTot/sim_steps
	
	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1
		
	user_data = {}
	user_data_2 = {}

	if species_list[0].find('Inert') == 0:
				
		for klabel,k in enumerate(objSpecies):

			if stdValue == 0.0:
				if objSpecies[klabel] == '1':
					#try:
					fitStartValue = False
					user_data[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'_std.csv',header=None)
	else:
		species_list = species_list[:len(species_list)]#'./experimental_data/'
		for k in range(0,len(species_list)):
			if stdValue == 0.0:
				if objSpecies[k] == '1':
					fitStartValue = False
					user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'_std.csv',header=None)

	if species_list[0].find('Inert') == 0:
				
		for klabel,k in enumerate(objSpecies):

			if objSpecies[klabel] == '1':
				try:
					fitStartValue = False
					user_data_2[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
				except:
					fitStartValue = True
					user_data_2[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
	else:
		species_list = species_list[:len(species_list)]#'./experimental_data/'
		for k in range(0,len(species_list)):
			if objSpecies[k] == '1':
				try:
					fitStartValue = False
					user_data_2[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
				except:
					fitStartValue = True
					user_data_2[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	
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

			time_step = []
			times = []
			values = []

			exp_time_step = user_data_2[k_new][0][user_data_2[k_new].shape[0]-2] - user_data_2[k_new][0][user_data_2[k_new].shape[0]-3]
			near_start = round(user_data_2[k_new].iloc[30,0],6)/(timeTot/sim_steps)
			if fitStartValue == True:
				fitStartTime = 0
				frequency = 10
			else: 
				fitStartTime = 30
				frequency = 10

			if stdValue == 0.0:
				for k in range(fitStartTime,int(sim_steps),frequency):
					time_step.append(k)
					times.append(k*(syn_time_step))
					#print(exp_time_step)

					if round(exp_time_step,5) == syn_time_step:
						values.append(user_data[k_new][1][k])
					else:
						values.append(find_experimental_point(k,exp_time_step))#k*(syn_time_step)
					
			else:
				for k in range(fitStartTime,int(sim_steps),frequency):
					time_step.append(k)
					times.append(k*(syn_time_step))
					values.append(stdValue)

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
				
				if approx_exp_n != n:
					high = math.ceil(approx_exp_n)
					low = int(approx_exp_n)
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


def generateGif(molecules,exp_loc,fit_loc,all_steps,constants,reactions,time_data,xscale,yscale):
	
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
			if yscale == 'normalized':
				ax.scatter(exp_data[k_names][0], exp_data[k_names][1]/(exp_data[k_names][1].max()),label="Exp. "+k_names,alpha=0.3)
			else:
				ax.scatter(exp_data[k_names][0], exp_data[k_names][1],label="Exp. "+k_names,alpha=0.3)
		for k_names in molecules:
			if yscale == 'normalized':
				ax.plot(sim_data[k_names][0], sim_data[k_names][1]/(exp_data[k_names][1].max()),label="Syn. "+k_names,ls='--')
			else:
				ax.plot(sim_data[k_names][0], sim_data[k_names][1],label="Syn. "+k_names,ls='--')

		ax.set_xlabel('Time (s)', fontsize=16)
		if yscale == 'normalized':
			ax.set_ylabel('Normalized Flow (1/s)', fontsize=16)
		else:
			ax.set_ylabel('Flow (nmol/s)', fontsize=16)

		props = dict(facecolor='white')

		ax.text(0.42, 0.95,'Iteration: '+str(step),transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props)
		peak_peak = 0

		ax.legend(title='Gas Species',loc='upper right')
		if xscale == 'log':
			ax.set_xscale('log')
			ax.set_xlim(0.002,1)
		if yscale == 'normalized':
			ax.set_ylim(0,1.5)
		
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
		pass
	else:  
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

def molecularProperties(gasSpecies,propValue,temperature=398):

	def thermoCalculation(T,dH,A,B,C,D,E,F,G,H):
		temp = T/1000

		h = dH + (A*temp + (B*(temp**2)/2) + (C*(temp**3)/3) + (D*(temp**4)/4) - (E/temp) + F - H)#/1000
		s = A*math.log(temp) + B*temp + C*(temp**2)/2 + D*(temp**3)/3 - E/(2*temp**2) + G

		return h,s/1000

	molecularValues = {}

	## Nist Chemistry WebBook was used to determine the shomate equation parameters
	## 

	# Oxygen Scrambling
	molecularValues['O218'] = {'mass':36,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O216'] = {'mass':32,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O2'] = {'mass':32,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O18O16'] = {'mass':34,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	molecularValues['O16O18'] = {'mass':34,'shomate':{'dH':0.0,'A':31.32234,'B':-20.23531,'C':57.86644,'D':-36.50624,'E':-0.007374,'F':-8.903471,'G':246.7945,'H':0.0,'Trange':[100,700]}}
	
	# Natural Gas
	molecularValues['CO'] = {'mass':28.01,'shomate':{'dH':-110.5271,'A':25.5679,'B':6.096130,'C':4.05465,'D':-2.67130,'E':0.131021,'F':-118.0089,'G':227.366,'H':-110.5271,'Trange':[298,1300]}}
	molecularValues['CO2'] = {'mass':44.01,'shomate':{'dH':-393.5224,'A':24.99735,'B':55.18696,'C':-33.69137,'D':7.948387,'E':-0.136638,'F':-403.6075,'G':228.2431,'H':-393.5224,'Trange':[298,1200]}}
	#molecularValues['CH2'] = {'mass':14.01}
	#molecularValues['CH3'] = {'mass':15.01}
	molecularValues['CH4'] = {'mass':16.01,'shomate':{'dH':-74.87310,'A':-0.703029,'B':108.4773,'C':-42.52157,'D':5.862788,'E':0.678565,'F':-76.84376,'G':158.7163,'H':-74.87310,'Trange':[298,1300]}}
	molecularValues['C2H4'] = {'mass':32,'shomate':{'dH':52.46694,'A':-6.387880,'B':184.4019,'C':-112.9718,'D':28.49593,'E':0.315540,'F':48.17332,'G':163.1568,'H':52.46694,'Trange':[298,1200]}}
	#molecularValues['C2H6'] = {'mass':34,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C3H4'] = {'mass':42.03,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C3H6'] = {'mass':42.03,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C3H8'] = {'mass':44.03,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C4H6'] = {'mass':54.04,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C4H8'] = {'mass':56.04,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C4H10'] = {'mass':58.04,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C5H10'] = {'mass':70.05,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	#molecularValues['C5H12'] = {'mass':72.05,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	molecularValues['H2S'] = {'mass':34,'shomate':{'dH':-20.502,'A':26.88412	,'B':18.67809,'C':3.434203,'D':-3.378702,'E':0.135882	,'F':-28.91211,'G':233.3747,'H':-20.502,'Trange':[298,1400]}}
	#molecularValues['H2O'] = {'mass':18,'shomate':{'dH':,'A':,'B':,'C':,'D':,'E':,'F':,'G':,'H':,'Trange':[]}}
	molecularValues['SO2'] = {'mass':64,'shomate':{'dH':-296.8422,'A':21.43049,'B':74.35094,'C':-57.75217,'D':16.35534,'E':0.086731,'F':-305.7688,'G':254.8872,'H':-296.8422,'Trange':[298,1200]}}

	# Atmospheric / Exhaust / Fertilizer
	molecularValues['N2'] = {'mass':28,'shomate':{'dH':0.0,'A':28.98641,'B':1.853978,'C':-9.647459,'D':16.63537,'E':0.000117,'F':-8.671914,'G':226.4168,'H':0.0,'Trange':[100,500]}}
	molecularValues['NO'] = {'mass':30,'shomate':{'dH':90.29114,'A':23.83491,'B':12.58878,'C':-1.139011,'D':-1.497459,'E':0.214194,'F':83.35783,'G':237.1219,'H':90.29114,'Trange':[298,1200]}}
	molecularValues['NO2'] = {'mass':46,'shomate':{'dH':33.09502,'A':16.10857,'B':75.89525,'C':-54.38740,'D':14.30777,'E':0.239423,'F':26.17464,'G':240.5386,'H':33.09502,'Trange':[298,1200]}}
	molecularValues['N2O'] = {'mass':44.01,'shomate':{'dH':82.04824,'A':27.67988,'B':51.14898,'C':-30.64454,'D':6.847911,'E':-0.157906,'F':71.24934,'G':238.6164,'H':82.04824,'Trange':[298,1400]}}
	molecularValues['NO3'] = {'mass':62,'shomate':{'dH':71.12800,'A':11.22316	,'B':166.3889,'C':-148.4458,'D':47.40598,'E':-0.176791,'F':61.00858	,'G':221.7679,'H':71.12800,'Trange':[298,1100]}}
	molecularValues['N2O5'] = {'mass':108,'shomate':{'dH':82.84320,'A':39.09663,'B':114.8006,'C':-81.97125,'D':21.77249,'E':-0.088738,'F':66.46786,'G':324.5776,'H':82.84320,'Trange':[298,1100]}}
	molecularValues['NH3'] = {'mass':17.03,'shomate':{'dH':-45.89806,'A':19.99563,'B':49.77119,'C':-15.37599,'D':1.921168,'E':0.189174,'F':-53.30667,'G':203.8591,'H':-45.89806,'Trange':[298,1400]}}
	molecularValues['HNO3'] = {'mass':63,'shomate':{'dH':-134.3060,'A':19.63229,'B':153.9599,'C':-115.8378,'D':32.87955,'E':-0.249114,'F':-146.8818,'G':247.7049,'H':-134.3060,'Trange':[298,1200]}}
	molecularValues['O3'] = {'mass':48,'shomate':{'dH':142.6740,'A':21.66157,'B':79.86001,'C':-66.02603,'D':19.58363,'E':-0.079251,'F':132.9407,'G':243.6406,'H':142.6740,'Trange':[298,1200]}}
	molecularValues['N2H4'] = {'mass':32}


	# Other diatomic / common species
	molecularValues['Cl2'] = {'mass':71,'shomate':{'dH':0.0,'A':33.05060,'B':12.22940,'C':-12.06510,'D':4.385330,'E':-0.159494,'F':-10.83480,'G':259.0290,'H':0.0,'Trange':[298,1000]}}
	molecularValues['Br2'] = {'mass':159.8,'shomate':{'dH':30.91001,'A':38.52723,'B':-1.976835,'C':1.526107,'D':-0.198398,'E':-0.185815,'F':18.87620,'G':291.4863,'H':30.91001,'Trange':[333,3400]}}
	molecularValues['I2'] = {'mass':253.8,'shomate':{'dH':62.42110,'A':37.79763	,'B':0.225453,'C':-0.912556,'D':1.034913,'E':-0.083826,'F':50.86865,'G':305.9199,'H':62.42110,'Trange':[458,2000]}}
	molecularValues['H2'] = {'mass':2.01,'shomate':{'dH':0.0,'A':33.066178,'B':-11.363417,'C':11.432816,'D':-2.772874,'E':-0.158558,'F':-9.980797,'G':172.707974,'H':0.0,'Trange':[298,1000]}}
	molecularValues['HCl'] = {'mass':36.5,'shomate':{'dH':-92.312,'A':32.12392,'B':-13.45805,'C':19.86852,'D':-6.853936,'E':-0.049672,'F':-101.6206,'G':228.6866,'H':-92.312,'Trange':[298,1200]}}

	molecularValues['He'] = {'mass':4,'shomate':{'dH':0.0,'A':20.78603,'B':4.8506e-10,'C':-1.5829e-10,'D':1.5251e-11,'E':3.1963e-11,'F':-6.197341,'G':151.3064,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Ne'] = {'mass':20.18,'shomate':{'dH':0.0,'A':20.78603,'B':4.8506e-10,'C':-1.5829e-10,'D':1.525102e-11,'E':3.1963e-11,'F':-6.197341,'G':171.48,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Ar'] = {'mass':40,'shomate':{'dH':0.0,'A':20.786,'B':2.825911e-7,'C':-1.46419e-7,'D':1.092131e-8,'E':-3.661371e-8,'F':-6.197350,'G':179.999,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Kr'] = {'mass':83.798,'shomate':{'dH':0.0,'A':20.78603,'B':4.850638e-10,'C':-1.582916e-10,'D':1.525102e-11,'E':3.196347e-11,'F':-6.197341,'G':189.239,'H':0.0,'Trange':[298,6000]}}
	molecularValues['Xe'] = {'mass':131.293,'shomate':{'dH':0.0,'A':20.786,'B':7.4493e-7,'C':-2.0494e-7,'D':1.066e-8,'E':2.500261e-8,'F':-6.197350,'G':194.8380,'H':0.0,'Trange':[298,6000]}}

	if gasSpecies in molecularValues.keys():
		if propValue == 'mass':
			return molecularValues[gasSpecies][propValue]
		elif propValue == 'freeEnergy':

			tempShomate = molecularValues[gasSpecies]['shomate']
			hnew,snew = thermoCalculation(temperature,tempShomate['dH'],tempShomate['A'],tempShomate['B'],tempShomate['C'],tempShomate['D'],tempShomate['E'],tempShomate['F'],tempShomate['G'],tempShomate['H'])
			
			return hnew-temperature*snew
	else:
		return ''




def genInput0():
	f= open("input_file_0.csv","w+")
	f.write("Reactor_Information,Zone1,Zone2,Zone3,,,,,,\n")
	f.write("Zone Length,2,0.25,2,,,,,,\n")
	f.write("Zone Void,0.4,0.4,0.4,,,,,,\n")
	f.write("Reactor Radius,1,,,,,,,,\n")
	f.write("Reactor Temperature,385.65,,,,,,,,\n")
	f.write("Mesh Size,200,,,,,,,,\n")
	f.write("Catalyst Mesh Density,0,,,,,,,,\n")
	f.write("Output Folder Name,new_results,,,,,,,,\n")
	f.write("Experimental Data Folder,none,,,,,,,,\n")
	f.write("Reference Diffusion Inert,16,,,,,,,,\n")
	f.write("Reference Diffusion Catalyst,16,,,,,,,,\n")
	f.write("Reference Temperature,385.5,,,,,,,,\n")
	f.write("Reference Mass,40,,,,,,,,\n")
	f.write("Advection Value,0.0,,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Feed_&_Surface_Composition,,,,,,,,,\n")
	f.write(",CO,O2,Inert-1,,,,,,\n")
	f.write("Intensity,1,1,1,,,,,,\n")
	f.write("Time,0,0.1,0,,,,,,\n")
	f.write("Mass,28,28,28,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write(",CO*,O2*,*,,,,,,\n")
	f.write("Initial Concentration,0,0,30,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Reaction_Information,,,,,,,,,\n")
	f.write("CO + * <-> CO*,2,1,,,,,,,\n")
	f.write("O2 + * <-> O2*,2,1,,,,,,,\n")
	f.close()







def genInput1():
	f= open("input_file.csv","w+")
	f.write("Reactor_Information,Zone1,Zone2,Zone3,,,,,,\n")
	f.write("Zone Length,2,0.25,2,,,,,,\n")
	f.write("Zone Void,0.4,0.4,0.4,,,,,,\n")
	f.write("Reactor Radius,1,,,,,,,,\n")
	f.write("Reactor Temperature,385.65,,,,,,,,\n")
	f.write("Mesh Size,200,,,,,,,,\n")
	f.write("Catalyst Mesh Density,0,,,,,,,,\n")
	f.write("Output Folder Name,new_results,,,,,,,,\n")
	f.write("Experimental Data Folder,none,,,,,,,,\n")
	f.write("Reference Diffusion Inert,16,,,,,,,,\n")
	f.write("Reference Diffusion Catalyst,16,,,,,,,,\n")
	f.write("Reference Temperature,385.5,,,,,,,,\n")
	f.write("Reference Mass,40,,,,,,,,\n")
	f.write("Advection Value,0.0,,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Feed_&_Surface_Composition,,,,,,,,,\n")
	f.write(",CO,,,,,,,,\n")
	f.write("Intensity,1,,,,,,,,\n")
	f.write("Time,0,,,,,,,,\n")
	f.write("Mass,32,,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write(",CO*,*,,,,,,,\n")
	f.write("Initial Concentration,0,30,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Reaction_Information,,,,,,,,,\n")
	f.write("CO + * <-> CO*,2,1,,,,,,,\n")
	f.close()

def genInput2():
	f= open("input_file_2.csv","w+")
	f.write("Reactor_Information,Zone1,Zone2,Zone3,,,,,,\n")
	f.write("Zone Length,2,0.25,2,,,,,,\n")
	f.write("Zone Void,0.4,0.4,0.4,,,,,,\n")
	f.write("Reactor Radius,1,,,,,,,,\n")
	f.write("Reactor Temperature,385.65,,,,,,,,\n")
	f.write("Mesh Size,200,,,,,,,,\n")
	f.write("Catalyst Mesh Density,0,,,,,,,,\n")
	f.write("Output Folder Name,new_results,,,,,,,,\n")
	f.write("Experimental Data Folder,none,,,,,,,,\n")
	f.write("Reference Diffusion Inert,16,,,,,,,,\n")
	f.write("Reference Diffusion Catalyst,16,,,,,,,,\n")
	f.write("Reference Temperature,385.5,,,,,,,,\n")
	f.write("Reference Mass,40,,,,,,,,\n")
	f.write("Advection Value,0.0,,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Feed_&_Surface_Composition,,,,,,,,,\n")
	f.write(",CO,,,,,,,,\n")
	f.write("Intensity,1,,,,,,,,\n")
	f.write("Time,0,,,,,,,,\n")
	f.write("Mass,32,,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write(",CO*,*,,,,,,,\n")
	f.write("Initial Concentration,0,30,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Reaction_Information,,,,,,,,,\n")
	f.write("CO + * <-> CO*,0,0,,,,,,,\n")
	f.close()

def genInput3():
	f= open("input_file_3.csv","w+")
	f.write("Reactor_Information,Zone1,Zone2,Zone3,,,,,,\n")
	f.write("Zone Length,2,0.25,2,,,,,,\n")
	f.write("Zone Void,0.4,0.4,0.4,,,,,,\n")
	f.write("Reactor Radius,1,,,,,,,,\n")
	f.write("Reactor Temperature,385.65,,,,,,,,\n")
	f.write("Mesh Size,200,,,,,,,,\n")
	f.write("Catalyst Mesh Density,0,,,,,,,,\n")
	f.write("Output Folder Name,new_results_3,,,,,,,,\n")
	f.write("Experimental Data Folder,none,,,,,,,,\n")
	f.write("Reference Diffusion Inert,16,,,,,,,,\n")
	f.write("Reference Diffusion Catalyst,16,,,,,,,,\n")
	f.write("Reference Temperature,385.5,,,,,,,,\n")
	f.write("Reference Mass,40,,,,,,,,\n")
	f.write("Advection Value,0.0,,,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Feed_&_Surface_Composition,,,,,,,,,\n")
	f.write(",CO,O2,CO2,Inert-1,,,,,\n")
	f.write("Intensity,1,0.5,0,1,,,,,\n")
	f.write("Time,0,0.1,0,0,,,,,\n")
	f.write("Mass,28,32,44,40,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write(",CO*,O*,*,,,,,,\n")
	f.write("Initial Concentration,0,0,30,,,,,,\n")
	f.write(",,,,,,,,,\n")
	f.write("Reaction_Information,,,,,,,,,\n")
	f.write("CO + * <-> CO*,2,1,,,,,,,\n")
	f.write("O2 + 2* <-> 2O*,2,1,,,,,,,\n")
	f.write("CO* + O* <-> CO2 + 2*,2,1,,,,,,,\n")
	f.close()


def genExample1():
	f= open("tapsolver_example.py","w+")
	f.write("from tapsolver import *\n")
	f.write("\n")
	f.write("# Run the simulation\n")
	f.write("run_tapsolver(timeFunc=0.5,input_file='./input_file.csv',add_noise=True,pulseNumber=10)\n")
	f.write("\n")
	f.write("# Generate an outlet flux graph\n")
	f.write("flux_graph(input_file='./input_file_2.csv',dispExper=True,dispAnalytic=False,pulse=1)\n")
	f.write("\n")
	f.close()

def genExample2():
	f= open("tapsolver_example_2.py","w+")
	f.write("from tapsolver import *\n")
	f.write("\n")
	f.write("# Run the simulation (Sythetic Data)\n")
	f.write("run_tapsolver(timeFunc=0.5,input_file='./input_file.csv',add_noise=True,pulseNumber=10)\n")
	f.write("\n")
	f.write("# Run the optimization (Fitting Parameters)\n")
	f.write("run_tapsolver(timeFunc=0.5,input_file='./input_file.csv',add_noise=True,pulseNumber=10)\n")
	f.write("\n")
	f.close()