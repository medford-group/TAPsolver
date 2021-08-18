
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import mechanism

def readCSV_mechanism(mechanism_data: mechanism, fileName):

	class newProcess():
		def __init__(mechanism_data):
			mechanism_data.processString = ''
			mechanism_data.forward = {}
			mechanism_data.backward = {}


	class processDetails():
		def __init__(mechanism_data):
			mechanism_data.k = {'value':None,'fd':'fixed'}
			mechanism_data.Ao = {'value':None,'fd':'fixed'}
			mechanism_data.Ea = {'value':None,'fd':'fixed'}
			mechanism_data.Ga = {'value':None,'fd':'fixed'}
			mechanism_data.dG = {'value':None,'fd':'fixed'}
			mechanism_data.link = {'variable':0}

	def parseValues(v1):
		
		v1Dict = {}
		if v1.find("!") < 0:
			v1Dict['value'] = float(v1)
			v1Dict['fd'] = 'dynamic'
		else:
			v1Dict['value'] = float(v1[:-1])
			v1Dict['fd'] = 'fixed'

		return v1Dict

	data = pd.read_csv(fileName,header=None)

	rows_1, cols_1 = np.where(data == 'Reactor_Information')
	rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
	rows_4, cols_4 = np.where(data == 'Reaction_Information')

	if data[0].str.contains('Linked Kinetics').any():
		rows_5, cols_5 = np.where(data == 'Linked Kinetics')
	
		if data[0].str.contains('Thermodynamic Constraints').any():
			rows_6, cols_6 = np.where(data == 'Thermodynamic Constraints')
			reaction_info = data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			linked_kinetics = data.iloc[1+rows_5[0]:(rows_6[0]-1),]
			thermo_constraints = data.iloc[1+rows_6[0]:,:]
	
		else:
			reaction_info = data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			linked_kinetics = data.iloc[1+rows_5[0]:,:]
			thermo_constraints = None

	elif data[0].str.contains('Thermodynamic Constraints').any():
		rows_6, cols_6 = np.where(data == 'Thermodynamic Constraints')   
		reaction_info = data.iloc[(1+rows_4[0]):(rows_6[0]-1),:]
		linked_kinetics = None
		thermo_constraints = data.iloc[1+rows_6[0]:,:]
			
	else:
		reaction_info = data.iloc[(1+rows_4[0]):,:]
		linked_kinetics = None
		thermo_constraints = None

	for j in range(0,len(reaction_info.index)):

		mechanism_data.elementaryProcesses[j] = newProcess()
			
		mechanism_data.elementaryProcesses[j].processString = reaction_info.iloc[j,0]
		mechanism_data.elementaryProcesses[j].forward = processDetails()
		mechanism_data.elementaryProcesses[j].backward = processDetails()

		if reaction_info.iloc[j,1].find("#") > 0:
			
			n1, n2 = reaction_info.iloc[j,1].split("#")
			if n1.find("{") == 0:
				mechanism_data.elementaryProcesses[j].forward.link['Ga'] = float(n1)
			else:
				mechanism_data.elementaryProcesses[j].forward.Ga = parseValues(n1)
			if n2.find("{") == 0:
				mechanism_data.elementaryProcesses[j].forward.link['dG'] = float(n2)
			else:
				mechanism_data.elementaryProcesses[j].forward.dG = parseValues(n2)
			
		elif reaction_info.iloc[j,1].find("$") > 0:

			n1, n2 = reaction_info.iloc[j,1].split("$")
			if n1.find("{") == 0:
				mechanism_data.elementaryProcesses[j].forward.link['Ao'] = float(n1)
			else:
				mechanism_data.elementaryProcesses[j].forward.Ao = parseValues(n1)
			if n2.find("{") == 0:
				mechanism_data.elementaryProcesses[j].forward.link['Ea'] = float(n2)
			else:
				mechanism_data.elementaryProcesses[j].forward.Ea = parseValues(n2)
			
		else:
			n1 = reaction_info.iloc[j,1]
			if n1.find("{") == 0:
				mechanism_data.elementaryProcesses[j].forward.Ao = parseValues(n1)
			else:
				mechanism_data.elementaryProcesses[j].forward.k = parseValues(n1)
			
		if str(reaction_info.iloc[j,2]) != 'nan':
				
			if reaction_info.iloc[j,2].find("$") > 0:
				n1, n2 = reaction_info.iloc[j,2].split("$")
				if n1.find("{") == 0:
					mechanism_data.elementaryProcesses[j].forward.link['Ao'] = float(n1)
				else:
					mechanism_data.elementaryProcesses[j].backward.Ao = parseValues(n1)
				if n2.find("{") == 0:
					mechanism_data.elementaryProcesses[j].forward.link['Ea'] = float(n1)
				else:	
					mechanism_data.elementaryProcesses[j].backward.Ea = parseValues(n2)
			
			else:
				n1 = reaction_info.iloc[j,2]
				if n1.find("{") == 0:
					mechanism_data.elementaryProcesses[j].forward.backward['k'] = float(n1)
				else:
					mechanism_data.elementaryProcesses[j].backward.k = parseValues(n1)
