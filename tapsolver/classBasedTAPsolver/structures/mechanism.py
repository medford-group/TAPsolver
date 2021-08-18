
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class mechanism():

	def __init__(self):
		
		self.elementaryProcesses = {}
		self.rateArray = None
		self.rateArray = None
		self.rateArray = None


	def displayProcesses(self):
		print('ELEMENTARY PROCESSES')
		
		for k in self.elementaryProcesses.keys():
			print('__________________________________')
			print('|'+'Reaction '+str(k)+': '+self.elementaryProcesses[k].processString)
			print('|')
			print('|'+'forward - ')
			print('|'+str(self.elementaryProcesses[k].forward.k))
			print('|'+str(self.elementaryProcesses[k].forward.Ao))
			print('|'+str(self.elementaryProcesses[k].forward.Ea))
			print('|'+str(self.elementaryProcesses[k].forward.Ga))
			print('|'+str(self.elementaryProcesses[k].forward.dG))
			print('|'+str(self.elementaryProcesses[k].forward.link))
			print('|')
			print('|'+'backward - ')
			print('|'+str(self.elementaryProcesses[k].backward.k))
			print('|'+str(self.elementaryProcesses[k].backward.Ao))
			print('|'+str(self.elementaryProcesses[k].backward.Ea))
			print('|'+str(self.elementaryProcesses[k].backward.Ga))
			print('|'+str(self.elementaryProcesses[k].backward.dG))
			print('|'+str(self.elementaryProcesses[k].backward.link))
			print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
		print('')
		print('')

	def addProcess(newReaction,):
		
		pass

	def readCSVInput(self, fileName):

		class newProcess():
			def __init__(self):
				self.processString = ''
				self.forward = {}
				self.backward = {}


		class processDetails():
			def __init__(self):
				self.k = {'value':None,'fd':'fixed'}
				self.Ao = {'value':None,'fd':'fixed'}
				self.Ea = {'value':None,'fd':'fixed'}
				self.Ga = {'value':None,'fd':'fixed'}
				self.dG = {'value':None,'fd':'fixed'}
				self.link = {'variable':0}

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

			self.elementaryProcesses[j] = newProcess()
			
			self.elementaryProcesses[j].processString = reaction_info.iloc[j,0]
			self.elementaryProcesses[j].forward = processDetails()
			self.elementaryProcesses[j].backward = processDetails()

			if reaction_info.iloc[j,1].find("#") > 0:
				
				n1, n2 = reaction_info.iloc[j,1].split("#")
				if n1.find("{") == 0:
					self.elementaryProcesses[j].forward.link['Ga'] = float(n1)
				else:
					self.elementaryProcesses[j].forward.Ga = parseValues(n1)
				if n2.find("{") == 0:
					self.elementaryProcesses[j].forward.link['dG'] = float(n2)
				else:
					self.elementaryProcesses[j].forward.dG = parseValues(n2)
			
			elif reaction_info.iloc[j,1].find("$") > 0:

				n1, n2 = reaction_info.iloc[j,1].split("$")
				if n1.find("{") == 0:
					self.elementaryProcesses[j].forward.link['Ao'] = float(n1)
				else:
					self.elementaryProcesses[j].forward.Ao = parseValues(n1)
				if n2.find("{") == 0:
					self.elementaryProcesses[j].forward.link['Ea'] = float(n2)
				else:
					self.elementaryProcesses[j].forward.Ea = parseValues(n2)
			
			else:
				n1 = reaction_info.iloc[j,1]
				if n1.find("{") == 0:
					self.elementaryProcesses[j].forward.Ao = parseValues(n1)
				else:
					self.elementaryProcesses[j].forward.k = parseValues(n1)
			
			if str(reaction_info.iloc[j,2]) != 'nan':
				
				if reaction_info.iloc[j,2].find("$") > 0:
					n1, n2 = reaction_info.iloc[j,2].split("$")
					if n1.find("{") == 0:
						self.elementaryProcesses[j].forward.link['Ao'] = float(n1)
					else:
						self.elementaryProcesses[j].backward.Ao = parseValues(n1)
					if n2.find("{") == 0:
						self.elementaryProcesses[j].forward.link['Ea'] = float(n1)
					else:	
						self.elementaryProcesses[j].backward.Ea = parseValues(n2)
				
				else:
					n1 = reaction_info.iloc[j,2]
					if n1.find("{") == 0:
						self.elementaryProcesses[j].forward.backward['k'] = float(n1)
					else:
						self.elementaryProcesses[j].backward.k = parseValues(n1)


	def meachanismReactants(self):

		reactants = []

		for k,i in enumerate(self.elementaryProcesses.keys()):
			tempReaction = self.elementaryProcesses[k].processString
			print(tempReaction)
			if '<->' in tempReaction:
				tempReaction = tempReaction.replace('<->','')
			else:
				tempReaction = tempReaction.replace('->','')
				
			tempReaction = tempReaction.replace('+','')
			tempReaction = tempReaction.split()
			reactants.append(tempReaction)

		distilledReactants = []
		print(reactants)
		for k in reactants:

			for j in k:
				if j[0].isdigit():
					j = j[1:]
				if j in distilledReactants:
					pass
				else:
					distilledReactants.append(j)
		print(distilledReactants)

	def reactionMatrix(self,initCond):
		


		insert_location = 0
		active_sites = 0


		for k,i in enumerate(reactants):
			if '*' not in i:
				reactants.insert(insert_location, reactants.pop(k))
				insert_location += 1
	
		rate_array = np.zeros((len(self.elementaryProcesses.keys()),len(initCond.gasses.keys())+len(initCond.surface_species.keys())))
	
		for k,n in enumerate(reactions):
			reactions[k] = reactions[k].replace('+','')
			if '<->' in reactions[k]:
				neg,pos = reactions[k].split('<->')
				new_neg = neg.split()
				new_pos = pos.split()
	
				for val in new_neg:
					into_array = -1
					new_val = val
					if val[0].isdigit():
						new_val = val[1:]
						into_array = into_array*int(val[0])
					rate_array[k,reactants.index(new_val)] = into_array
	
				for val in new_pos:
					into_array = 1
					new_val = val
					if val[0].isdigit():
						new_val = val[1:]
						into_array = into_array*int(val[0])
					rate_array[k,reactants.index(new_val)] = into_array
			else:
				neg,pos = reactions[k].split('->')
				new_neg = neg.split()
				new_pos = pos.split()
	
				for val in new_neg:
					into_array = -1
					new_val = val
					if val[0].isdigit():
						new_val = val[1:]
						into_array = into_array*int(val[0])
					rate_array[k,reactants.index(new_val)] = into_array
	
				for val in new_pos:
					into_array = 1
					new_val = val
					if val[0].isdigit():
						new_val = val[1:]
						into_array = into_array*int(val[0])
					rate_array[k,reactants.index(new_val)] = into_array
	
		return rate_array, reactions, reactants
