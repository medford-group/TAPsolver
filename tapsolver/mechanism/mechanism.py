# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 12:45:02 2023

@author: ayonge
"""
#from .mechanism_constructor import *
#from .mechanism_constructor import *

# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
import sys
import os
from .process import *

import warnings
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

class mechanism():
	processes = None
	def __init__(self):
		
		self._processes = {}
		self._matrix = None
		self._reactions = []
		self._reactants = []
		self._links = {}


	def add_process(self,value):
		new_value = len(self._processes)
		self._processes[new_value] = value
	
	@property
	def processes(self):
		return self._processes
	@processes.setter
	def processes(self,value):
		self._processes[value[0]] = process(value[1])
		mechanism_constructor(self)

	@property
	def matrix(self):
		return self._matrix
	@matrix.setter
	def matrix(self,value):
		self._matrix = value

	@property
	def reactions(self):
		return self._reactions
	@reactions.setter
	def reactions(self,value):
		self._reactions = value

	@property
	def reactants(self):
		return self._reactants
	@reactants.setter
	def reactants(self,value):
		self._reactants = value

	@property
	def links(self):
		return self._links
	@links.setter
	def links(self,value):
		self._links = value

	def write_mech_excel(self,file_name,rewrite=False,use_data=False):
		if os.path.exists(file_name+'.xls') == True:
			if rewrite == True:
				pass
			else:
				print('Mechanism file exists. Change file name or include argument "rewrite = True" to rewrite.')
				print('Not rewriting file and going to next step')
				return

		variable_names = ['reaction','kf','kb','Aof','Aob','Eaf','Eab','dGf','dGb','Gaf','Gab','use_f','use_b','link_f','link_b','b_ub','b_lb','f_ub','f_lb']
		mech_define = np.zeros((len(self.reactions),len(variable_names)))

		df = pd.DataFrame(mech_define,columns=variable_names)
		
		for lnum, l in enumerate(self.reactions):
			df['reaction'][lnum] = l
			df['use_f'][lnum] = 'k'
			df['use_b'][lnum] = 'k'


		if use_data == True:
			for k in range(0,df.shape[0]):
				df['kf'][k] = self.processes[k].f.k
				df['Aof'][k] = self.processes[k].f.Ao
				df['Eaf'][k] = self.processes[k].f.Ea
				df['dGf'][k] = self.processes[k].f.Ga
				df['Gaf'][k] = self.processes[k].f.dG
				df['link_f'][k] = self.processes[k].f.link
				df['use_f'][k] = self.processes[k].f.use
				df['f_lb'][k] = self.processes[k].f.lower_bound
				df['f_ub'][k] = self.processes[k].f.upper_bound

				df['kb'][k] = self.processes[k].b.k  
				df['Aob'][k] = self.processes[k].b.Ao  
				df['Eab'][k] = self.processes[k].b.Ea  
				df['dGb'][k] = self.processes[k].b.Ga  
				df['Gab'][k] = self.processes[k].b.dG 
				df['link_b'][k] = self.processes[k].b.link 
				df['use_b'][k] = self.processes[k].b.use 
				df['b_lb'][k] = self.processes[k].b.lower_bound  
				df['b_ub'][k] = self.processes[k].b.upper_bound 

		filepath = file_name+'.xls'

		df.to_excel(filepath)

		return

	def read_mech_excel(self,file_name):
		df = pd.read_excel(file_name+'.xls',index_col=0)
		for k in range(0,df.shape[0]):
			self.processes[k].reaction = df['reaction'][k]

			self.processes[k].f.k =  df['kf'][k]
			self.processes[k].f.Ao = df['Aof'][k]
			self.processes[k].f.Ea = df['Eaf'][k]
			self.processes[k].f.Ga = df['dGf'][k]
			self.processes[k].f.dG = df['Gaf'][k]
			self.processes[k].f.link = df['link_f'][k]
			self.processes[k].f.use = df['use_f'][k]
			self.processes[k].f.lower_bound = df['f_lb'][k]
			self.processes[k].f.upper_bound = df['f_ub'][k]

			self.processes[k].b.k =  df['kb'][k]
			self.processes[k].b.Ao = df['Aob'][k]
			self.processes[k].b.Ea = df['Eab'][k]
			self.processes[k].b.Ga = df['dGb'][k]
			self.processes[k].b.dG = df['Gab'][k]
			self.processes[k].b.link = df['link_b'][k]
			self.processes[k].b.use = df['use_b'][k]
			self.processes[k].b.lower_bound = df['b_lb'][k]
			self.processes[k].b.upper_bound = df['b_ub'][k]

	def kinetics_conversion(self,file_name,scaling_factor):
		variable_names = ['reaction','kf','kb']
		mech_define = np.zeros((len(self.reactions),len(variable_names)))

		df = pd.DataFrame(mech_define,columns=variable_names)

		for lnum, l in enumerate(self.reactions):
			df['reaction'][lnum] = l
			temp = l.split('<->')
			
			df['kf'][lnum] = self.processes[lnum].f.k*scaling_factor**(temp[0].count('+'))
			df['kb'][lnum] = self.processes[lnum].b.k*scaling_factor**(temp[1].count('+'))

		filepath = file_name+'.xls'

		df.to_excel(filepath)

	def store_matrix(self,file_name):
		df = pd.DataFrame(self.matrix,columns=self.reactants,index=self.reactions)
		filepath = file_name+'.xls'
		df.to_excel(filepath)

def mechanism_constructor(mechanism_data: mechanism):
	
	"""
	
	This function generates the rate_array (or stoichiometric matrix) from the list of elementary reactions provided by the user.
	
	Args:
		
		gasses (dict of ): The

	"""

	insert_location = 0

	def stoich_function(new_var):
		
		if len(new_var) > 2:
			if new_var[:2].isdigit():
				temp_val = val[2:]
				temp_coeff = int(val[:2])
				return temp_val, temp_coeff
			elif new_var[0].isdigit():
				temp_val = val[1:]
				temp_coeff = int(val[0])
				return temp_val, temp_coeff
			else:
				return new_var, 1
		else:
			if new_var[0].isdigit():
				temp_val = val[1:]
				temp_coeff = int(val[0])
				return temp_val, temp_coeff
			else:
				return new_var, 1

	reactants = mechanism_reactants(mechanism_data)

	for k,i in enumerate(reactants):
		if '*' not in i:
			reactants.insert(insert_location, reactants.pop(k))
			insert_location += 1

	rate_array = np.zeros((len(mechanism_data.processes.keys()),len(reactants)))

	reactions = []

	for k,n in enumerate(mechanism_data.processes.keys()):
		reactions.append(mechanism_data.processes[k].reaction)
		mechanism_data.processes[k].reaction = mechanism_data.processes[k].reaction.replace('+','')
		if '<->' in mechanism_data.processes[k].reaction:
			neg,pos = mechanism_data.processes[k].reaction.split('<->')
			new_neg = neg.split()
			new_pos = pos.split()

			for val in new_neg:
				into_array = -1
				new_val = val
				new_val, coeff = stoich_function(new_val)
				rate_array[k,reactants.index(new_val)] = into_array*coeff

			for val in new_pos:
				into_array = 1
				new_val = val
				new_val, coeff = stoich_function(new_val)
				rate_array[k,reactants.index(new_val)] = into_array*coeff

		else:
			neg,pos = mechanism_data.processes[k].reaction.split('->')
			new_neg = neg.split()
			new_pos = pos.split()

			for val in new_neg:
				into_array = -1
				new_val = val
				new_val, coeff = stoich_function(new_val)
				rate_array[k,reactants.index(new_val)] = into_array*coeff

			for val in new_pos:
				into_array = 1
				new_val = val
				new_val, coeff = stoich_function(new_val)
				rate_array[k,reactants.index(new_val)] = into_array*coeff

	mechanism_data.matrix = rate_array
	mechanism_data.reactions = reactions
	mechanism_data.reactants = reactants
	
	return mechanism_data

def mechanism_reactants(mechanism_data: mechanism):

	reactants = []

	for k,i in enumerate(mechanism_data.processes.keys()):
		tempReaction = mechanism_data.processes[k].reaction
		if '<->' in tempReaction:
			tempReaction = tempReaction.replace('<->','')
		else:
			tempReaction = tempReaction.replace('->','')
				
		tempReaction = tempReaction.replace('+','')
		tempReaction = tempReaction.split()
		reactants.append(tempReaction)

	distilledReactants = []
	
	for k in reactants:
		for j in k:
			if len(j) > 2:
				if j[:2].isdigit():
					j = j[2:]
				elif j[0].isdigit():
					j = j[1:]
			if j in distilledReactants:
				pass
			else:
				distilledReactants.append(j)

	insert_location = 0
	active_sites = 0
	
	for k,i in enumerate(distilledReactants):
		#if '*' not in i and '^' not in i and '@' not in i and '#' not in i:
		if '*' not in i:
			distilledReactants.insert(insert_location, distilledReactants.pop(k))
			insert_location += 1

		for k2 in range(0,10):
			if i == '*'+"("+str(k2)+")":
				distilledReactants.pop(k)
				if k2 > active_sites-1:
					active_sites = k2+1
	
	active_site_names = []
	
	for k in range(0,active_sites):
		active_site_names.append('*'+"("+str(k)+")")
	
	for kj in range(0,active_sites):
		distilledReactants.append(active_site_names[kj])

	return distilledReactants