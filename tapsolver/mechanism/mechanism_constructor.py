# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from .mechanism import mechanism
from .mechanism_reactants import mechanism_reactants
import sys

def mechanism_constructor(mechanism_data: mechanism):
	
	"""
	
	This function generates the rate_array (or stoichiometric matrix) from the list of elementary reactions provided by the user.
	
	Args:
		
		gasses (dict of ): The

	"""

	insert_location = 0

	def stoich_function(new_var):
		
		if len(new_var) > 2:
			if new_var[:1].isdigit():
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