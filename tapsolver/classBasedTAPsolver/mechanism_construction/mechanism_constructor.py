
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import mechanism
from .mechanism_reactants import mechanism_reactants

def mechanism_constructor(mechanism_data: mechanism):
	
	insert_location = 0
	active_sites = 0

	reactants = mechanism_reactants(mechanism_data)

	for k,i in enumerate(reactants):
		if '*' not in i:
			reactants.insert(insert_location, reactants.pop(k))
			insert_location += 1

	rate_array = np.zeros((len(mechanism_data.elementaryProcesses.keys()),len(reactants)))

	reactions = []

	for k,n in enumerate(mechanism_data.elementaryProcesses.keys()):
		reactions.append(mechanism_data.elementaryProcesses[k].processString)
		mechanism_data.elementaryProcesses[k].processString = mechanism_data.elementaryProcesses[k].processString.replace('+','')
		if '<->' in mechanism_data.elementaryProcesses[k].processString:
			neg,pos = mechanism_data.elementaryProcesses[k].processString.split('<->')
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
			neg,pos = mechanism_data.elementaryProcesses[k].processString.split('->')
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

	mechanism_data.rate_array = rate_array
	mechanism_data.reactions = reactions
	mechanism_data.reactants = reactants

	return mechanism_data