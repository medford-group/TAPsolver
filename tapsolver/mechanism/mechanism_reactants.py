
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from .mechanism import mechanism
import sys

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
				if j[:1].isdigit():
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