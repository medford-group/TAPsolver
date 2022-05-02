
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from mechanism import mechanism
import sys

def mechanism_reactants(mechanism_data: mechanism):

	reactants = []

	for k,i in enumerate(mechanism_data.elementary_processes.keys()):
		tempReaction = mechanism_data.elementary_processes[k].processString
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
			if j[0].isdigit():
				j = j[1:]
			if j in distilledReactants:
				pass
			else:
				distilledReactants.append(j)

	insert_location = 0
	active_sites = 0

	for k,i in enumerate(distilledReactants):
		if '*' not in i and '^' not in i and '@' not in i and '#' not in i:
			distilledReactants.insert(insert_location, distilledReactants.pop(k))
			insert_location += 1

		if i == '*':
			distilledReactants.pop(k)
			if active_sites == 0:
				active_sites += 1
		if i == '^':
			distilledReactants.pop(k)
			if active_sites == 1:
				active_sites += 1
		if i == '@':
			distilledReactants.pop(k)
			if active_sites == 2:
				active_sites += 1
		if i == '#':
			distilledReactants.pop(k)
			if active_sites == 3:
				active_sites += 1

	active_site_names = ['*','^','@','#']
	
	for kj in range(0,int(active_sites)):
		distilledReactants.insert(len(distilledReactants), active_site_names[kj])

	return distilledReactants