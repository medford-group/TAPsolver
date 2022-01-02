
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import initial_conditions


def display_surface(initial_conditions_data: initial_conditions):
	print('REACTION/PRODUCT GASSES')
	for k in initial_conditions_data.surfaceSpecies.keys():
		print('_________________')
		print('|'+k)
		print('|Conecentration: '+str(initial_conditions_data.surfaceSpecies[k].concentration))
		print('~~~~~~~~~~~~~~~~~')
	if len(initial_conditions_data.surfaceSpecies) == 0:
		print('_________________')
		print('| None')
		print('~~~~~~~~~~~~~~~~~')
	print('')
	print('')