
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import initial_conditions

def display_gasses(initial_conditions_data: initial_conditions):

	print('REACTION/PRODUCT GASSES')

	for k in initial_conditions_data.gasses.keys():
		print('_________________')
		print('|'+k)
		print('|'+'Intensity: '+str(initial_conditions_data.gasses[k].intensity))
		print('|'+'Delay: '+str(initial_conditions_data.gasses[k].delay))
		print('|'+'mass: '+str(initial_conditions_data.gasses[k].mass))#['intensity']
		print('~~~~~~~~~~~~~~~~~')

	if len(initial_conditions_data.gasses) == 0:
		print('_________________')
		print('| None')
		print('~~~~~~~~~~~~~~~~~')

	print('')
	print('')
	print('INERT GASSES')

	for k in initial_conditions_data.inertGasses.keys():
		print('_________________')
		print('|'+k)
		print('|Intensity: '+str(initial_conditions_data.inertGasses[k].intensity))
		print('|Delay: '+str(initial_conditions_data.inertGasses[k].delay))
		print('|mass: '+str(initial_conditions_data.inertGasses[k].mass))#['intensity']
		print('~~~~~~~~~~~~~~~~~')

	if len(initial_conditions_data.inertGasses) == 0:
		print('_________________')
		print('| None')
		print('~~~~~~~~~~~~~~~~~')

	print('')
	print('')