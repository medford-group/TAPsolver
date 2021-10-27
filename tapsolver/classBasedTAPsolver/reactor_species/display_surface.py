
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import reactor_species


def display_surface(reactor_species_data: reactor_species):
	print('Surface species/adspecies')
	for k in reactor_species_data.adspecies.keys():
		print('_________________')
		print('|'+k)
		print('|Conecentration: '+str(reactor_species_data.adspecies[k].concentration))
		print('~~~~~~~~~~~~~~~~~')
	if len(reactor_species_data.adspecies) == 0:
		print('_________________')
		print('| None')
		print('~~~~~~~~~~~~~~~~~')
	print('')
	print('')