
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from .reactor_species import reactor_species

def display_gasses(reactor_species_data: reactor_species):

	print('REACTION/PRODUCT GASSES')
	if len(reactor_species_data.gasses) == 0:
		print('_________________')
		print('| None')
		print('~~~~~~~~~~~~~~~~~')
	else:
		for k in reactor_species_data.gasses.keys():
			print('_________________')
			print('|'+k)
			print('|General_Parameters')
			print('|-'+'mass: '+str(reactor_species_data.gasses[k].mass))
			print('|TAP_Parameters')
			print('|-'+'Intensity: '+str(reactor_species_data.gasses[k].intensity))
			print('|-'+'Delay: '+str(reactor_species_data.gasses[k].delay))
			print('|-'+'Inert Diffusion: '+str(reactor_species_data.gasses[k].inert_diffusion))
			print('|-'+'Catalyst Diffusion: '+str(reactor_species_data.gasses[k].catalyst_diffusion))
			print('|Alternative_Reactor_Parameters')
			print('|-'+'Initial Concentration: '+str(reactor_species_data.gasses[k].initial_concentration))
			print('|-'+'Inlet Concentration: '+str(reactor_species_data.gasses[k].inlet_concentration))
			print('~~~~~~~~~~~~~~~~~')

	print('')
	print('')
	print('INERT GASSES')

	if len(reactor_species_data.inert_gasses) == 0:
		print('_________________')
		print('| None')
		print('~~~~~~~~~~~~~~~~~')

	for k in reactor_species_data.inert_gasses.keys():
			print('_________________')
			print('|'+k)
			print('|General_Parameters')
			print('|-'+'mass: '+str(reactor_species_data.inert_gasses[k].mass))
			print('|TAP_Parameters____')
			print('|-'+'Intensity: '+str(reactor_species_data.inert_gasses[k].intensity))
			print('|-'+'Delay: '+str(reactor_species_data.inert_gasses[k].delay))
			print('|-'+'Inert Diffusion: '+str(reactor_species_data.inert_gasses[k].inert_diffusion))
			print('|-'+'Catalyst Diffusion: '+str(reactor_species_data.inert_gasses[k].catalyst_diffusion))
			print('|Alternative_Reactor_Parameters')
			print('|-'+'Initial Concentration: '+str(reactor_species_data.inert_gasses[k].initial_concentration))
			print('|-'+'Inlet Concentration: '+str(reactor_species_data.inert_gasses[k].inlet_concentration))
			print('~~~~~~~~~~~~~~~~~')


	print('')
	print('')