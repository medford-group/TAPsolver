
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import reactor_species
from reactor_species import define_gas, define_adspecies
def readCSV_reactor_species(fileName):

	"""
	This function converts the input file initial condition information into an initial condition object.

	Args:
		input_file (str): Directory path and name of the input file.

	Returns:
		An initial condition object with reactor parameters from the input_file.

	"""

	reactor_species_data = reactor_species()

	data = pd.read_csv(fileName,header=None,dtype=object)
	
	rows_1, cols_1 = np.where(data == 'Reactor_Information')
	rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
	rows_4, cols_4 = np.where(data == 'Reaction_Information')

	reactor_info = data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = data.iloc[1+rows_2[0]:rows_4[0]-1,:]

	gas_species = [x for x in feed_surf_info.iloc[0,1:] if type(x) == str]
	surface_species = [x for x in feed_surf_info.iloc[5,1:] if type(x) == str]

	for jnum,j in enumerate(gas_species):

		if j.find('Inert') == 0:
			reactor_species_data.inertGasses[j] = define_gas()
			reactor_species_data.inertgasses[j].intensity = float(feed_surf_info.iloc[1,1+jnum])
			reactor_species_data.inertgasses[j].delay = float(feed_surf_info.iloc[2,1+jnum])
			reactor_species_data.inertgasses[j].mass = float(feed_surf_info.iloc[3,1+jnum])

		else:
			reactor_species_data.gasses[j] = define_gas()
			reactor_species_data.gasses[j].intensity = float(feed_surf_info.iloc[1,1+jnum])
			reactor_species_data.gasses[j].delay = float(feed_surf_info.iloc[2,1+jnum])
			reactor_species_data.gasses[j].mass = float(feed_surf_info.iloc[3,1+jnum])


	reactor_species_data.inert_diffusion = float(reactor_info.iloc[8,1])
	reactor_species_data.catalyst_diffusion = float(reactor_info.iloc[9,1])
	reactor_species_data.reference_temperature = float(reactor_info.iloc[10,1])
	reactor_species_data.reference_mass = float(reactor_info.iloc[11,1])
	reactor_species_data.advection = float(reactor_info.iloc[12,1])

	for jnum,j in enumerate(surface_species):
		reactor_species_data.surfaceSpecies[j] = define_adspecies()
		if jnum == len(surface_species)-1:
			reactor_species_data.surfaceSpecies[j].conectration = feed_surf_info.iloc[6,1+jnum]

	return reactor_species_data