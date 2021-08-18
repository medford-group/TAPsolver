
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import initial_conditions

def readCSV_ic(initial_conditions_data: initial_conditions, fileName):

	data = pd.read_csv(fileName,header=None)
	rows_1, cols_1 = np.where(data == 'Reactor_Information')
	rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
	rows_4, cols_4 = np.where(data == 'Reaction_Information')

	feed_surf_info = data.iloc[1+rows_2[0]:rows_4[0]-1,:]

	gas_species = [x for x in feed_surf_info.iloc[0,1:] if type(x) == str]
	surface_species = [x for x in feed_surf_info.iloc[5,1:] if type(x) == str]

	for jnum,j in enumerate(gas_species):


		class specificGas():
			def __init__(initial_conditions_data):
				initial_conditions_data.mass = 0
				initial_conditions_data.intensity = 0
				initial_conditions_data.delay = 0


		class specificAdspecies():
			def __init__(initial_conditions_data):
				initial_conditions_data.concentration = 0

		if j.find('Inert') == 0:
			initial_conditions_data.inertGasses[j] = specificGas()
			initial_conditions_data.inertgasses[j].intensity = float(feed_surf_info.iloc[1,1+jnum])
			initial_conditions_data.inertgasses[j].delay = float(feed_surf_info.iloc[2,1+jnum])
			initial_conditions_data.inertgasses[j].mass = float(feed_surf_info.iloc[3,1+jnum])

		else:
			initial_conditions_data.gasses[j] = specificGas()
			initial_conditions_data.gasses[j].intensity = float(feed_surf_info.iloc[1,1+jnum])
			initial_conditions_data.gasses[j].delay = float(feed_surf_info.iloc[2,1+jnum])
			initial_conditions_data.gasses[j].mass = float(feed_surf_info.iloc[3,1+jnum])


	for jnum,j in enumerate(surface_species):
		initial_conditions_data.surfaceSpecies[j] = specificAdspecies()
		if jnum == len(surface_species)-1:
			initial_conditions_data.surfaceSpecies[j].conectration = feed_surf_info.iloc[6,1+jnum]