
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
import sys
from structures import reactor

def readCSV_reactor(fileName):

	"""
	This function converts the input file reactor information into a reactor object.

	Args:
		fileName (str): Directory path and name of the input file.

	Returns:
		A reactor object with reactor parameters from the fileName.

	"""

	reactor_data = reactor()

	data = pd.read_csv(fileName,header=None,dtype=object)

	rows_1, cols_1 = np.where(data == 'Reactor_Information')
	rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
	rows_4, cols_4 = np.where(data == 'Reaction_Information')
	
	reactor_info = data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 

	reactor_data.zone_lengths = {'zone0': float(reactor_info.iloc[0,1]), 'zone1':  float(reactor_info.iloc[0,2]), 'zone2': float(reactor_info.iloc[0,3])}
	reactor_data.zone_voids = {'zone0': float(reactor_info.iloc[1,1]), 'zone1': float(reactor_info.iloc[1,2]), 'zone2': float(reactor_info.iloc[1,3])}
	reactor_data.reactor_radius = float(reactor_info.iloc[2,1])
	reactor_data.temperature = float(reactor_info.iloc[3,1])
	reactor_data.mesh = int(reactor_info.iloc[4,1])
	reactor_data.catalyst_mesh_density = int(reactor_info.iloc[5,1])
	reactor_data.output_name = reactor_info.iloc[6,1]
	reactor_data.inert_diffusion = float(reactor_info.iloc[8,1])
	reactor_data.catalyst_diffusion = float(reactor_info.iloc[9,1])
	reactor_data.reference_temperature = float(reactor_info.iloc[10,1])
	reactor_data.reference_mass = float(reactor_info.iloc[11,1])
	reactor_data.advection = float(reactor_info.iloc[12,1])

	return reactor_data