
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import reactor

def readCSV_reactor(reactor_data: reactor, fileName):

	data = pd.read_csv(fileName,header=None)

	rows_1, cols_1 = np.where(data == 'Reactor_Information')
	rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
	rows_4, cols_4 = np.where(data == 'Reaction_Information')

	reactor_info = data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 

	reactor_data.length = [float(reactor_info.iloc[0,1]), float(reactor_info.iloc[0,2]), float(reactor_info.iloc[0,3])]
	reactor_data.void = [float(reactor_info.iloc[1,1]), float(reactor_info.iloc[1,2]), float(reactor_info.iloc[1,3])]
	reactor_data.radius = float(reactor_info.iloc[2,1])
	reactor_data.temp = float(reactor_info.iloc[3,1])
	reactor_data.mesh = int(reactor_info.iloc[4,1])
	reactor_data.meshDensity = int(reactor_info.iloc[5,1])
	reactor_data.outputName = reactor_info.iloc[6,1]
	reactor_data.inertDiff = float(reactor_info.iloc[8,1])
	reactor_data.catalystDiff = float(reactor_info.iloc[9,1])
	reactor_data.refTemp = float(reactor_info.iloc[10,1])
	reactor_data.refMass = float(reactor_info.iloc[11,1])
	reactor_data.advection = float(reactor_info.iloc[12,1])