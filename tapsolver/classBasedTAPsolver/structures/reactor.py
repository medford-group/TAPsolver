
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class reactor():
	
	"""
	
	This class acts as a container for all of the reactor parameters
	
	Args:
		**TAPsap Variables**

		zone_lengths (dict): The zone lengths for the reactor.

		zone_void (dict): The bed porosity within each zone.

		zone_diffusion (dict): The diffusion coefficient for each zone.

		zone_residence_time (dict): The residence time within each zone.

		reactor_radius (float): The radius of the reactor.

		catalyst_weight (float): The catalyst weight.

		mol_per_pulse (float): The mol per pulse.

		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		**TAPsolver Specific Variables**

		

	"""

	def __init__(self):
		#self.length = [2.80718,0.17364,2.80718]
		
		#self.void = [0.4,0.4,0.4]
		
		#self.radius = 1.0
		
		#self.temp = 385.65
		
		#self.mesh = 200
		
		self.meshDensity = 4
		
		self.outputName = 'exp_new'
		
		self.inertDiff = 16
		
		self.catalystDiff = 16
		
		self.refTemp = 385.6
		
		self.refMass = 40
		
		self.advection = 0

	def readCSVInput(self, fileName):

		data = pd.read_csv(fileName,header=None)

		rows_1, cols_1 = np.where(data == 'Reactor_Information')
		rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(data == 'Reaction_Information')

		reactor_info = data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 

		self.length = [float(reactor_info.iloc[0,1]), float(reactor_info.iloc[0,2]), float(reactor_info.iloc[0,3])]
		self.void = [float(reactor_info.iloc[1,1]), float(reactor_info.iloc[1,2]), float(reactor_info.iloc[1,3])]
		self.radius = float(reactor_info.iloc[2,1])
		self.temp = float(reactor_info.iloc[3,1])
		self.mesh = int(reactor_info.iloc[4,1])
		self.meshDensity = int(reactor_info.iloc[5,1])
		self.outputName = reactor_info.iloc[6,1]
		self.inertDiff = float(reactor_info.iloc[8,1])
		self.catalystDiff = float(reactor_info.iloc[9,1])
		self.refTemp = float(reactor_info.iloc[10,1])
		self.refMass = float(reactor_info.iloc[11,1])
		self.advection = float(reactor_info.iloc[12,1])

	def lengthFractions(self):

		return [self.length[0]/sum(self.length),self.length[1]/sum(self.length),self.length[2]/sum(self.length)]

	def reactorCenterFraction(self):
		
		return (self.length[0] + self.length[0]/2)/sum(self.length)