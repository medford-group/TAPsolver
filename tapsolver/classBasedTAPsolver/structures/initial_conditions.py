
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class initial_conditions():

	"""
	
	This class acts as a container for all of the assumed initial conditions of an experiment.
	
	Args:
		
		gasses (dict of ): The  

	"""

	def __init__(self):
		self.gasses = {}
		self.inertGasses = {}
		self.surfaceSpecies = {}
		self.advection = 0