
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class initial_conditions():

	def __init__(self):
		self.gasses = {}
		self.inertGasses = {}
		self.surfaceSpecies = {}
		self.advection = 0