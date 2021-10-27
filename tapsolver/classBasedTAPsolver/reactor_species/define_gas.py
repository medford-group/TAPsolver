
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
import math as mp
from structures import reactor_species

class define_gas():
	def __init__(self, mass = 0):
		self.inert_diffusion = 0
		self.catalyst_diffusion = 0
		self.intensity = 0
		self.delay = 0

		self.initial_concentration = 0
		self.inlet_concentration = 0
		self.mass = mass

	@property
	def mass(self):
		return self._mass
	@mass.setter
	def mass(self,value):
		
		if value < 0:
			raise ValueError("Mass must be positive (non-negative)")
		
		if self.inert_diffusion > 0 and self.catalyst_diffusion > 0:
			prior_mass = self.mass
			self._mass = value
			self.inert_diffusion = self.inert_diffusion*mp.sqrt(prior_mass)/mp.sqrt(self._mass)
			self.catalyst_diffusion = self.catalyst_diffusion*mp.sqrt(prior_mass)/mp.sqrt(self._mass)
		else:
			self._mass = value