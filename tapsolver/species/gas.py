
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
import math as mp
import copy
#from structures import reactor_species

class gas():
	def __init__(self,new_mass=40):
		self._inert_diffusion = 0
		self._catalyst_diffusion = 0
		self._intensity = 0
		self._delay = 0
		self._noise = 0 
		self._sigma = 0.1
		self._temperature_used = 0

		self._in_conc = 0
		self._mass = new_mass

	@property
	def mass(self):
		return self._mass
	@mass.setter
	def mass(self,value):
		
		if (type(value) == float) or (type(value) == int):

			if value < 0:
				raise ValueError("Mass must be positive ")
			
			prior_mass = copy.deepcopy(self._mass)
			self._mass = value
			self.inert_diffusion = self.inert_diffusion*mp.sqrt(prior_mass)/mp.sqrt(self.mass)
			self.catalyst_diffusion = self.catalyst_diffusion*mp.sqrt(prior_mass)/mp.sqrt(self.mass)
			#print(self._inert_diffusion)
			#print(self._catalyst_diffusion)
		else:
			self._mass = value

	#@property
	#def delay(self):
	#	return self._delay
	#@delay.setter
	#def delay(self,a):
	#	self._delay = a

	@property
	def delay(self):
		return self._delay
	@delay.setter
	def delay(self,a):
		self._delay = a

	@property
	def intensity(self):
		return self._intensity
	@intensity.setter
	def intensity(self,a):
		self._intensity = a

	@property
	def noise(self):
		return self._noise
	@noise.setter
	def noise(self,a):
		self._noise = a

	@property
	def sigma(self):
		return self._sigma
	@sigma.setter
	def sigma(self,a):
		self._sigma = a

	@property
	def temperature_used(self):
		return self._temperature_used
	@temperature_used.setter
	def temperature_used(self,a):
		self._temperature_used = a
