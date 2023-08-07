
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved


import pandas as pd
import numpy as np
import math as mp
from .gas import gas
from .adspecies import adspecies

class species():

	"""
	
	This class acts as a container for all of the assumed initial conditions of an experiment.
	
	Args:
		
		gasses (dict of ): The  

	"""

	def __init__(self):#, inert_diffusion = 16, catalyst_diffusion = 16, reference_temperature = 385.6, reference_mass = 40, temperature = 385.6):
		self.gasses = {}
		self.inert_gasses = {}
		self.adspecies = {}
		self.inert_diffusion = 16
		self.catalyst_diffusion = 16
		#self.reference_diffusions = {0:16,1:16,2:16} 
		self.reference_temperature = 385.6
		self.reference_mass = 40
		self.temperature = 385.6
		self.advection = 0
		self.reference_pulse_size = 1

	def initialize(self,b):
		for a in b:
			if '*' in a:
				s = adspecies()
				s.conc = 0
				self.add_adspecies(a,s)

			if '*' not in a:
				s = gas()
				s.in_conc = 0
				self.add_gas(a,s)

	def add_gas(self,name='', gas_data = gas):
		def calculate_diffusion_coefficient(rd, rm, rt, temp, ms):
			return rd*(mp.sqrt(rm*temp)/mp.sqrt(rt*ms))

		if name not in self.gasses:
			self.gasses[name] = gas_data
			try:
				self.gasses[name].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[name].mass)
				self.gasses[name].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[name].mass)
			except:
				self.gasses[name].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self._reference_mass, self._reference_temperature, self._temperature[self._gasses[name].temperature_used], self._gasses[name].mass)
				self.gasses[name].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self._reference_mass, self._reference_temperature, self._temperature[self._gasses[name].temperature_used], self._gasses[name].mass)
		else:
			pass
			#print('Gas already defined in dictionary.')

	def add_inert_gas(self,name='', gas_data = gas):
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))

		if name not in self.inert_gasses:
			self.inert_gasses[name] = gas_data
			try:
				pass
				#self.inert_gasses[name].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self._reference_mass, self._reference_temperature, self._temperature, self._inert_gasses[name].mass)
				#self.inert_gasses[name].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self._reference_mass, self._reference_temperature, self._temperature, self._inert_gasses[name].mass)
			except:
				pass
				#self.inert_gasses[name].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self._reference_mass, self._reference_temperature, self._temperature[self.inert_gasses[name].temperature_used], self.inert_gasses[name].mass)
				#self.inert_gasses[name].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self._reference_mass, self._reference_temperature, self._temperature[self.inert_gasses[name].temperature_used], self.inert_gasses[name].mass)
			
		else:
			print('Gas already defined in dictionary.')

	def add_adspecies(self,name='', adspecies_data = adspecies):
		
		if name not in self.adspecies:
			self.adspecies[name] = adspecies_data
		else:
			print('Gas already defined in dictionary.')

	@property
	def inert_diffusion(self):
		return self._inert_diffusion
	@inert_diffusion.setter
	def inert_diffusion(self,value):
		
		if (type(value) == float) or (type(value) == int):
			
			def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
				return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

			if value <= 0:
				raise ValueError("Reference inert diffusion must be positive (non-negative")
			self._inert_diffusion = value
			for j in self.gasses:
				try:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
				except:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
			for j in self.inert_gasses:
				try:
					self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self._inert_gasses[j].mass)
				except:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self._inert_diffusion, self.reference_mass, self.reference_temperature, self.temperature[self._inert_gasses[j].temperature_used], self._inert_gasses[j].mass)
			
		else:
			self._inert_diffusion = value

	@property
	def catalyst_diffusion(self):
		return self._catalyst_diffusion
	@catalyst_diffusion.setter
	def catalyst_diffusion(self,value):		
		if (type(value) == float) or (type(value) == int):
			def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
				return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

			if value <= 0:
				raise ValueError("Reference catalyst diffusion must be positive (non-negative")
			self._catalyst_diffusion = value
			
			for j in self.gasses:
				try:
					self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
				except:
					self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
			for j in self.inert_gasses:
				try:
					self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature, self._inert_gasses[j].mass)
				except:
					self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self._catalyst_diffusion, self.reference_mass, self.reference_temperature, self.temperature[self._inert_gasses[j].temperature_used], self._inert_gasses[j].mass)
					
		else:
			self._catalyst_diffusion = value

	@property
	def reference_temperature(self):
		return self._reference_temperature
	@reference_temperature.setter
	def reference_temperature(self,value):
		if (type(value) == float) or (type(value) == int):

			def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
				return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	

			#if type(value) != Constant: 
			if value <= 0:
				raise ValueError("Temperature must be positive (non-negative)")
			self._reference_temperature = value
			for j in self.gasses:
				try:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self._reference_temperature, self.temperature, self.gasses[j].mass)
					self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self._reference_temperature, self.temperature, self.gasses[j].mass)
				except:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self._reference_temperature, self.temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self._reference_temperature, self.temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
			
			for j in self.inert_gasses:
				try:
					self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self._reference_temperature, self.temperature, self.inert_gasses[j].mass)
					self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self._reference_temperature, self.temperature, self.inert_gasses[j].mass)
				except:
					self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self._reference_temperature, self.temperature[self.gasses[j].temperature_used], self.inert_gasses[j].mass)
					self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self._reference_temperature, self.temperature[self._inert_gasses[j].temperature_used], self.inert_gasses[j].mass)
					
		else:
			self._reference_temperature = value

	@property
	def reference_mass(self):
		return self._reference_mass
	@reference_mass.setter
	def reference_mass(self,value):		
		if (type(value) == float) or (type(value) == int) or (type(value) == list):
			def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
				return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	
			if value <= 0:
				raise ValueError("Species mass must be positive (non-negative")
			self._reference_mass = value
			for j in self.gasses:
				try:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self._reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
					self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self._reference_mass, self.reference_temperature, self.temperature, self.gasses[j].mass)
				except:
					self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self._reference_mass, self.reference_temperature, self.temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
					self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self._reference_mass, self.reference_temperature, self.temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
			for j in self.inert_gasses:
				try:
					self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self._reference_mass, self.reference_temperature, self.temperature, self.inert_gasses[j].mass)
					self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self._reference_mass, self.reference_temperature, self.temperature, self.inert_gasses[j].mass)
				except:
					self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self._reference_mass, self.reference_temperature, self.temperature[self.inert_gasses[j].temperature_used], self.inert_gasses[j].mass)
					self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self._reference_mass, self.reference_temperature, self.temperature[self.inert_gasses[j].temperature_used], self.inert_gasses[j].mass)
					
		else:
			self._reference_mass = value

	@property
	def temperature(self):
		return self._temperature
	@temperature.setter
	def temperature(self,value):
		def calculate_diffusion_coefficient(reference_diffusion, reference_mass, reference_temperature, temperature, mass):
			return reference_diffusion*(mp.sqrt(reference_mass*temperature)/mp.sqrt(reference_temperature*mass))	
		
		self._temperature = value
		
		for j in self.gasses:
			
			if type(self._temperature) != dict:
				self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self._temperature, self.gasses[j].mass)
				self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self._temperature, self.gasses[j].mass)
			
			elif type(self._temperature) != {}:
				self.gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self._temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
				self.gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self._temperature[self.gasses[j].temperature_used], self.gasses[j].mass)
		
		for j in self.inert_gasses:
			
			if type(self._temperature) != dict:
				self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self._temperature, self.inert_gasses[j].mass)
				self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self._temperature, self.inert_gasses[j].mass)
			
			elif type(self._temperature) != {}:
				self.inert_gasses[j].inert_diffusion = calculate_diffusion_coefficient(self.inert_diffusion, self.reference_mass, self.reference_temperature, self._temperature[self.inert_gasses[j].temperature_used], self.inert_gasses[j].mass)
				self.inert_gasses[j].catalyst_diffusion = calculate_diffusion_coefficient(self.catalyst_diffusion, self.reference_mass, self.reference_temperature, self._temperature[self.inert_gasses[j].temperature_used], self.inert_gasses[j].mass)
				
		
