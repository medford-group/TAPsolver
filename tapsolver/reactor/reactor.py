
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import numpy as np
import sys

class reactor():
	
	"""
	
	This class acts as a container for all of the reactor parameters
	
	Args:
		**TAPsap Variables**

		zone_lengths (dict): The zone lengths for the reactor.

		zone_void (dict): The bed porosity within each zone.

		reactor_radius (float): The radius of the reactor.
		
		*zone_diffusion (dict): The diffusion coefficient for each zone.

		*zone_residence_time (dict): The residence time within each zone.

		*catalyst_weight (float): The catalyst weight.

		*mol_per_pulse (float): The mol per pulse.

		~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		**TAPsolver Specific Variables**

		mesh (int): The number of cells in the mesh used to solve the PDEs.

		catalyst_mesh_density (int): The number of times the number of mesh cells should be doubled in the catalyst zone.

		temperature (float): The temperature of the reactor (in K).

		output_name (str): The name of the output folder for the results of the analysis. 

		inert_diffusion (float): The reference value of the diffusion coefficient in the inert zone of the reactor.

		catalyst_diffusion (float): The reference value of the diffusion coefficient in the catalyst zone of the reactor. 

		reference_temperature (float): The temperature (in K) that the inert_diff and inert_cat values were determined.

		reference_mass (float): The mass (in a.m.u.) that inert_diff and inert_cat values were determined.

		advection (float): The rate of advection (velocity) of the gas species. One advection value for all gasses.


	"""

	def __init__(self):
		
		self._lengths = {0: 3, 1: 0.06, 2: 3}
		self._voids = {0: 0.4, 1: 0.4, 2: 0.4}
		self._radius = 1.0
		self._length = self._lengths[0] + self._lengths[1] + self._lengths[2]
		self._fractions = {0: self._lengths[0]/self._length,1: self._lengths[1]/self._length, 2: self._lengths[2]/self._length}
		self._cat_center = (self._lengths[0] + self._lengths[1]/2)/self._length
		self._cs_area = 3.14159*(self._radius**2)
		self._t_volumes = {0: self._lengths[0]*self._cs_area,1:self._lengths[1]*self._cs_area,2:self._lengths[2]*self._cs_area}
		self._v_volumes = {0: self._lengths[0]*self._cs_area*self._voids[0],1:self._lengths[1]*self._cs_area*self._voids[1],2:self._lengths[2]*self._cs_area*self._voids[2]}
		self._p_volumes = {0: self._lengths[0]*self._cs_area*(1-self._voids[0]),1:self._lengths[1]*self._cs_area*(1-self._voids[1]),2:self._lengths[2]*self._cs_area*(1-self._voids[2])}
		self._reference_diffusions = {0:16,1:16,2:16} 
	@property
	def lengths(self):
		return self._lengths
	
	@lengths.setter
	def lengths(self,value):
		if type(value) == dict:

			for j in value.keys():
				if value[j] < 0:
					raise ValueError("Zone length dimensions must all be positive (non-negative ")
				else:
					self._lengths[j] = value[j]

		elif type(value) == list:
			self._lengths[value[0]] = value[1] 
		
		if len(self._lengths) > len(self._voids):
			current_void = len(self._voids)
			while len(self._lengths) != len(self._voids):
				self._voids[current_void] = 0.4
				current_void += 1
			
		self._length = 0
		for j in self._lengths.keys():
			self._length += self._lengths[j]



		self._fractions = {}
		for j in self._lengths.keys():
			self._fractions[j] = self._lengths[j]/self._length

		if len(self._lengths.keys()) == 3:
			self._cat_center = (value[0] + value[1]/2)/self._length
		else:
			self._cat_center = 'multiple catalyst zones'
		
		self._t_volumes = {}
		self._v_volumes = {}
		self._p_volumes = {}
		self._reference_diffusions = {}

		for j in self._lengths.keys():
			self._t_volumes[j] = self._lengths[j]*self._cs_area
			self._v_volumes[j] = self._lengths[j]*self._cs_area*self._voids[j]
			self._p_volumes[j] = self._lengths[j]*self._cs_area*(1-self._voids[j])
			self._reference_diffusions[j] = 16

	@property
	def voids(self):
		return self._voids
	@voids.setter
	def voids(self,value):
		if type(value) == dict:
			if len(value.keys()) > len(self._lengths.keys()):
				print("Must add lengths before adding additional voids.")
				sys.exit()
			for j in value.keys():
				if (value[j] < 0) or (value[j] > 1):
					raise ValueError("Zone length dimensions must all be positive (non-negative ")
				else:
					self._voids[j] = value[j]

		elif type(value) == list:
			self._lengths[value[0]] = value[1] 
			if len(value.keys()) > len(self._lengths.keys()):
				print("Must add lengths before adding additional voids.")
				sys.exit()
		## Edit : 
		#self._voids = voids
		
		for j in self._fractions.keys():
			
			
			self._t_volumes = self._lengths[j]*self._cs_area
			self._v_volumes = self._lengths[j]*self._cs_area*self._voids[j]
			self._p_volumes = self._lengths[j]*self._cs_area*(1-self._voids[j])
			self._reference_diffusions[j] = 16

	@property
	def length(self):
		return self._length
	@length.setter
	def length(self,value):
		if value < 0:
			raise ValueError("Length must be greater than 0.")
		scaling_factor = value/self._length
		self._length = value
		
		for j in self._fractions.keys():
			self._lengths[j] = self._lengths[j]*scaling_factor
			
			self._t_volumes = self._lengths[j]*self._cs_areaπ
			self._v_volumes = self._lengths[j]*self._cs_area*self._voids[j]
			self._p_volumes = self._lengths[j]*self._cs_area*(1-self._voids[j])
			self._reference_diffusions[j] = 16

	@property
	def fractions(self):
		return self._fractions

	@property
	def cat_center(self):
		return self._cat_center

	@property
	def radius(self):
		return self._radius
	@radius.setter
	def radius(self,value):
		if value <= 0:
			raise ValueError("Reactor radius must be positive (non-negative ")
		
		self._radius = value
		self._cs_area = (self._radius**2)*3.14159

		for j in self._fractions.keys():
			self._t_volumes = self._lengths[j]*self._cs_area
			self._v_volumes = self._lengths[j]*self._cs_area*self._voids[j]
			self._p_volumes = self._lengths[j]*self._cs_area*(1-self._voids[j])
			self._reference_diffusions[j] = 16

	@property
	def cs_area(self):
		return self._cs_area
	@cs_area.setter
	def cs_area(self,value):
		if value < 0:
			raise ValueError("Cross sectional area must be positive (non-negative ")
		
		self._cs_area = value
		self._radius = np.sqrt(self._cs_area/3.14159)
		for j in self._fractions.keys():
			self._t_volumes = self._lengths[j]*self._cs_area
			self._v_volumes = self._lengths[j]*self._cs_area*self._voids[j]
			self._p_volumes = self._lengths[j]*self._cs_area*(1-self._voids[j])
			self._reference_diffusions[j] = 16

	@property
	def t_volumes(self):
		return self._t_volumes

	@property
	def v_volumes(self):
		return self._v_volumes

	@property
	def p_volumes(self):
		return self._p_volumes

	@property 
	def reference_diffusions(self):
		return self._reference_diffusions
	@reference_diffusions.setter
	def reference_diffusions(self,value):
		if type(value) == dict:
			for j in value.keys():
				if value[j] < 0:
					raise ValueError("Zone length dimensions must all be positive (non-negative ")
				else:
					self._reference_diffusions[j] = value[j]
		elif type(value) == list:
			self._reference_diffusions[value[0]] = value[1] 