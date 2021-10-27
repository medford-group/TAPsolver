
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved


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

	def __init__(self, zone_lengths = {0: 2.80718, 1: 0.17364, 2: 2.80718}, zone_voids = {0: 0.4, 1: 0.4, 2: 0.4}, reactor_radius = 1.0, temperature = 385.65, total_length = None, catalyst_center_fraction=0):
		
		self.zone_lengths = zone_lengths
		self.zone_voids = zone_voids
		self.reactor_radius = reactor_radius
		self.temperature = temperature

	@property
	def zone_lengths(self):
		return self._zone_lengths
	@zone_lengths.setter
	def zone_lengths(self,value):
		if value[0] <= 0 or value[1] < 0 or value[2] <= 0:
			raise ValueError("Zone length dimensions must all be positive (non-negative ")
		self._zone_lengths = {0: value[0], 1: value[1], 2: value[2]}
		self._total_length = value[0] + value[1] + value[2]
		self._length_fractions = [value[0]/self.total_length,value[1]/self.total_length,value[2]/self.total_length]
		self.catalyst_center_fraction = (value[0] + value[1]/2)/self.total_length

	@property
	def total_length(self):
		return self._total_length
	@total_length.setter
	def total_length(self, value):
		self._total_length = self.zone_lengths[0] + self.zone_lengths[1] + self.zone_lengths[2]

	@property
	def length_fractions(self):
		return self._length_fractions
	@length_fractions.setter
	def length_fractions(self,value):
		self._length_fractions = [self.zone_lengths[0]/self.total_length,self.zone_lengths[1]/self.total_length,self.zone_lengths[2]/self.total_length]

	@property
	def catalyst_center_fraction(self):
		return self._catalyst_center_fraction
	@catalyst_center_fraction.setter
	def catalyst_center_fraction(self,value):
		self._catalyst_center_fraction = (self.zone_lengths[0] + self.zone_lengths[1]/2)/self.total_length

	@property
	def reactor_radius(self):
		return self._reactor_radius
	@reactor_radius.setter
	def reactor_radius(self,value):
		if value <= 0:
			raise ValueError("Reactor radius must be positive (non-negative ")
		self._reactor_radius = value
		self._cross_sectional_radius = (value**2)*3.14159

	@property
	def cross_sectional_radius(self):
		return self._cross_sectional_radius
	@cross_sectional_radius.setter
	def cross_sectional_radius(self,value):
		self._cross_sectional_radius = (self._reactor_radius**2)*3.14159
