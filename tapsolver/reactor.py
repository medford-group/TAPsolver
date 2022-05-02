
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

	def __init__(self, zone_lengths = {0: 3, 1: 0.06, 2: 3}, zone_voids = {0: 0.4, 1: 0.4, 2: 0.4}, reactor_radius = 1.0, total_length = None,catalyst_locations=[0,1,0], catalyst_center_fraction=0):
		
		self.zone_lengths = zone_lengths
		self.zone_voids = zone_voids
		self.reactor_radius = reactor_radius
		self.catalyst_locations = catalyst_locations
		
	@property
	def zone_lengths(self):
		return self._zone_lengths
	@zone_lengths.setter
	def zone_lengths(self,value):
		for j in value.keys():
			if value[j] < 0:#if value[0] <= 0 or value[1] < 0 or value[2] <= 0:
				raise ValueError("Zone length dimensions must all be positive (non-negative ")
		self._zone_lengths = {0: value[0], 1: value[1], 2: value[2]}
		new_value = 0
		for j in value.keys():
			new_value += value[j]
		self._total_length = new_value
		self._length_fractions = []
		for j in value.keys():
			self._length_fractions.append(value[j]/self.total_length)
		self._catalyst_center_fraction = []
		temp_current_length = 0
		for j in value.keys():
			self._catalyst_center_fraction.append((temp_current_length + value[j]/2)/self.total_length)
			temp_current_length += value[j]

	@property
	def total_length(self):
		return self._total_length
	@total_length.setter
	def total_length(self, value):
		for j in value.keys():
			print(j)
			print(self._total_length)
			self._total_length += self.zone_lengths[j]# + self.zone_lengths[1] + self.zone_lengths[2]

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
		if value == (float or int):
			if value <= 0:
				raise ValueError("Reactor radius must be positive (non-negative ")
			self._reactor_radius = value
			self._cross_sectional_radius = (value**2)*3.14159
		else:
			self._reactor_radius = value
			self._cross_sectional_radius = (value**2)*3.14159

	@property
	def cross_sectional_radius(self):
		return self._cross_sectional_radius
	@cross_sectional_radius.setter
	def cross_sectional_radius(self,value):
		self._cross_sectional_radius = (self._reactor_radius**2)*3.14159