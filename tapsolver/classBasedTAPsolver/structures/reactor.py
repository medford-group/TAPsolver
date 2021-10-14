
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

	def __init__(self):
		self.zone_lengths = {'zone0': 2.80718, 'zone1': 0.17364, 'zone2': 2.80718}
		
		self.zone_voids = {'zone0': 0.4, 'zone1': 0.4, 'zone2': 0.4}
		
		self.reactor_radius = 1.0
		
		self.temperature = 385.65
		
		self.mesh = 200
		
		self.catalyst_mesh_density = 4
		
		self.output_name = 'exp_new'
		
		self.inert_diffusion = 16
		
		self.catalyst_diffusion = 16
		
		self.reference_temperature = 385.6
		
		self.reference_mass = 40
		
		self.advection = 0

		self.mesh_size = 200

	"""
	
	This function just converts the input lengths to fractions.

	Args:
		reactor (class Reactor): The reactor information.

	Returns:
		zone_fractions (list of floats): The length fraction of each zone in the reactor.

	"""

	def lengthFractions(self):

		return [self.length[0]/sum(self.length),self.length[1]/sum(self.length),self.length[2]/sum(self.length)]

	"""
	
	This function calculates the central point of the catalyst zone.

	Args:
		reactor (class Reactor): The reactor information.

	Return:
		zone_fractions (float): The central point of the catalyst zone.

	"""

	def reactorCenterFraction(self):
		
		return (self.length[0] + self.length[0]/2)/sum(self.length)


#### Include details

	def cross_sectional_area(self):
		return (self.reactor_radius**2)*3.14159

	def total_length(self):
		return self.zone_lengths['zone0'] + self.zone_lengths['zone1'] + self.zone_lengths['zone2']

	def mesh_step_size(self):
		return self.total_length()/mesh_size

	def point_volume(self):
		return self.mesh_step_size()*cross_sectional_area()*self.zone_voids['zone0']