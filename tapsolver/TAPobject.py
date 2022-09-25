
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class TAPobject():

	"""
	
	Description: This class acts as a container for the three necessary components 
	to run a TAP simulation: the reactor, mechanism and reactor_species 
	objects. It also lets the user specify the experimental data set, 
	which can be used for the inverse problem or simple data visualization.

	The user can also easily adjust simulation precision and inverse 
	specifications.

	Although this class can appear overloaded with parameters, it is much 
	easier and accessable than specifying and processing user preferences 
	on a per function basis (i.e. the former implementation in TAPsolver).
	
	Attributes:
		
		reactor (reactor object): The reactor details.

		mechanism (mechanism object): The elementary processes and their 
		associated kinetics.

		reactor_species (reactor_species object): The details of the 
		reactor_species + transport	+ initial conditions.

		experimental_data (experimental_data object): The experimental or
		synthetic data desired for analysis (or visualization).
		
		mesh (int): The number of cells in the FEM simulation.

		catalyst_mesh_density (int): The number of times the number of 
		catalyst mesh cells should be doubled. For example, if the number
		of catalyst cells is initially 4 in the catalyst zone and the 
		catalyst_mesh_density is 3, the new number of cells in the catalyst
		zone will be 4^(catalyst_mesh_density) = 4^(3) = 64. The number of 
		cells in each inert zone will remain the same. 

		output_name (str): The name of the directory you want to generate
		for storing any generated information.

		data_name (str): the file location of the experimental data. 

		store_flux_data (boolean): Store the outlet flux data or not. 

		store_catalyst_data (boolean): Store the catalyst zone data or not 
		(i.e. gas and adspecies concentrations).

		catalyst_data_type (str): The type of catalyst zone data to be stored. 

		gas_noise (boolean): Specify if noise should be included during a simulation.

		surface_noise (boolean): Specify if noise should be included in surface data during a simulation. 

		catalyst_data_type (boolean): define the way in which catalyst data is stored. 

		objective (str): The objective function to be used during the inverse
		problem.

		gasses_objective (list): The reactive gasses that should be included 
		in the objective function.

		inert_gasses_objective (list): The inert gasses that should be included 
		in the objective function.

		adspecies_objective (list): The adspecies that should be included 
		in the objective function.

		parameters_of_interest (list): specify the list of parameters to be optimized.

		finite_difference_trans_sensitivity (boolean): Calculate the dynamic sensitivity matrix.

		optimize (boolean): Fit the parameters_of_interest to the objective_function.

		uncertainty (boolean): Calculate the uncertainty around the local minimum.

		objective_return (boolean): Print the value of the objective function.

		  

	"""

	def __init__(self):
		
		# REQUIRED FOR FORWARD AND INVERSE ANALYSIS
		self.reactor = None
		self.mechanism = None
		self.reactor_species = None
		
		# REQUIRED FOR INVERSE ANALYSIS
		self.experimental_data = None

		# Simulation precision preferences
		self.mesh = 400
		self.catalyst_mesh_density = 3
		self.parameter_scale = 0

		# Data storage preferences
		self.output_name = 'exp_new'
		self.derivative_name = 'exp_new'
		self.data_name = None
		self.store_flux_data = True
		self.store_thin_data = True
		self.gas_noise = True
		self.surface_noise = True
		self.catalyst_data_type = 'single_point'

		# Objective function preferences
		self.objective = ''
		self.gasses_objective = []
		self.inert_gasses_objective = []
		self.adspecies_objective = []
		self.thermodynamic_constraints = False
		self.thermodynamic_alpha = 1
		self.thermodynamic_free_energy = 0
		self.thermo_equations = ['']

		# Objective function preferences
		self.parameters_of_interest = []

		# Sensitivity Analysis
		self.tangent_linear_sensitivity = False
		self.finite_difference_trans_sensitivty = False
		self.adjoint_sensitivitiy = False
		self.optimize = False
		self.uncertainty = False
		self.objective_return = False 
		self.ftol = 1e-22
		self.gtol = 1e-22

		# Flux graph name
		self.pulses_graphed = 1
		self.display_analytical = False
		self.scaled_graph = False
		self.display_objective = False
		self.show_graph = False
		self.store_graph = False

		self.data_storage_frequency = 1

		self.optimization_method = 'L-BFGS-B'