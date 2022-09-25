from .TAPobject import TAPobject
from .reactor_species import reactor_species
from .reactor import reactor
from .mechanism import mechanism
from .define_gas import define_gas
from .define_adspecies import define_adspecies
#from structures import TAPobject,reactor_species,reactor,mechanism, define_gas, define_adspecies
import jsonpickle
import json
import sys

def read_TAPobject(file_name):

	loaded_TAPobject = TAPobject()

	with open(file_name) as f:
		data = json.loads(f.read())
	f.close()

	sameObject = jsonpickle.decode(data)
	sameObject2 = sameObject["1"]

	# REQUIRED FOR FORWARD AND INVERSE ANALYSIS
	loaded_TAPobject.reactor = reactor()
	loaded_TAPobject.reactor.zone_lengths = {0:sameObject2.reactor.zone_lengths['0'], 1: sameObject2.reactor.zone_lengths['1'], 2: sameObject2.reactor.zone_lengths['2']}
	loaded_TAPobject.reactor.zone_voids = {0:sameObject2.reactor.zone_voids['0'], 1: sameObject2.reactor.zone_voids['2'], 2: sameObject2.reactor.zone_voids['2']}
	loaded_TAPobject.reactor.reactor_radius = sameObject2.reactor.reactor_radius	
	
	loaded_TAPobject.mechanism = mechanism()
	loaded_TAPobject.mechanism.rate_array = sameObject2.mechanism.rate_array
	loaded_TAPobject.mechanism.reactants = sameObject2.mechanism.reactants 
	loaded_TAPobject.mechanism.reactions = sameObject2.mechanism.reactions
	
	for jnum,j in enumerate(sameObject2.mechanism.elementary_processes):
		loaded_TAPobject.mechanism.elementary_processes[jnum] = sameObject2.mechanism.elementary_processes[j]

	for jnum,j in enumerate(sameObject2.mechanism.kinetic_links):
		loaded_TAPobject.mechanism.kinetic_links[jnum] = sameObject2.mechanism.kinetic_links[j]		

	loaded_TAPobject.reactor_species = sameObject2.reactor_species

########################################

	for j in sameObject2.reactor_species.gasses:
		new_gas = define_gas()
		new_gas.inert_diffusion = sameObject2.reactor_species.gasses[j].inert_diffusion
		new_gas.catalyst_diffusion = sameObject2.reactor_species.gasses[j].catalyst_diffusion
		new_gas.intensity = sameObject2.reactor_species.gasses[j].intensity
		new_gas.delay = sameObject2.reactor_species.gasses[j].delay
		new_gas.noise = sameObject2.reactor_species.gasses[j].noise
		new_gas.sigma = sameObject2.reactor_species.gasses[j].sigma
		new_gas.temperature_used = sameObject2.reactor_species.gasses[j].temperature_used
		new_gas.initial_concentration = sameObject2.reactor_species.gasses[j].initial_concentration
		new_gas.inlet_concentration = sameObject2.reactor_species.gasses[j].inlet_concentration
		new_gas.mass = sameObject2.reactor_species.gasses[j].mass

		loaded_TAPobject.reactor_species.add_gas(j,new_gas)

	for j in sameObject2.reactor_species.inert_gasses:
		new_gas = define_gas()
		new_gas.inert_diffusion = sameObject2.reactor_species.inert_gasses[j].inert_diffusion
		new_gas.catalyst_diffusion = sameObject2.reactor_species.inert_gasses[j].catalyst_diffusion
		new_gas.intensity = sameObject2.reactor_species.inert_gasses[j].intensity
		new_gas.delay = sameObject2.reactor_species.inert_gasses[j].delay
		new_gas.noise = sameObject2.reactor_species.inert_gasses[j].noise
		new_gas.sigma = sameObject2.reactor_species.inert_gasses[j].sigma
		new_gas.temperature_used = sameObject2.reactor_species.inert_gasses[j].temperature_used
		new_gas.initial_concentration = sameObject2.reactor_species.inert_gasses[j].initial_concentration
		new_gas.inlet_concentration = sameObject2.reactor_species.inert_gasses[j].inlet_concentration
		new_gas.mass = sameObject2.reactor_species.inert_gasses[j].mass

		loaded_TAPobject.reactor_species.add_inert_gas(j,new_gas)

	for j in sameObject2.reactor_species.adspecies:
		new_adspecies = define_adspecies()
		new_adspecies.concentration = sameObject2.reactor_species.adspecies[j].concentration
		new_adspecies.noise = sameObject2.reactor_species.adspecies[j].noise
		new_adspecies.sigma = sameObject2.reactor_species.adspecies[j].sigma
		loaded_TAPobject.reactor_species.add_adspecies(j,new_adspecies)

	loaded_TAPobject.reactor_species.inert_diffusion = sameObject2.reactor_species.inert_diffusion
	loaded_TAPobject.reactor_species.catalyst_diffusion = sameObject2.reactor_species.catalyst_diffusion
	loaded_TAPobject.reactor_species.reference_temperature = sameObject2.reactor_species.reference_temperature
	loaded_TAPobject.reactor_species.reference_mass = sameObject2.reactor_species.reference_mass
	print('read temperature')
	print(type(sameObject2.reactor_species.temperature))
	print(sameObject2.reactor_species.temperature)
	if type(sameObject2.reactor_species.temperature) == dict:
		temp_temp_dictionary = {}
		for j in sameObject2.reactor_species.temperature.keys():
			temp_temp_dictionary[int(j)] = sameObject2.reactor_species.temperature[j]
		loaded_TAPobject.reactor_species.temperature = temp_temp_dictionary
	else:
		loaded_TAPobject.reactor_species.temperature = sameObject2.reactor_species.temperature
	loaded_TAPobject.reactor_species.advection = sameObject2.reactor_species.advection
	loaded_TAPobject.ftol = sameObject2.ftol
	loaded_TAPobject.gtol = sameObject2.gtol
	loaded_TAPobject.reactor_species.reference_pulse_size = sameObject2.reactor_species.reference_pulse_size

	########################################

	# REQUIRED FOR INVERSE ANALYSIS
	loaded_TAPobject.experimental_data = sameObject2.experimental_data

	# Simulation precision preferences
	loaded_TAPobject.mesh = sameObject2.mesh
	loaded_TAPobject.data_storage_frequency = sameObject2.data_storage_frequency
	loaded_TAPobject.catalyst_mesh_density = sameObject2.catalyst_mesh_density
	loaded_TAPobject.optimization_method = sameObject2.optimization_method
	
	loaded_TAPobject.parameter_scale = sameObject2.parameter_scale

	# Data storage preferences
	loaded_TAPobject.output_name = sameObject2.output_name
	loaded_TAPobject.derivative_name = sameObject2.derivative_name
	loaded_TAPobject.data_name = sameObject2.data_name
	loaded_TAPobject.store_flux_data = sameObject2.store_flux_data
	loaded_TAPobject.store_thin_data = sameObject2.store_thin_data
	loaded_TAPobject.gas_noise = sameObject2.gas_noise
	loaded_TAPobject.surface_noise = sameObject2.surface_noise
	loaded_TAPobject.catalyst_data_type = sameObject2.catalyst_data_type

	# Objective function preferences
	loaded_TAPobject.objective = sameObject2.objective
	loaded_TAPobject.gasses_objective = sameObject2.gasses_objective # Provide the names of the gasses to include in these objectives
	loaded_TAPobject.inert_gasses_objective = sameObject2.inert_gasses_objective
	loaded_TAPobject.adspecies_objective = sameObject2.adspecies_objective
	loaded_TAPobject.thermodynamic_constraints = sameObject2.thermodynamic_constraints

	# Objective function preferences
	loaded_TAPobject.parameters_of_interest = sameObject2.parameters_of_interest

	# Sensitivity Analysis
	loaded_TAPobject.tangent_linear_sensitivity = sameObject2.tangent_linear_sensitivity
	loaded_TAPobject.finite_difference_trans_sensitivty = sameObject2.finite_difference_trans_sensitivty
	loaded_TAPobject.adjoint_sensitivitiy = sameObject2.adjoint_sensitivitiy
	loaded_TAPobject.optimize = sameObject2.optimize
	loaded_TAPobject.uncertainty = sameObject2.uncertainty
	loaded_TAPobject.objective_return = sameObject2.uncertainty

	# Flux graph name
	loaded_TAPobject.pulses_graphed = sameObject2.pulses_graphed
	loaded_TAPobject.display_analytical = sameObject2.display_analytical
	loaded_TAPobject.scaled_graph = sameObject2.scaled_graph
	loaded_TAPobject.display_objective = sameObject2.display_objective
	loaded_TAPobject.show_graph = sameObject2.show_graph
	loaded_TAPobject.store_graph = sameObject2.store_graph

	return loaded_TAPobject