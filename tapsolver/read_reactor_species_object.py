from reactor_species import reactor_species
from define_gas import define_gas
from define_adspecies import define_adspecies
#from structures import reactor_species, define_gas, define_adspecies
import jsonpickle
import json
import sys

def read_reactor_species_object(file_name):
	loaded_reactor_species = reactor_species()

	with open(file_name) as f:
		data = json.loads(f.read())
	f.close()

	sameObject = jsonpickle.decode(data)
	sameObject2 = sameObject["1"]

	for j in sameObject2.gasses:
		print(j)
		new_gas = define_gas()
		new_gas.inert_diffusion = sameObject2.gasses[j].inert_diffusion
		new_gas.catalyst_diffusion = sameObject2.gasses[j].catalyst_diffusion
		new_gas.intensity = sameObject2.gasses[j].intensity
		new_gas.delay = sameObject2.gasses[j].delay
		new_gas.initial_concentration = sameObject2.gasses[j].initial_concentration
		new_gas.inlet_concentration = sameObject2.gasses[j].inlet_concentration
		new_gas.mass = sameObject2.gasses[j].mass

		loaded_reactor_species.add_gas(j,new_gas)

	for j in sameObject2.inert_gasses:
		print(j)
		new_gas = define_gas()
		new_gas.inert_diffusion = sameObject2.inert_gasses[j].inert_diffusion
		new_gas.catalyst_diffusion = sameObject2.inert_gasses[j].catalyst_diffusion
		new_gas.intensity = sameObject2.inert_gasses[j].intensity
		new_gas.delay = sameObject2.inert_gasses[j].delay
		new_gas.initial_concentration = sameObject2.inert_gasses[j].initial_concentration
		new_gas.inlet_concentration = sameObject2.inert_gasses[j].inlet_concentration
		new_gas.mass = sameObject2.inert_gasses[j].mass

		loaded_reactor_species.add_inert_gas(j,new_gas)

	for j in sameObject2.adspecies:
		print(j)
		new_adspecies = define_adspecies()
		new_adspecies.concentration = sameObject2.adspecies[j].concentration

		loaded_reactor_species.add_adspecies(j,new_adspecies)



	loaded_reactor_species.inert_diffusion = sameObject2.inert_diffusion
	loaded_reactor_species.catalyst_diffusion = sameObject2.catalyst_diffusion
	loaded_reactor_species.reference_temperature = sameObject2.reference_temperature
	loaded_reactor_species.reference_mass = sameObject2.reference_mass
	loaded_reactor_species.temperature = sameObject2.temperature
	loaded_reactor_species.advection = sameObject2.advection
	loaded_reactor_species.reference_pulse_size = sameObject2.reference_pulse_size

	return loaded_reactor_species