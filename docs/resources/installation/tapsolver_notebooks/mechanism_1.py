from tapsolver import *

# Define the reactor
new_reactor = reactor()

species = reactor_species()

# Define the gasses
CO = define_gas()
O2 = define_gas()
CO2 = define_gas()

CO.mass = 28
CO.intensity = 1
CO.sigma = 0.2

O2.mass = 32
O2.intensity = 1
O2.sigma = 0.2

CO2.mass = 44
CO2.intensity = 0
CO2.sigma = 0.2

species.add_gas('CO',CO)
species.add_gas('O2',O2)
species.add_gas('CO2',CO2)

Ar = define_gas()

Ar.mass = 40
Ar.intensity = 1
Ar.noise = 0.2

species.add_inert_gas('Ar',Ar)

s = define_adspecies()
s.concentration = 0
species.add_adspecies('O*',s)

s = define_adspecies()
s.concentration = 100
species.add_adspecies('*',s)

mech = mechanism()

mech.elementary_processes[0] = elementary_process('O2 + 2* <-> 2O*')
mech.elementary_processes[0].forward.k = 1
mech.elementary_processes[0].backward.k = 1

mech.elementary_processes[1] = elementary_process('CO + O* <-> CO2 + *')
mech.elementary_processes[1].forward.k = 1
mech.elementary_processes[1].backward.k = 0

for j in mech.elementary_processes:
	mech.elementary_processes[j].forward.use = 'k'
	try:
		mech.elementary_processes[j].backward.use = 'k'
	except:
		pass

mechanism_constructor(mech)

simulation = TAPobject()
simulation.reactor = new_reactor
simulation.reactor_species = species
simulation.mechanism = mech

simulation.output_name = 'fit_new'
simulation.data_name = './exp_new/TAP_experimental_data.json'
simulation.gasses_objective = ['CO','O2','CO2']
simulation.parameters_of_interest = ['TAPobject_data.mechanism.elementary_processes[0].forward.k','TAPobject_data.mechanism.elementary_processes[0].backward.k','TAPobject_data.mechanism.elementary_processes[1].forward.k','TAPobject_data.mechanism.elementary_processes[1].backward.k']
simulation.optimize = False

save_object(simulation,'./mech_1.json')
simulation2 = update_parameters(simulation)
forward_problem(0.3,1,simulation2)

flux_graph(simulation2)
