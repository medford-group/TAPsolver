from tapsolver import *
print('check reactor')
new_reactor = reactor()
print('pass')
print('check species')
species = reactor_species()
# Define the gasses
CO = define_gas()
CO.mass = 28
CO.intensity = 1
species.add_gas('CO',CO)
s = define_adspecies()
s.concentration = 0
species.add_adspecies('CO*',s)
s = define_adspecies()
s.concentration = 100
species.add_adspecies('*',s)
print('pass')
print('check mechanism')
mech = mechanism()
mech.elementary_processes[0] = elementary_process('CO + * <-> CO*')
mech.elementary_processes[0].forward.k = 1
mech.elementary_processes[0].backward.k = 1
for j in mech.elementary_processes:
	mech.elementary_processes[j].forward.use = 'k'
	try:
		mech.elementary_processes[j].backward.use = 'k'
	except:
		pass
mechanism_constructor(mech)
print('pass')
print('check tapobject')
simulation = TAPobject()
simulation.reactor = new_reactor
simulation.reactor_species = species
simulation.mechanism = mech
print('pass')
print('check forward')
forward_problem(pulse_time=0.01,pulse_number=1,TAPobject_data_original=simulation)
print('pass')
print('')
print('All good!')