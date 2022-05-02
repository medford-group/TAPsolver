import jsonpickle
from tapsolver import *

new_reactor = reactor()

#new_reactor.zone_lengths = {0:1,1:0.1,2:.3,3:0.1,4:1}
#new_reactor.zone_voids = {0:0.4,1:0.4,2:0.4,3:0.4,3:0.4}
#new_reactor.catalyst_locations = [0,1,0,1,0]

new_reactor.zone_lengths = {0:1.15,1:0.2,2:1.15}
new_reactor.zone_voids = {0:0.4,1:0.4,2:0.4}
new_reactor.catalyst_locations = [0,1,0]

#new_reactor.zone_lengths = {0:1.85,1:0.36,2:1,3:0.36,4:1.85}
#new_reactor.zone_voids = {0:0.4,1:0.4,2:0.4,3:0.4,3:0.4}
#new_reactor.catalyst_locations = [0,1,0,1,0]

### Add the reactor gasses (reactive and inert)
new_reactor_species = reactor_species()
new_reactor_species.inert_diffusion = 25.04
new_reactor_species.catalyst_diffusion = 25.04
new_reactor_species.reference_temperature = 700
new_reactor_species.reference_mass = 40
new_reactor_species.temperature = 833
new_reactor_species.advection = 0

C3H8 = define_gas()
C3H8.mass = 44.1
C3H8.intensity = 6.870029803076434
C3H8.delay = 0.0
C3H8.sigma = 0.37
C3H8.noise = 0.37
new_reactor_species.add_gas('C3H8',C3H8)

H2 = define_gas()
H2.mass = 2
H2.intensity = 0.0
H2.delay = 0.0
H2.sigma = 1
H2.noise = 1
new_reactor_species.add_gas('H2',H2)
#
C3H6 = define_gas()
C3H6.mass = 42.08
C3H6.intensity = 0.0
C3H6.delay = 0.0
C3H6.sigma = 1
C3H6.noise = 1
new_reactor_species.add_gas('C3H6',C3H6)

Ar = define_gas()
Ar.mass = 40.1
Ar.intensity = 6.870029803076434
Ar.delay = 0.0
Ar.noise = 1
Ar.sigma = 1
new_reactor_species.add_inert_gas('Ar',Ar)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H8*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('H2*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H6*',s)

s = define_adspecies()
s.concentration = 240
new_reactor_species.add_adspecies('*',s)

new_mechanism = mechanism()

# Mechanism #1
new_mechanism.elementary_processes[0] = elementary_process('C3H8 + * <-> C3H8*')
new_mechanism.elementary_processes[1] = elementary_process('C3H8* + * <-> C3H6* + H2*')
new_mechanism.elementary_processes[2] = elementary_process('C3H6* <-> C3H6 + *')
new_mechanism.elementary_processes[3] = elementary_process('H2* <-> H2 + *')

### Mechanism #2
#new_mechanism.elementary_processes[0] = elementary_process('C3H8 + * <-> C3H8*')
#new_mechanism.elementary_processes[] = elementary_process('C3H8* + O* <-> C3H7* + OH*')
#new_mechanism.elementary_processes[] = elementary_process('C3H7* + O* <-> C3H6* + OH*')
#new_mechanism.elementary_processes[] = elementary_process('C3H6* <-> C3H6 + *')
#new_mechanism.elementary_processes[] = elementary_process('2OH* <-> H2 + 2O*')

new_mechanism.elementary_processes[0].forward.k = 0.0001
new_mechanism.elementary_processes[0].backward.k = 0.0001
new_mechanism.elementary_processes[1].forward.k = 0.0001
new_mechanism.elementary_processes[1].backward.k = 0.0001
new_mechanism.elementary_processes[2].forward.k = 0.0001
new_mechanism.elementary_processes[2].backward.k = 0.0001
new_mechanism.elementary_processes[3].forward.k = 0.0001
new_mechanism.elementary_processes[3].backward.k = 0.0001

for j in new_mechanism.elementary_processes:
	new_mechanism.elementary_processes[j].forward.use = 'k'
	try:
		new_mechanism.elementary_processes[j].backward.use = 'k'
	except:
		pass

mechanism_constructor(new_mechanism)

pdh_1 = TAPobject()
pdh_1.reactor = new_reactor
pdh_1.reactor_species = new_reactor_species
pdh_1.mechanism = new_mechanism

pdh_1.show_graph = False
pdh_1.output_name = 'exp_new'

pdh_1.data_name = None
print(pdh_1.reactor.total_length)
forward_problem(0.5,1,pdh_1)