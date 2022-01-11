from tapsolver import *
import numpy as np
import os

new_reactor = reactor()

new_reactor.zone_lengths = {0:1.8675,1:0.2,2:1.8675}
new_reactor.zone_voids = {0:0.4,1:0.4,2:0.4}
new_reactor.reactor_radius = 0.19

### Add the reactor gasses (reactive and inert)
new_reactor_species = reactor_species()
new_reactor_species.inert_diffusion = 27
new_reactor_species.catalyst_diffusion = 27
new_reactor_species.reference_temperature = 700
new_reactor_species.reference_mass = 40
new_reactor_species.temperature = 700
new_reactor_species.advection = 0

C3H8 = define_gas()
C3H8.mass = 44.1
C3H8.intensity = 1.0
C3H8.delay = 0.0
C3H8.sigma = 0.17
C3H8.noise = 0.17
new_reactor_species.add_gas('C3H8',C3H8)

O2 = define_gas()
O2.mass = 32
O2.intensity = 1.0
O2.delay = 0.0
O2.sigma = 0.17
O2.noise = 0.17
new_reactor_species.add_gas('O2',O2)
#
C3H6 = define_gas()
C3H6.mass = 42.08
C3H6.intensity = 0.0
C3H6.delay = 0.0
C3H6.sigma = 0.17
C3H6.noise = 0.17
new_reactor_species.add_gas('C3H6',C3H6)

H2O = define_gas()
H2O.mass = 18
H2O.intensity = 0.0
H2O.delay = 0
H2O.sigma = 0.17
H2O.noise = 0.17
new_reactor_species.add_gas('H2O',H2O)

CO2 = define_gas()
CO2.mass = 44.01
CO2.intensity = 0.0
CO2.delay = 0
CO2.sigma = 0.17
CO2.noise = 0.17
new_reactor_species.add_gas('CO2',CO2)

Ar = define_gas()
Ar.mass = 44.1
Ar.intensity = 1
Ar.delay = 0.0
Ar.noise = 0.17
Ar.sigma = 0.17
new_reactor_species.add_inert_gas('Ar',Ar)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H8*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O^',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H6*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('H*',s)

s = define_adspecies()
s.concentration = 60
new_reactor_species.add_adspecies('*',s)

s = define_adspecies()
s.concentration = 60
new_reactor_species.add_adspecies('^',s)

new_mechanism = mechanism()

new_mechanism.elementary_processes[0] = elementary_process('C3H8 + * <-> C3H8*') # Adsorption
new_mechanism.elementary_processes[1] = elementary_process('O2 + 2^ <-> 2O^') # Adsorption
new_mechanism.elementary_processes[2] = elementary_process('C3H6 + * <-> C3H6*') # Adsorption
new_mechanism.elementary_processes[3] = elementary_process('C3H8* + O^ <-> C3H6* + H2O')
new_mechanism.elementary_processes[4] = elementary_process('C3H8* + 9O^ <-> 3CO2 + 4H2O')
new_mechanism.elementary_processes[5] = elementary_process('C3H6* + 9O^ <-> 3CO2 + 3H2O')

# Free energy figure
new_mechanism.elementary_processes[0].forward.dG = -0.75
new_mechanism.elementary_processes[0].forward.Ga = 1.45
new_mechanism.elementary_processes[1].forward.dG = -1.5
new_mechanism.elementary_processes[1].forward.Ga = 1.35
new_mechanism.elementary_processes[2].forward.dG = -.35
new_mechanism.elementary_processes[2].forward.Ga = 1.1
new_mechanism.elementary_processes[3].forward.dG = -.61
new_mechanism.elementary_processes[3].forward.Ga = 1.5#0.56
new_mechanism.elementary_processes[4].forward.dG = -.89
new_mechanism.elementary_processes[4].forward.Ga = 1.09 #1.29
new_mechanism.elementary_processes[5].forward.dG = -0.89
new_mechanism.elementary_processes[5].forward.Ga = 0.65

## Make a function
for j in new_mechanism.elementary_processes:
	new_mechanism.elementary_processes[j].forward.use = 'G'
	#new_mechanism.elementary_processes[j].backward.use = 'G'
try:
	new_mechanism.elementary_processes[j].backward.use = 'G'
except:
	pass

mechanism_constructor(new_mechanism)

opdh_2 = TAPobject()
opdh_2.mechanism = new_mechanism
opdh_2.reactor_species = new_reactor_species
opdh_2.reactor = new_reactor
opdh_2.show_graph = False

opdh_2.parameter_scale = 60

opdh_2.output_name = 'exp_new_2'
