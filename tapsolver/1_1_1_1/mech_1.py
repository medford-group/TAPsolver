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
C3H8.delay = 0
C3H8.sigma = 0.1
C3H8.noise = 0.1
new_reactor_species.add_gas('C3H8',C3H8)

#O2 = define_gas()
#O2.mass = 32
#O2.intensity = 1.0
#O2.delay = 0
#O2.sigma = 0.1
#O2.noise = 0.1
#new_reactor_species.add_gas('O2',O2)
#
#C3H6 = define_gas()
#C3H6.mass = 42.08
#C3H6.intensity = 0.0
#C3H6.delay = 0
#C3H6.sigma = 0.1
#C3H6.noise = 0.1
#new_reactor_species.add_gas('C3H6',C3H6)
#
#H2O = define_gas()
#H2O.mass = 18
#H2O.intensity = 0.0
#H2O.delay = 0
#H2O.sigma = 0.1
#H2O.noise = 0.1
#new_reactor_species.add_gas('H2O',H2O)
#
#CO2 = define_gas()
#CO2.mass = 44.01
#CO2.intensity = 0.0
#CO2.delay = 0
#CO2.sigma = 0.0
#CO2.noise = 0.0
#new_reactor_species.add_gas('CO2',CO2)

Ar = define_gas()
Ar.mass = 44.1
Ar.intensity = 1
Ar.delay = 0.0
Ar.noise = 0.1
Ar.sigma = 0.1
new_reactor_species.add_inert_gas('Ar',Ar)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H8*',s)

#s = define_adspecies()
#s.concentration = 0
#new_reactor_species.add_adspecies('O*',s)
#
#s = define_adspecies()
#s.concentration = 0
#new_reactor_species.add_adspecies('C3H6*',s)
#
#s = define_adspecies()
#s.concentration = 0
#new_reactor_species.add_adspecies('H2O*',s)
#
#s = define_adspecies()
#s.concentration = 0
#new_reactor_species.add_adspecies('CO2*',s)
#
#s = define_adspecies()
#s.concentration = 0
#new_reactor_species.add_adspecies('H*',s)

s = define_adspecies()
s.concentration = 120
new_reactor_species.add_adspecies('*',s)

new_mechanism = mechanism()

new_mechanism.elementary_processes[0] = elementary_process('C3H8 + * <-> C3H8*')
#new_mechanism.elementary_processes[1] = elementary_process('O2 + 2* <-> 2O*')
#new_mechanism.elementary_processes[2] = elementary_process('C3H6 + * <-> C3H6*')
#new_mechanism.elementary_processes[3] = elementary_process('H2O + * <-> H2O*')
#new_mechanism.elementary_processes[4] = elementary_process('O182 + 2* <-> 2O18*')
#new_mechanism.elementary_processes[5] = elementary_process('CO2 + * <-> CO2*')
#new_mechanism.elementary_processes[6] = elementary_process('C3H8* + 2* <-> C3H6* + 2H*')
#new_mechanism.elementary_processes[7] = elementary_process('O* + 2H* <-> H2O*')
#new_mechanism.elementary_processes[8] = elementary_process('C3H8* + 10O* <-> 3CO2* + 4H2O*')
#new_mechanism.elementary_processes[9] = elementary_process('C3H6* + 9O* <-> 3CO2* + 3H2O*')

new_mechanism.elementary_processes[0].forward.k = 1#2*2.7e-19
new_mechanism.elementary_processes[0].backward.k = 1#2.7e-19
#new_mechanism.elementary_processes[0].forward.dG = -.2#2*2.7e-19
#new_mechanism.elementary_processes[0].forward.Ga = 0.4#2.7e-19

# T = K
# k = ev / K
# h = ev * s
# reactions / s
# (6.022e15)
#print((((8.617e-5)*700)/((6.582e-16)*(120)))*np.exp(-new_mechanism.elementary_processes[0].forward.Ga/((8.617e-5)*700)))
#print((((8.617e-5)*700)/((6.582e-16)))*np.exp(-(new_mechanism.elementary_processes[0].forward.Ga-new_mechanism.elementary_processes[0].forward.dG)/((8.617e-5)*700)))

#sys.exit()
#new_mechanism.elementary_processes[0].backward.k = 0
#new_mechanism.elementary_processes[1].forward.dG = 0
#new_mechanism.elementary_processes[1].forward.Ga = 10
#new_mechanism.elementary_processes[2].forward.dG = 0
#new_mechanism.elementary_processes[2].forward.Ga = 10
#new_mechanism.elementary_processes[3].forward.dG = 0
#new_mechanism.elementary_processes[3].forward.Ga = 10
#new_mechanism.elementary_processes[4].forward.dG = 0
#new_mechanism.elementary_processes[4].forward.Ga = 10
#new_mechanism.elementary_processes[5].forward.dG = 0
#new_mechanism.elementary_processes[5].forward.Ga = 10
#new_mechanism.elementary_processes[6].forward.dG = 0
#new_mechanism.elementary_processes[6].forward.Ga = 10
#new_mechanism.elementary_processes[7].forward.dG = 0 
#new_mechanism.elementary_processes[7].forward.Ga = 10
#new_mechanism.elementary_processes[8].forward.dG = 0 
#new_mechanism.elementary_processes[8].forward.Ga = 10
#new_mechanism.elementary_processes[9].forward.dG = 0 
#new_mechanism.elementary_processes[9].forward.Ga = 10

#print(new_mechanism.elementary_processes[0].forward.use)
#print(new_mechanism.elementary_processes[0].backward.use)

## Make a function
for j in new_mechanism.elementary_processes:
	new_mechanism.elementary_processes[j].forward.use = 'k'
	#new_mechanism.elementary_processes[j].backward.use = 'G'
try:
	new_mechanism.elementary_processes[j].backward.use = 'k'
except:
	pass

mechanism_constructor(new_mechanism)

opdh_1 = TAPobject()
opdh_1.mechanism = new_mechanism
opdh_1.reactor_species = new_reactor_species
opdh_1.reactor = new_reactor
opdh_1.show_graph = True

opdh_1.parameter_scale = 120

opdh_1.output_name = 'exp_new'

#forward_problem(1,1,opdh_1)
#flux_graph(opdh_1)
#sys.exit()

#save_object(opdh_1,'./quick_test.json')
##forward_problem(1,1,opdh_1)
#opdh_2 = read_TAPobject('./quick_test.json')
#forward_problem(1,1,opdh_2)

C3H8_int = [0.8,1.2]#[0.0,0.2,0.4,0.6,.0.8,1,1.2,1.4,1.6,1.8,2]
o2_int = [0.8,1.2]#[0.0,0.2,0.4,0.6,.0.8,1,1.2,1.4,1.6,1.8,2]
o2_delay = [0.0,0.2]#[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
temperature = [700,750]#[573,623,673,723,773,823]

for jnum1,j1 in enumerate(C3H8_int):
	for jnum2,j2 in enumerate(o2_int):
		for jnum3,j3 in enumerate(o2_delay):
			for jnum4,j4 in enumerate(temperature):

				try:
					os.mkdir('./'+str(jnum1)+'_'+str(jnum2)+'_'+str(jnum3)+'_'+str(jnum4))
				except OSError:
					pass

				#opdh_1.intensity = j1
				#opdh_1.intensity = j2
				#opdh_1.delay = j3
				#opdh_1.reactor_species.temperature = j4

				save_object(opdh_1,'./'+str(jnum1)+'_'+str(jnum2)+'_'+str(jnum3)+'_'+str(jnum4)+'/mech_1_run.json')

#opdh_1.data_name = './TAP_experimental_data.json'

#OCM_tapobject.gasses_objective = ['']
#OCM_tapobject.optimize = True # , #new_
#OCM_tapobject.parameters_of_interest = ["TAPobject_data.mechanism."]
