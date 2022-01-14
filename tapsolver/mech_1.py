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
Ar.mass = 40.1
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
new_reactor_species.add_adspecies('O*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('C3H6*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('H*',s)

s = define_adspecies()
s.concentration = 120
new_reactor_species.add_adspecies('*',s)

new_mechanism = mechanism()

new_mechanism.elementary_processes[0] = elementary_process('C3H8 + * <-> C3H8*') # Adsorption
new_mechanism.elementary_processes[1] = elementary_process('O2 + 2* <-> 2O*') # Adsorption
new_mechanism.elementary_processes[2] = elementary_process('C3H6 + * <-> C3H6*') # Adsorption
#new_mechanism.elementary_processes[3] = elementary_process('H2O + * <-> H2O*') # Adsorption
#new_mechanism.elementary_processes[4] = elementary_process('CO2 + * <-> CO2*') # Adsorption
new_mechanism.elementary_processes[3] = elementary_process('C3H8* + O* <-> C3H6* + H2O')
#new_mechanism.elementary_processes[5] = elementary_process('O* + 2H* <-> H2O*')
new_mechanism.elementary_processes[4] = elementary_process('C3H8* + 9O* <-> 3CO2 + 4H2O')
new_mechanism.elementary_processes[5] = elementary_process('C3H6* + 9O* <-> 3CO2 + 3H2O')

# Random
###new_mechanism.elementary_processes[0].forward.dG = -0.75
###new_mechanism.elementary_processes[0].forward.Ga = 1.45
###new_mechanism.elementary_processes[1].forward.dG = -1.5
###new_mechanism.elementary_processes[1].forward.Ga = 1.35
###new_mechanism.elementary_processes[2].forward.dG = -.35
###new_mechanism.elementary_processes[2].forward.Ga = 1.1
####new_mechanism.elementary_processes[3].forward.dG = -.4
####new_mechanism.elementary_processes[3].forward.Ga = 0.8
####new_mechanism.elementary_processes[4].forward.dG = -.4
####new_mechanism.elementary_processes[4].forward.Ga = 1
###new_mechanism.elementary_processes[3].forward.dG = -.45
###new_mechanism.elementary_processes[3].forward.Ga = 0.7
####new_mechanism.elementary_processes[5].forward.dG = -.38
####new_mechanism.elementary_processes[5].forward.Ga = 0.35
###new_mechanism.elementary_processes[4].forward.dG = -.8
###new_mechanism.elementary_processes[4].forward.Ga = 0.5
###new_mechanism.elementary_processes[5].forward.dG = -1.5
###new_mechanism.elementary_processes[5].forward.Ga = 0.4

# Free energy figure
new_mechanism.elementary_processes[0].forward.dG = -0.35
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

#new_mechanism.elementary_processes[0].forward.dG = -0.75
#new_mechanism.elementary_processes[0].forward.Ga = 10
#new_mechanism.elementary_processes[1].forward.dG = -1.5
#new_mechanism.elementary_processes[1].forward.Ga = 10
#new_mechanism.elementary_processes[2].forward.dG = -.35
#new_mechanism.elementary_processes[2].forward.Ga = 10
#new_mechanism.elementary_processes[3].forward.dG = -.61
#new_mechanism.elementary_processes[3].forward.Ga = 10#0.56
#new_mechanism.elementary_processes[4].forward.dG = -.89
#new_mechanism.elementary_processes[4].forward.Ga = 10 #1.29
#new_mechanism.elementary_processes[5].forward.dG = -0.89
#new_mechanism.elementary_processes[5].forward.Ga = 10


## Make a function
for j in new_mechanism.elementary_processes:
	new_mechanism.elementary_processes[j].forward.use = 'G'
try:
	new_mechanism.elementary_processes[j].backward.use = 'G'
except:
	pass

mechanism_constructor(new_mechanism)

opdh_1 = TAPobject()
opdh_1.reactor = new_reactor
opdh_1.reactor_species = new_reactor_species
opdh_1.mechanism = new_mechanism

opdh_1.show_graph = False
opdh_1.parameter_scale = 120
opdh_1.output_name = 'exp_new'

#temp_gasses = opdh_1.reactor_species.gasses
#temp_inert_gasses = opdh_1.reactor_species.inert_gasses
#temp_adspecies = opdh_1.reactor_species.adspecies
#temp_reactants = opdh_1.mechanism.reactants
##opdh_2 = new_experiments(opdh_1,'-B',original_gasses=temp_gasses,original_inert_gasses=temp_inert_gasses,original_adspecies=temp_adspecies)
##opdh_3 = new_experiments(opdh_2,'-C',original_gasses=temp_gasses,original_inert_gasses=temp_inert_gasses,original_adspecies=temp_adspecies)
##opdh_2.reactor_species.temperature = {0:700,1:700,2:700}
#opdh_2.reactor_species.gasses['C3H8'].intensity = 1
#opdh_2.reactor_species.gasses['O2'].intensity = 1
#opdh_2.reactor_species.gasses['C3H8-B'].intensity = 1
#opdh_2.reactor_species.gasses['O2-B'].intensity = 1
#mechanism_constructor(new_mechanism)
#print(opdh_3.reactor_species.gasses)
#print(opdh_3.reactor_species.inert_gasses)
#display_surface(opdh_2.reactor_species)
##print(opdh_3.reactor_species.reactants)
#display_gasses(opdh_2.reactor_species)
#sys.exit()
#display_surface(opdh_3.reactor_species)
#display_elementary_processes(opdh_3.mechanism)
#sys.exit()
#2.41 min
#8.45 min
#18.139 min
opdh_1.parameters_of_interest = ['TAPobject_data.mechanism.elementary_processes[0].forward.dG'],#'TAPobject_data.mechanism.elementary_processes[1].forward.dG','TAPobject_data.mechanism.elementary_processes[2].forward.dG','TAPobject_data.mechanism.elementary_processes[3].forward.dG','TAPobject_data.mechanism.elementary_processes[4].forward.dG','TAPobject_data.mechanism.elementary_processes[5].forward.dG','TAPobject_data.mechanism.elementary_processes[0].forward.Ga','TAPobject_data.mechanism.elementary_processes[1].forward.Ga','TAPobject_data.mechanism.elementary_processes[2].forward.Ga','TAPobject_data.mechanism.elementary_processes[3].forward.Ga','TAPobject_data.mechanism.elementary_processes[4].forward.Ga','TAPobject_data.mechanism.elementary_processes[5].forward.Ga']

#opdh_1.parameters_of_interest = ['TAPobject_data.mechanism.elementary_processes[0].forward.Ga','TAPobject_data.mechanism.elementary_processes[1].forward.Ga','TAPobject_data.mechanism.elementary_processes[2].forward.Ga','TAPobject_data.mechanism.elementary_processes[3].forward.Ga','TAPobject_data.mechanism.elementary_processes[4].forward.Ga','TAPobject_data.mechanism.elementary_processes[5].forward.Ga']
##opdh_1.finite_difference_trans_sensitivty = True
#finite_difference_trans_sensitivty 
##transient_sensitivity(0.5, 1, opdh_1)
##sys.exit()
#forward_problem(0.5,1,opdh_1)
opdh_1.data_name = 'exp_new/TAP_experimental_data.json'
opdh_1.gasses_objective = ['C3H8','O2','C3H6','H2O','CO2']
opdh_1.output_name = 'fit_new'
opdh_1.optimize = True
#forward_problem(0.5,1,opdh_1)
opdh_1.optimize = True
opdh_1.uncertainty = True
forward_problem(0.1,1,opdh_1)
#flux_graph(opdh_1)
sys.exit()
C3H8_int = [0,0.5,1,1.5,2]
o2_int = [0,0.5,1,1.5,2]
o2_delay = [0,0.5,1]
temperature = [650,700,750]

for jnum1,j1 in enumerate(C3H8_int):
	for jnum2,j2 in enumerate(o2_int):
		for jnum3,j3 in enumerate(o2_delay):
			for jnum4,j4 in enumerate(temperature):

				try:
					os.mkdir('./'+str(jnum1)+'_'+str(jnum2)+'_'+str(jnum3)+'_'+str(jnum4))
				except OSError:
					pass

				opdh_1.reactor_species.gasses['C3H8'].intensity = j1
				opdh_1.reactor_species.gasses['O2'].intensity = j2
				opdh_1.reactor_species.gasses['O2'].delay = j3
				opdh_1.reactor_species.temperature = j4

				save_object(opdh_1,'./'+str(jnum1)+'_'+str(jnum2)+'_'+str(jnum3)+'_'+str(jnum4)+'/mech_1.json')

#opdh_1.data_name = './TAP_experimental_data.json'

#OCM_tapobject.gasses_objective = ['']
#OCM_tapobject.optimize = True # , #new_
#OCM_tapobject.parameters_of_interest = ["TAPobject_data.mechanism."]
