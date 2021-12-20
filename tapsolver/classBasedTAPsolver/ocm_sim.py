from fenics import *
from fenics_adjoint import *
from structures import *
from reactor_species import *
from file_io import *
from mechanism_construction import *
from visualization import *
from forward_problem import *
import jsonpickle
import json
import sys
import copy

f = open('./2MnSiO2.json')
experimental_data = json.load(f)
f.close()
tapsolver_data = {}

conversion = ['Ar','O162','O18O16','O182','time']
for jnum,j in enumerate(experimental_data.keys()):
	tapsolver_data[conversion[jnum]] = {}
	tapsolver_data[conversion[jnum]][0] = []
	for k in experimental_data[j]:
		
		tapsolver_data[conversion[jnum]][0].append(experimental_data[j][str(k)]) 
save_object(tapsolver_data,'./'+'./'+'/TAP_experimental_data.json')
new_data = read_experimental_data_object('./'+'./'+'/TAP_experimental_data.json')


#with open('./2MnSiO2.json', 'r', encoding='utf-8') as f:
#	dumped2 = f
#f.close()
#print(dumped2)

#sameObject = jsonpickle.decode(dumped2)
#print(sameObject.reactor_species.gasses)

new_reactor = reactor()

new_reactor.zone_lengths = {0:1.8675,1:0.07664,2:1.8675}
new_reactor.zone_voids = {0:0.4,1:0.4,2:0.4}
new_reactor.reactor_radius = 0.19

### Add the reactor gasses (reactive and inert)

new_reactor_species = reactor_species()
new_reactor_species.inert_diffusion = 14.5
new_reactor_species.catalyst_diffusion = 14.5
new_reactor_species.reference_temperature = 385.5
new_reactor_species.reference_mass = 40.1
new_reactor_species.temperature = 973
new_reactor_species.advection = 0
new_reactor_species.reference_pulse_size = 1

O182 = define_gas()
O182.mass = 36
O182.intensity = 0.16
O182.delay = 0
O182.sigma = 0.02
O182.noise = 0.0
new_reactor_species.add_gas('O182',O182)

O18O16 = define_gas()
O18O16.mass = 34
O18O16.intensity = 0.16
O18O16.delay = 0.0
O18O16.sigma = 0.0028
O18O16.noise = 0.0
new_reactor_species.add_gas('O18O16',O18O16)

O162 = define_gas()
O162.mass = 32
O162.intensity = 0.16
O162.delay = 0.0
O162.sigma = 0.0028
O18O16.noise = 0.0
new_reactor_species.add_gas('O162',O162)

Ar = define_gas()
Ar.mass = 40
Ar.intensity = 0.16
Ar.delay = 0.0
Ar.noise = 0.0
Ar.sigma = 0.02
new_reactor_species.add_inert_gas('Ar',Ar)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O18*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O16*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O18O16*',s)

s = define_adspecies()
s.concentration = 120
new_reactor_species.add_adspecies('*',s)

new_mechanism = mechanism()

new_mechanism.elementary_processes[0] = elementary_process('O182 + 2* <-> 2O18*')
new_mechanism.elementary_processes[1] = elementary_process('O18* <-> O16*')
new_mechanism.elementary_processes[2] = elementary_process('O18* + O16* <-> O18O16* + *')
new_mechanism.elementary_processes[3] = elementary_process('O18O16* <-> O18O16 + *')
new_mechanism.elementary_processes[4] = elementary_process('2O16* <-> O162 + 2*')

new_mechanism.elementary_processes[0].forward.k = 0
new_mechanism.elementary_processes[0].backward.k = 0
new_mechanism.elementary_processes[1].forward.k = 0
new_mechanism.elementary_processes[1].backward.k = 0
new_mechanism.elementary_processes[2].forward.k = 0
new_mechanism.elementary_processes[2].backward.k = 0
new_mechanism.elementary_processes[3].forward.k = 0
new_mechanism.elementary_processes[3].backward.k = 0
new_mechanism.elementary_processes[4].forward.k = 0
new_mechanism.elementary_processes[4].backward.k = 0

mechanism_constructor(new_mechanism)

OCM_tapobject = TAPobject()
OCM_tapobject.mechanism = new_mechanism
OCM_tapobject.reactor_species = new_reactor_species
OCM_tapobject.reactor = new_reactor

OCM_tapobject.data_name = './TAP_experimental_data.json'
OCM_tapobject.show_graph = False
OCM_tapobject.output_name = 'exp_new'
OCM_tapobject.gasses_objective = ['Ar']#['O182','O162','O18O16']
OCM_tapobject.optimize = True
OCM_tapobject.parameters_of_interest = ["TAPobject_data.reactor_species.reference_mass","TAPobject_data.reactor_species.inert_gasses['Ar'].intensity"]##['TAPobject_data.mechanism.elementary_processes[0].forward.k']
forward_problem(0.5,1,OCM_tapobject)
flux_graph(OCM_tapobject)
sys.exit()