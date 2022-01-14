from tapsolver import *
import numpy as np
import os
import jsonpickle
import json
import sys
import copy

#### Step 1
# Names = 2MnSiO2.json , 1.5Mn5Na2WO4SiO2.json , 1Mn5Na2WO4SiO2.json , 2Mn5Na2WO4SiO2.json , 3Mn5Na2WO4SiO2.json , 4Mn5Na2WO4SiO2.json , 5Na2WO4SiO2.json , 5WOxSiO2.json
#f = open('./2MnSiO2.json')
f = open('./1.5Mn5Na2WO4SiO2.json')
try:
	experimental_data = json.load(f)['data']
	conversion = ['time','Ar','O162','O18O16','O182']
except:
	experimental_data = json.load(f)
	conversion = ['Ar','O162','O18O16','O182','time']
f.close()
tapsolver_data = {}

for jnum,j in enumerate(experimental_data.keys()):
	tapsolver_data[conversion[jnum]] = {}
	tapsolver_data[conversion[jnum]][0] = []
	try:
		for k in experimental_data[j]:
			tapsolver_data[conversion[jnum]][0].append(experimental_data[j][str(k)]) 
	except:
		tapsolver_data[conversion[jnum]][0] = experimental_data[j]
save_object(tapsolver_data,'./'+'./'+'/TAP_experimental_data.json')
new_data = read_experimental_data_object('./'+'./'+'/TAP_experimental_data.json')

new_reactor = reactor()

new_reactor.zone_lengths = {0:1.8675,1:0.2,2:1.8675}
new_reactor.zone_voids = {0:0.4,1:0.4,2:0.4}
new_reactor.reactor_radius = 0.19

new_reactor_species = reactor_species()
new_reactor_species.inert_diffusion = 27
new_reactor_species.catalyst_diffusion = 27
new_reactor_species.reference_temperature = 700
new_reactor_species.reference_mass = 29.72#24.01
new_reactor_species.temperature = 700
new_reactor_species.advection = 0

O182 = define_gas()
O182.mass = 36
O182.intensity = 0.118
O182.delay = 0
O182.sigma = 0.02
O182.noise = 0.0
new_reactor_species.add_gas('O182',O182)

O18O16 = define_gas()
O18O16.mass = 34
O18O16.intensity = 0
O18O16.delay = 0.0
O18O16.sigma = 0.0028
O18O16.noise = 0.0
new_reactor_species.add_gas('O18O16',O18O16)

O162 = define_gas()
O162.mass = 32
O162.intensity = 0
O162.delay = 0.0
O162.sigma = 0.0028
O18O16.noise = 0.0
new_reactor_species.add_gas('O162',O162)

Ar = define_gas()
Ar.mass = 40
Ar.intensity = 0.1
Ar.delay = 0.0
Ar.noise = 0.0
Ar.sigma = 0.02
new_reactor_species.add_inert_gas('Ar',Ar)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O182*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O18*',s)

s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O162*',s)

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

new_mechanism.elementary_processes[0] = elementary_process('O182 + * <-> O182*')
new_mechanism.elementary_processes[1] = elementary_process('O182* + * <-> 2O18*')
new_mechanism.elementary_processes[2] = elementary_process('O162 + * <-> O162*')
new_mechanism.elementary_processes[3] = elementary_process('O162* + * <-> 2O16*')
new_mechanism.elementary_processes[4] = elementary_process('O18O16 + * <-> O18O16*')
new_mechanism.elementary_processes[5] = elementary_process('O18O16* + * <-> O18* + O16*')
new_mechanism.elementary_processes[6] = elementary_process('O18* <-> O16*')

new_mechanism.kinetic_links = {0:0.001,1:0.001,2:0.001,3:0.001,4:0.001}

new_mechanism.elementary_processes[0].forward.link = 0
new_mechanism.elementary_processes[0].backward.link = 1
new_mechanism.elementary_processes[1].forward.link = 2
new_mechanism.elementary_processes[1].backward.link = 3
new_mechanism.elementary_processes[2].forward.link = 0
new_mechanism.elementary_processes[2].backward.link = 1
new_mechanism.elementary_processes[3].forward.link = 2
new_mechanism.elementary_processes[3].backward.link = 3
new_mechanism.elementary_processes[4].forward.link = 0
new_mechanism.elementary_processes[4].backward.link = 1
new_mechanism.elementary_processes[5].forward.link = 2
new_mechanism.elementary_processes[5].backward.link = 3
new_mechanism.elementary_processes[6].forward.link = 4
new_mechanism.elementary_processes[6].backward.link = 4


# Make a function
for j in new_mechanism.elementary_processes:
	new_mechanism.elementary_processes[j].forward.use = 'link'
	try:
		new_mechanism.elementary_processes[j].backward.use = 'link'
	except:
		pass

mechanism_constructor(new_mechanism)

OCM_tapobject = TAPobject()
OCM_tapobject.mechanism = new_mechanism
OCM_tapobject.reactor_species = new_reactor_species
OCM_tapobject.reactor = new_reactor

OCM_tapobject.data_name = './TAP_experimental_data.json'
OCM_tapobject.show_graph = False
OCM_tapobject.output_name = 'exp_new'

OCM_tapobject.gasses_objective = ['Ar']
OCM_tapobject.optimize = True # , #new_
OCM_tapobject.parameters_of_interest = ["TAPobject_data.reactor_species.inert_gasses['Ar'].intensity","TAPobject_data.reactor_species.reference_mass"]
#OCM_tapobject.objective_function = True
##20.324451446533203,0.10808474063873305
#def min_function(x):
#        print(x[0])
#        print(x[1])
#        OCM_tapobject.reactor_species.reference_mass = x[0]
#        OCM_tapobject.reactor_species.inert_gasses['Ar'].intensity = x[1]
#        obj_new = forward_problem(0.5,1,OCM_tapobject)
#        print(obj_new)
#        return obj_new
#print('test 3')
#x0 = np.array([23,0.11])
#res = minimize(min_function,x0,method='nelder-mead')
###
#sys.exit()

### Step 2
OCM_tapobject.gasses_objective = ['Ar']
OCM_tapobject.optimize = True # , #new_
OCM_tapobject.parameters_of_interest = ["TAPobject_data.reactor_species.inert_gasses['Ar'].intensity","TAPobject_data.reactor_species.reference_mass"]

### Step 3
#OCM_tapobject.gasses_objective = ['Ar']
#OCM_tapobject.optimize = True
#OCM_tapobject.parameters_of_interest = ['TAPobject_data.reactor_species.reference_mass']

### Step 4
OCM_tapobject.gasses_objective = ['O182','O162','O18O16']#['Ar']
OCM_tapobject.optimize = True
OCM_tapobject.parameters_of_interest = ['TAPobject_data.mechanism.kinetic_links[0]','TAPobject_data.mechanism.kinetic_links[1]','TAPobject_data.mechanism.kinetic_links[2]','TAPobject_data.mechanism.kinetic_links[3]','TAPobject_data.mechanism.kinetic_links[4]']
#OCM_tapobject.parameters_of_interest = ['TAPobject_data.mechanism.elementary_processes[0].forward.k','TAPobject_data.mechanism.elementary_processes[0].backward.k','TAPobject_data.mechanism.elementary_processes[1].forward.k','TAPobject_data.mechanism.elementary_processes[2].forward.k','TAPobject_data.mechanism.elementary_processes[3].forward.k','TAPobject_data.mechanism.elementary_processes[4].forward.k','TAPobject_data.mechanism.elementary_processes[5].forward.k']
#OCM_tapobject.parameters_of_interest = ['TAPobject_data.mechanism.elementary_processes[0].forward.k','TAPobject_data.mechanism.elementary_processes[0].backward.k','TAPobject_data.mechanism.elementary_processes[1].forward.k','TAPobject_data.mechanism.elementary_processes[1].backward.k','TAPobject_data.mechanism.elementary_processes[2].forward.k','TAPobject_data.mechanism.elementary_processes[2].backward.k','TAPobject_data.mechanism.elementary_processes[3].forward.k','TAPobject_data.mechanism.elementary_processes[3].backward.k','TAPobject_data.mechanism.elementary_processes[4].forward.k','TAPobject_data.mechanism.elementary_processes[4].backward.k','TAPobject_data.mechanism.elementary_processes[5].forward.k','TAPobject_data.mechanism.elementary_processes[5].backward.k']


### Step 5
#OCM_tapobject_2 = update_parameters(OCM_tapobject)

forward_problem(0.5,1,OCM_tapobject)
#flux_graph(OCM_tapobject)