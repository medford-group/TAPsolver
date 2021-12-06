from fenics import *
from fenics_adjoint import *
from structures import *
from reactor_species import *
from file_io import *
from mechanism_construction import *
from forward_problem import *
import jsonpickle
import json
import sys
import copy

# Add data storage
# Add thin zone storage
# Add data read
# Add object storage
# Add objective function reading
# Add control specifications
# Read / make new objects

# 1. Edit TAP object
# 2. Store in the new folder
# 3. Edit TAP object

testGen1 = reactor()
#testGen1 = readCSV_reactor('./input_file_2.csv')
new_reactor_species = reactor_species()
new_reactor_species.advection = 0

CO = define_gas()
CO.mass = 36
CO.intensity = 1
CO.delay = 0.0
new_reactor_species.add_gas('CO',CO)

O2 = define_gas()
O2.mass = 36
O2.intensity = 1
O2.delay = 0.0
new_reactor_species.add_gas('O2',O2)

CO2 = define_gas()
CO2.mass = 36
CO2.intensity = 1
CO2.delay = 0.0
new_reactor_species.add_gas('CO2',CO2)

argon = define_gas()
argon.mass = 36
argon.intensity = 1
argon.delay = 0.0
new_reactor_species.add_inert_gas('argon',argon)

helium = define_gas()
helium.mass = 36
helium.intensity = 1
helium.delay = 0.0
new_reactor_species.add_inert_gas('helium',helium)


#jsonStr = json.dumps(vars(new_reactor_species))

#sys.exit()

s = define_adspecies()
s.concentration = 30
new_reactor_species.add_adspecies('CO*',s)

s = define_adspecies()
s.concentration = 20
new_reactor_species.add_adspecies('O*',s)
new_mechanism = mechanism()

s = define_adspecies()
s.concentration = 10
new_reactor_species.add_adspecies('*',s)

new_mechanism.elementary_processes[0] = elementary_process('CO + * <-> CO*')
new_mechanism.elementary_processes[1] = elementary_process('O2 + 2* <-> 2O*')
new_mechanism.elementary_processes[2] = elementary_process('CO* + O* <-> 2* + CO2')

new_mechanism.elementary_processes[0].forward.k = 1
new_mechanism.elementary_processes[0].backward.k = 1
new_mechanism.elementary_processes[1].forward.k = 1
new_mechanism.elementary_processes[1].backward.k = 1
new_mechanism.elementary_processes[2].forward.k = 1
new_mechanism.elementary_processes[2].backward.k = 1


mechanism_constructor(new_mechanism)


TAP_test = TAPobject()
TAP_test.mechanism = new_mechanism
TAP_test.reactor_species = new_reactor_species
TAP_test.reactor = testGen1

save_object(TAP_test,'./TAP_test.json')
TAP_test = read_TAPobject('./TAP_test.json')
#sys.exit()

#dumped = jsonpickle.encode({1: testGen1})
#save_object(testGen1,'./reactor_gen_test.json')
#save_object(new_mechanism,'./new_mechanism_test.json')
#save_object(new_reactor_species,'./reactor_species_test.json')
#load_reactor = read_reactor_object('./reactor_gen_test.json')
#
#load_species = read_reactor_species_object('./reactor_species.json')
#load_mechanism = read_mechanism_object('./mechanism.json')

#TAP_test = TAPobject()
#TAP_test.mechanism = load_mechanism 
#TAP_test.reactor_species = load_species
#TAP_test.reactor = load_reactor

#helium.intensity = Constant(helium.intensity)
#helium.intensity = Control(helium.intensity)
#TAP_test_2 = new_experiments(TAP_test,'-B')
#print(display_gasses(TAP_test_2.reactor_species))
#print(TAP_test_2)
TAP_test.gasses_objective = ['CO','O2','CO2']
list_of_variables = ['TAPobject_data.reactor_species.gasses["CO"].mass','TAPobject_data.reactor_species.temperature']
#TAP_test_2 = TAP_test.copy()
#sys.exit()
#TAP_test_2.reactor_species.temperature[1] = 350

#print(TAP_test_2.reactor_species.gasses['CO'].catalyst_diffusion)
#print(TAP_test_2.reactor_species.gasses['CO-B'].catalyst_diffusion)

#print(TAP_test_2.reactor_species.inert_gasses['helium'].inert_diffusion)
#print(TAP_test_2.reactor_species.inert_gasses['helium-B'].inert_diffusion)

#sys.exit()
#TAP_test_2.gasses_objective = ['CO','O2','CO2']

for k in range(2,5):
	TAP_test.parameters_of_interest = ['TAPobject_data.mechanism.elementary_processes[0].forward.k']
	TAP_test.output_name = './exp_new/flux_data_'+str(k)+'.json'
	TAP_test.reactor_species.temperature = k*365.5#{0:365.5,1:k*365.5}
	forward_problem(0.02,1,TAP_test)
	sys.exit()

sys.exit()

#with open('./reactor.json', 'r', encoding='utf-8') as f:
#	dumped2 = f
#f.close()

sys.exit()

sys.exit()
print(dumped)
sameObject = jsonpickle.decode(dumped)['1']
print(sameObject.reactor_species.gasses)




def json_load_file(filename, use_jsonpickle = True):
	f = open(filename)
json_str = f.read()

TAP_test.mesh = 200
TAP_test.catalyst_mesh_density = 5
forward_problem(1,1,TAP_test)
sys.exit()
print(display_gasses(new_reactor_species))



import pandas as pd

from mechanism_construction import *
from forward_problem import *
from initial_conditions import *
from reference_files import *


import os
import pickle

class sampling():
	pass

class visualization():

	def __init__(self):
		self.reac = reactor()
		self.mech = mechanism()
		self.iniCon = initial_conditions()
		self.expData = experimentalData()

class tapsolver():

	def __init__(self):
		self.reac = reactor()
		self.mech = mechanism()
		self.iniCon = initial_conditions()
		self.expData = experimentalData()
		self.processNotes = None

#testGen1 = load_example('inl_reactor_1.obj')
testGen1 = readCSV_reactor('./input_file_2.csv')
testGen2 = readCSV_mechanism('./input_file_2.csv')
testGen3 = readCSV_ic('./input_file_2.csv')

print(testGen1.cross_sectional_area())
print(testGen1.reactor_radius)
testGen1.reactor_radius = 2
print(testGen1.cross_sectional_area())
print(testGen1.reactor_radius)
#forward_problem(1,testGen1,testGen2,testGen3)

# elementary_process, elementary_process_details
sys.exit()

testGen_test = load_example('parallel_1.obj')
display_elementary_processes(testGen_test)

new_reactor = pickle.dump(testGen3,open('./parallel_2.obj','wb'))
new_reactor2 = pickle.load(open('./parallel_2.obj','rb'))
print(new_reactor2)

testGen3 = mechanism_constructor(testGen3)
f_new, domain = construct_f_equation(testGen3,testGen,testGen2)
display_elementary_processes(testGen3)

testGen3.displayProcesses()
testGen4 = experimentalData()
testGen4.readCSVInput('./input_file_2.csv',testGen.gasses)

testGen4 = tapsolver()
print('')
print('')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Predefined reactors available include:')
print(' - INL_reactor_1.obj ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Predefined mechanisms available include:')
print(' - irreversible_adsorption.obj')
print(' - reversible_adsorption.obj')
print(' - irreversible_dissociative_adsorption.obj')
print(' - reversible_dissociative_adsorption.obj')
print(' - series_1.obj')
print(' - series_2.obj')
print(' - series_3.obj')
print(' - parallel_1.obj')
print(' - parallel_2.obj')
print(' - multi_site_1.obj')
print(' - multi_site_2.obj')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('')
print('')
