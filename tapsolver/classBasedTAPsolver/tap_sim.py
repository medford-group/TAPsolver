from structures import *
from reactor_species import *
from file_io import *
from mechanism_construction import *
from forward_problem import *
import json
import sys

testGen1 = readCSV_reactor('./input_file_2.csv')
new_reactor_species = reactor_species()
new_reactor_species.advection = 0
CO = define_gas()
CO.mass = 28
CO.intensity = 3
CO.delay = 0.1
new_reactor_species.add_gas('CO',CO)
O2 = define_gas()
O2.mass = 36
O2.intensity = 1
O2.delay = 0.3
new_reactor_species.add_gas('O2',O2)
CO2 = define_gas()
CO2.mass = 36
CO2.intensity = 1
CO2.delay = 0.3
new_reactor_species.add_gas('CO2',CO2)
new_reactor_species.gasses['CO'].mass = 2
#new_reactor_species.reference_temeprature = 450
#print(display_gasses(new_reactor_species))
argon = define_gas()
argon.mass = 36
argon.intensity = 1
argon.delay = 0.3
new_reactor_species.add_inert_gas('argon',argon)
s = define_adspecies()
s.concentration = 100
new_reactor_species.add_adspecies('*',s)
s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('CO*',s)
s = define_adspecies()
s.concentration = 0
new_reactor_species.add_adspecies('O*',s)
new_mechanism = mechanism()

new_mechanism.elementary_processes[0] = elementary_process()			
new_mechanism.elementary_processes[0].processString = '* + CO <-> CO*'
new_mechanism.elementary_processes[0].forward = elementary_process_details()
new_mechanism.elementary_processes[0].backward = elementary_process_details()

new_mechanism.elementary_processes[1] = elementary_process()			
new_mechanism.elementary_processes[1].processString = '2* + O2 <-> 2O*'
new_mechanism.elementary_processes[1].forward = elementary_process_details()
new_mechanism.elementary_processes[1].backward = elementary_process_details()

new_mechanism.elementary_processes[2] = elementary_process()			
new_mechanism.elementary_processes[2].processString = 'O* + CO* <-> 2* + CO2'
new_mechanism.elementary_processes[2].forward = elementary_process_details()
new_mechanism.elementary_processes[2].backward = elementary_process_details()

mechanism_constructor(new_mechanism)
new_mechanism.elementary_processes[0].forward.k["value"] = 4
new_mechanism.elementary_processes[0].backward.k["value"] = 1
new_mechanism.elementary_processes[1].forward.k["value"] = 4
new_mechanism.elementary_processes[1].backward.k["value"] = 1
new_mechanism.elementary_processes[2].forward.k["value"] = 4
new_mechanism.elementary_processes[2].backward.k["value"] = 1
display_elementary_processes(new_mechanism)

TAP_test = TAPobject()
TAP_test.mechanism = new_mechanism
TAP_test.reactor_species = new_reactor_species
TAP_test.reactor = testGen1
forward_problem(1,TAP_test)
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
forward_problem(1,testGen1,testGen2,testGen3)

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