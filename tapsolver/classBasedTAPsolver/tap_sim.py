from structures import *
from file_io import *
import sys
import pandas as pd

from mechanism_construction import *
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

testGen2 = load_example('inl_reactor_1.obj')
examples()
testGen2 = readCSV_reactor('./input_file_2.csv')
testGen3 = readCSV_mechanism('./input_file_2.csv')
testGen = readCSV_ic('./input_file_2.csv')

testGen_test = load_example('parallel_1.obj')
display_elementary_processes(testGen_test)
#sys.exit()

new_reactor = pickle.dump(testGen3,open('./parallel_2.obj','wb'))
new_reactor2 = pickle.load(open('./parallel_2.obj','rb'))
print(new_reactor2)
sys.exit()

testGen3 = mechanism_constructor(testGen3)
f_new, domain = construct_f_equation(testGen3,testGen,testGen2)
display_elementary_processes(testGen3)
sys.exit()
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