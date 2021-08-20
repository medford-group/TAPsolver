from structures import *
from file_io import *
import sys
from mechanism_construction import *
from initial_conditions import *

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


testGen = readCSV_ic('./input_file_2.csv')
testGen2 = readCSV_reactor('./input_file_2.csv')
testGen3 = readCSV_mechanism(testGen3,'./input_file_2.csv')

testGen3 = mechanism_constructor(testGen3)
f_new, domain = construct_f_equation(testGen3,testGen,testGen2)
print(f_new)
sys.exit()
testGen3.displayProcesses()
testGen4 = experimentalData()
testGen4.readCSVInput('./input_file_2.csv',testGen.gasses)

testGen4 = tapsolver()