import mpmath
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import imageio
import csv
import ast
import shutil
import sys
import os
import scipy
import pip
import pkg_resources
import ufl
import re

from structures import reactor, mechanism, experiment, initial_conditions







class sampling():
	pass
class visualization():

	def __init__(self):
		self.reac = reactor()
		self.mech = mechanism()
		self.iniCon = initial_conditions()
		self.expData = experimentalData()

	#def fluxGraph():
	#def fitGif():
	#def pulseGif():
	#def concentrationDistributionGif():
	#def concentrationDistributionPlot():
	#

class tapsolver():

	def __init__(self):
		self.reac = reactor()
		self.mech = mechanism()
		self.iniCon = initial_conditions()
		self.expData = experimentalData()
		self.processNotes = None

	def generatePDEs(self):
		F = ''

	def confirmConsistency(self):
		pass





	#def simulation():
	#def sensitivity():
	#def fit():
	#def uq():

testGen = initial_conditions()
testGen.readCSVInput('./input_file_2.csv')
testGen.showGasses()
testGen.showSurface()
testGen2 = reactor()
testGen2.readCSVInput('./input_file_2.csv')
testGen3 = mechanism()
testGen3.readCSVInput('./input_file_2.csv')
testGen3.meachanismReactants()
sys.exit()
testGen3.displayProcesses()
testGen4 = experimentalData()
testGen4.readCSVInput('./input_file_2.csv',testGen.gasses)

testGen4 = tapsolver()