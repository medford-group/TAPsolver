
import pandas as pd
import numpy as np

pd.read_excel('tmp.xlsx', index_col=0)


# for each different gas
class gasOptions():
    def __init__(self, *args, **kwargs):
        self.name = ""
        self.amu = 40
        self.gain = 9
        self.baselineStart = 0.9
        self.baselineEnd = 1.0
        self.baselineType = "mean"
        self.calibrationCoef = 1
        self.pumpProbeSpacing = 0
        self.massCorrection = 0
        self.fluxSmoothing = 0
        self.smoothingType = "None"
        self.outliers = "None"
        self.__dict__.update(kwargs)
        availableProperties = ["name", "amu", "gain", "baselineStart", "baselineEnd", "baselineType", "calibrationCoef", "pumpProbeSpacing", "massCorrection", "fluxSmoothing", "smoothingType", "outliers"]
        for i in range(len(args)):
            setattr(self, availableProperties[i], args[i])
    def return_dict(self):
        return(self.__dict__)

class reactorOptions():
    def __init__(self, *args, **kwargs):
        self.inertZone1Length = 0.02
        self.inertZone2Length = 0.02
        self.catalystBedLength = 0.00075
        self.crossSectionalArea = 1.14009E-05
        self.catalystWeight = 1
        self.bedPorosity = 0.4
        self.molPerM0Inert = 1.63E-09
        self.__dict__.update(kwargs)
        availableProperties = ["inertZone1Length", "inertZone2Length", "catalystBedLength", "crossSectionalArea", "catalystWeight", "bedPorosity", "molPerM0Inert"]
        for i in range(len(args)):
            setattr(self, availableProperties[i], args[i])
    def return_dict(self):
        return(self.__dict__)

class experimentOptions():
    def __init__(self, *args, **kwargs):
        self.timeStart = 0
        self.timeEnd = 1
        self.diffusionCoef = 0
        self.inert = "AMU40"
        self.reactant = "AMU40"
        self.products = "AMU40"
        self.rateProcedure = "y"
        self.yProcSmoothing = 0
        self.pulseRemoved = 0
        self.numTime = 100
        self.numCores = 4
        self.numPulses = 5
        self.__dict__.update(kwargs)
        availableProperties = ["timeStart", "timeEnd", "diffusionCoef", "inert", "reactant", "products", "rateProcedure", "yProcSmoothing", "pulseRemoved", "numTime", "numCores", "numPulses"]
        for i in range(len(args)):
            setattr(self, availableProperties[i], args[i])
    def return_dict(self):
        return(self.__dict__)

class simulationOptions():
    def __init__(self, *args, **kwargs):
        self.meshSize = 360
        self.outputFolderName = "results"
        self.experimentalFolderName = "none"
        self.noise = False
        self.theta = 1
        self.solverMethod = "None"
        self.storeOutletFlux = True
        self.storeGraph = False
        self.displayExperimentalData = False
        self.displayGraph = False
        self.sensitivityAnalysis = False
        self.fitParameters = False
        self.optimizationMethod = "CG"
        self.objectivePoints = 1
        self.rrmAnalysis = False
        self.mkmAnalysis = True
        self.petalPlots = False
        self.reactionsTest = ["A + * -> A*"]
        self.massList = [40,40]
        self.pulseTime = [0,0]
        self.pulseRatio = [1,1]
        self.initialSurfaceComposition = [0,1000]
        self.kForward = [200]
        self.kBackward = [0]
        self.__dict__.update(kwargs)
        availableProperties = ["meshSize", "outputFolderName", "experimentalFolderName", "noise", "theta", "solverMethod", "storeOutletFlux", "storeGraph", "displayExperimentalData", "displayGraph", "sensitivityAnalysis", "fitParameters", "optimizationMethod", "objectivePoints", "rrmAnalysis", "mkmAnalysis","petalPlots","reactionsTest", "massList", "pulseTime", "pulseRatio","initialSurfaceComposition","kForward", "kBackward"]
        for i in range(len(args)):
            setattr(self, availableProperties[i], args[i])
    def return_dict(self):
        return(self.__dict__)


class tapParameters():
    def __init__(self, *args, **kwargs):
        self.reactor = reactorOptions()
        self.experiment = experimentOptions()
        self.simulation = simulationOptions()
        self.analysis = ["Start"] # this a list of string to keep track of what analysis the user performed
        self.__dict__.update(kwargs)
        availableProperties = ["reactor", "experiment", "simulation","analysis"]
        for i in range(len(args)):
            setattr(self, availableProperties[i], args[i])
    def return_dict(self):
        self.reactor = self.reactor.return_dict()
        self.experiment = self.experiment.return_dict()
        self.simulation = self.simulation.return_dict()
        return(self.__dict__)

class gasObj():
    def __init__(self, *args, **kwargs):
        self.name = ""
        self.matrix = [] # List per pulse number: each item is a flux
        self.moments = [] # List per pulse number: list within for each different moment type
        self.temperature = [] # Array of temperatures per pulse number
        self.originalBaseline = [] # Array of baseline values of original flux per pulse
        self.timeIndexFlux = [] # timing associated with 
        self.options = gasOptions(self.name)
        self.__dict__.update(kwargs)
        availableProperties = ["name", "matrix", "moments", "temperature", "originalBaseline", "timeIndexFlux"]
        for i in range(len(args)): 
            setattr(self, availableProperties[i], args[i])
    def return_dict(self):
        self.options = self.options.return_dict()
        return(self.__dict__)

class tapObj():
    def __init__(self, gasName):
        self.gasName = gasName
        self.parameters = tapParameters()
        # for each gas create a gasObj
        for name in gasName:
            setattr(self, name, gasObj(name))
    def setSimulation(self):
        self.parameters.simulation.pulseDuration = self.parameters.experiment.timeEnd
        self.parameters.simulation.timeSteps = self.parameters.experiment.numTime
        self.parameters.simulation.reactorLength = self.parameters.reactor.inertZone1Length + self.parameters.reactor.inertZone2Length + self.parameters.reactor.catalystBedLength
        self.parameters.simulation.catalystFraction = self.parameters.reactor.catalystBedLength / self.parameters.simulation.reactorLength
        self.parameters.simulation.reactorRadius = np.sqrt(self.parameters.reactor.crossSectionalArea / np.pi)
        
    def return_dict(self):
        self.parameters = self.parameters.return_dict()
        return(self.__dict__)



from fenics import *
from fenics_adjoint import *
from pyadjoint.enlisting import Enlist
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import imageio
import pandas as pd
import numpy as np
import math
import os
import json

def read_input():
	"""
	Returns a dictionary of simulation parameters from 
	"""
	user_data = pd.read_csv('./input_file.csv',header=None)
	rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
	rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
	rows_4, cols_4 = np.where(user_data == 'Reaction_Information')
	reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = user_data.iloc[1+rows_2[0]:rows_3[0]-1,:]
	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]
	reaction_info = user_data.iloc[1+rows_4[0]:,:]
	reactor_kinetics_input = {}
	for k in range(0,len(reactor_info.index)):
		try:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
	for k in range(0,len(data_storage.index)):
		try:
			reactor_kinetics_input[data_storage.iloc[k,0]] = float(data_storage.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[data_storage.iloc[k,0]] = data_storage.iloc[k,1]
	for k in range(0,len(feed_surf_info.index)):
		try:
			reactor_kinetics_input[feed_surf_info.iloc[k,0]] = float(feed_surf_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[feed_surf_info.iloc[k,0]] = feed_surf_info.iloc[k,1]
	reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
	reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()
	kinetic_parameters = {}
	for j in range(0,len(reaction_info.index)):
		kinetic_parameters['kf'+str(j)] = float(reaction_info.iloc[j,1])
		if str(reaction_info.iloc[j,2]) != 'nan':
			kinetic_parameters['kb'+str(j)] = float(reaction_info.iloc[j,2])
		else:
			pass
	kin_in = kinetic_parameters.copy()
	return reactor_kinetics_input,kinetic_parameters,kin_in

tempData = read_input()

with open("inputFileRoss.txt", "w") as outfile:
    json.dump(tempData, outfile)

def readNew():
    with open("inputFileRoss.txt") as jsonFile:
        data = json.load(jsonFile)
    reactor_kinetics_input = data[0]
    kinetic_parameters = data[1]
    kin_in = data[2]
    reactor_kinetics_input["Experimental Data Folder"] = np.nan
    #if type(reactor_kinetics_input['reactions_test']) != "list" :
    #    reactor_kinetics_input['reactions_test'] = [reactor_kinetics_input['reactions_test']]
    return reactor_kinetics_input,kinetic_parameters,kin_in