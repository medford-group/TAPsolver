from tapsolver_working import *
#import libmuqModelling
#from catmap import *

#import tapsolver
#from tap_sim import *
#from reac_odes import *
#from fenics import *
#import os
#from shutil import copyfile
#import sys

# Defining the input file
#define_reactor(reactor_name='./reactorInl.csv')

#'./o2_mech_1.csv'
#'./o2_mech_2.csv'
#'./o2_mech_3.csv'

#run_batch_reactor(timeFunc=2,includeNoise=False,inputFile='./input_file_batch.csv')
#batchGraph(input_file='./input_file_batch.csv',dispExper='none')
#fit_batch(timeFunc=2,inputFile='./input_file_batch.csv')
#batch_fit()

#sys.exit()
run_tapsolver(timeFunc=1,includeNoise=True,store_thin_func=False,inputFile='./input_file_thermoConst.csv',pulseNumber=1,input_form='new')
#run_tapsolver(timeFunc=0.3,includeNoise=True,store_thin_func=False,inputFile='./input_file_thermoConst.csv',pulseNumber=1,input_form='new')
#fit_parameters(0.1,inputFile='./input_file_thermoConst.csv',input_form='new')

fluxGraph(input_file='./input_file_thermoConst.csv',analytical=True,inputForm='new',pulse=1)
#fluxGraph(input_file='./input_file_thermoConst.csv',dispExper='none',objectivePoints=True,analytical=False,inputForm='new',pulse=1)
#uncert_quant(0.1,input_file='./input_file_thermoConst.csv',inputForm='new')

#design_experiment(0.1,inputFile = './input_file_thermoConst.csv',input_form='new',altVariables=['pulses'])

#input_construction(reactor_name='./o2_reactor.csv',reaction_name='./o2_mech_1.csv',new_name='./input_file_o2m1.csv')
#vary_Input('Output Folder Name','test',input_file='./input_file.csv')

# Running Forward Problem
#sys.exit()
#design_experiment(0.4,inputFile = './input_file_thermoConst.csv',altVariables=['kf2']) # ['pulses','time']
#run_tapsolver(timeFunc=0.4,store_thin_func=False,inputFile='./input_file_o2m1.csv',pulseNumber=1,input_form='new')

# Algorithmic Differentation Features
#run_sensitivity(0.1,'total')
#fitting_gif()
#uncert_quant(0.1)

# Visualization
#fluxGraph(input_file='./input_file_thermoConst.csv',inputForm='new')
#fluxGraph(input_file='./input_file_o2m1.csv',dispExper='none',analytical=False,inputForm='new',pulse=1)
#pulsingVisualization(input_file = './input_file_inl.csv',fileName = './output.gif',inputForm='new')
#concDistPlot(input_file='./input_file.csv',input_form = 'new')
