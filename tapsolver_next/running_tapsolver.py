from tap_sim import *
from reac_odes import *
import os
from shutil import copyfile
import sys

# Defining the input file
define_reactor(reactor_name='./inlTap_1.csv')
input_construction(reactor_name='./inlTap_1.csv',reaction_name='./reactionList2.csv',new_name='./input_file_inl.csv')
vary_Input('Output Folder Name','test',input_file='./input_file.csv')

# Running Forward Problem
run_tapsolver(timeFunc=1,inputFile='./input_file_inl.csv',pulseNumber=100,input_form='new')

# Algorithmic Differentation Features
run_sensitivity(0.1,'total')
fitting_gif()
uncert_quant(0.1)

# Visualization
fluxGraph(1,input_file='./input_file_inl.csv',inputForm='new',pulse=list(range(95,100)))
pulsingVisualization(input_file = './input_file_inl.csv',fileName = './output.gif',inputForm='new')
concnDistPlot()
