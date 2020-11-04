from tap_sim import *
import sys

#vary_Input('Display Graph','TRUE',input_file='./input_file.csv')
#vary_Input('Output Folder Name','test',input_file='./input_file.csv')
#vary_Input('Time Steps',300,input_file='./input_file.csv')


run_tapsolver(0.3,input_file='./input_file.csv')


print('done')
sys.exit()

#perform_simulation(storeThin='sparse' or 'robust'):
#
#generate_graph(display = False, analytical = False, scale = False, time_scale = None):
#
#catalyst_zone_gif()
#
#fit_parameters(method='L-BFGS-B',):
#
#fit_inert():
#
#transient_sens():
#
#sensitivity_analysis():
#
#design_of_experiment():
#
#fitting_gif(objective = False, expData = False, ):
