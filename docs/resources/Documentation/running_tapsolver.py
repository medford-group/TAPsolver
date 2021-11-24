from tapsolver import *

run_tapsolver(1,input_file='./input_file_exp.csv',add_noise=True,pulseNumber=1)
run_tapsolver(1,input_file='./input_file_fit.csv',add_noise=False)
flux_graph(input_file = './input_file_fit.csv',pulse=None,dispExper=True,dispAnalytic=False,dispObjective=False,show_graph=True,store_graph=False,output_name='./flux.png')

run_sensitivity(0.1,sigma=0.1,input_file = './input_file_fit.csv')
fit_tap(0.1,sigma=0.1,input_file = './input_file_fit.csv')
run_uncertainty(0.06,sigma=0.1,input_file = './input_file_fit.csv')
