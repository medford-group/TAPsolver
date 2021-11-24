
* [Overview](https://github.com/medford-group/TAPsolver/tree/master)
* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Documentation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/Documentation)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion)

# TAPsolver Crash Course

An example of how to run TAPsolver (with predefined input files) is presented here. Approaches for varying the input files will be added later. The two input files necessary to run these TAPsolver functions are provided in this folder (i.e. input_file_exp.csv and input_file_fit.csv) and should be downloaded to the directory where the processes will be run.

## Getting started and running a simulation (the forward problem):

All available functions written in TAPsolver are accesable through the importing of the tapsolver package as follows:

	from tapsolver import *

With the two input files in the directory where this script is being written, a basic forward simulation can be run with:

	run_tapsolver(1,input_file='./input_file_exp.csv',add_noise=True,pulseNumber=1)

where the "1" represents the amount of time to simulate each pulse for (in units of seconds), "input_file" is the file you wish to use to run the simulation, "add_noise" specifies if you would like to include normally distrubuted noise around the outlet flux curves, and the "pulseNumber" details how many pulses the user would like to include in the simulation.

So, the user can exclude noise by changing "add_noise" to False and run a multi-pulse (state-altering) experiment by changing "pulseNumber" to 10. To visualize the outlet fluxes, the following command can be utilized:

flux_graph(input_file = './input_file_exp.csv',pulse=[1,2,3,4,5],output_name='./flux.png')

where "pulse" represents which pulses you wish to include in the visual (if variable is excluded in function call, it defaults to the first pulse) and "output_name" is the file name you wish to save the figure as (if variable is excluded in function call, it defaults to './flux.png'). Some additional parameters to adjust what you're viewing include:

	flux_graph(input_file = './input_file_exp.csv',dispAnalytic=True)

which plots the analytical solution for pure diffusion (the gas species were not reacting at all).

	flux_graph(input_file = './input_file_exp.csv',store_graph=False)

which means no plot will be saved upon generation.

	flux_graph(input_file = './input_file_exp.csv',dispExper=False)

which shows the experimental data (if it specified in the input file).

You can also generate a gif showing the pulsing process (instead of viewing all pulses simultaneously) with the command:

	pulse_gif(input_file = './input_file_exp.csv',output_name = './output.gif')

## Analyzing the data (the inverse problem)

The input_file "./input_file_exp.csv" is used to generate our synthetic set of data. We assume that this is a black box where we can manipulate the inputs (temperature, intensity, delay, etc.) and get out outlet fluxes (with noise if desired). In most circumstances, we will not be able to observe the microkinetic model specified in the input_file and it is our job to determine what it is.

With this in mind, we write the additional input_file "./input_file_fit.csv" with some proposed mechanism (in which we have no structural uncertainty whatsoever, a circumstance that is likely to never happen :) ). To fit the proposed model to the experimental data, we can run the command: 

	fit_tap(0.1,sigma=0.1,input_file = './input_file_fit.csv')

with the commands '0.1' representing the window of time to fit the parameter over (0.1 seconds), 'sigma' being the assumed noise distribution (1 S.D. = 0.1 nmol/s) and the input_file used with the propsed mechanism included. 

An additionally useful method for analyzing the proposed model is to quantify the uncertainty through a Hessian-analysis (i.e. looking at the second order derivatives at the local minimum to determine the confidence in your kinetic parameters). This is accomplished through the command:

	run_uncertainty(0.06,sigma=0.1,input_file = './input_file_fit.csv')

