# TAPsolver

## Installation

For those new to python, [conda!](https://docs.conda.io/projects/conda/en/latest/use
r-guide/install/) is recommended. Conda makes it easier to install packages and man
age environments than the basic python installation.

	conda create -n fenicsproject -c conda-forge fenics

	source activate fenicsproject

	conda install -c conda-forge dolfin-adjoint

	pip install matplotlib

## Running TAPsolver

With the five core files downloaded in the working directory and with the fenicsproject environment activated, run the following command:

	python tap_sim.py

This will call the simulation process and store all desired information in the appropriate folders. 

An explanation of the input file and how to properly manipulate it can be found [hhere](https://github.com/medford-group/TAPsolver/tree/master/docs/outline/input_file). 

This is the current format of TAPsolver and is subject to change. Steps will be taken to merge this process with 

