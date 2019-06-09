# TAPsolver

## Installing Packages

TAPsolver is based around the python package FEniCS. Installation of this package can be challenging and it is recommended to install it on a linux operating system. To install all the necessary packages for using TAPsolver, enter the following commands in your terminal:

conda create -n fenicsproject -c conda-forge fenics

source activate fenicsproject

conda install -c conda-forge dolfin-adjoint

pip install matplotlib

For those new to python, [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) is recommended. Conda makes it easier to install packages and manage environments than the basic python installation.

## File structure

The current version of [TAPsolver](https://github.com/medford-group/TAPsolver/tree/master/tapsolver/v1) consists of five files, including tap_sim.py, func_sim.py, reac_odes.py, vari_form.py and input_file.csv. How these files are related to one another are shown in Figure 1.



When the simulator is run, an output directory will be generated with all the relevent information stored in appropriate directories. 

Two forms of TAPsolver currently exist. One for FEniCS version 2017 and the other for FEniCS version 2018. 

## User Input File Structure

A csv file is used as the input to the simulator. Four main componenents of TAPSolver input are described and include "Reactor Information", "Feed and Surface Composition", "Data Storage Options", and "Reaction Information". The terms in each of these components are discussed below. 

# Reactor Information



# Feed and Surface Composition

# Data Storage Options

# Reaction Information

## Examples

# Simulating Data

# Fitting Data
