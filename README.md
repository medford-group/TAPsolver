# TAPsolver

## Installing Packages

TAPsolver is based around the python package FEniCS. Installation of this package can be challenging and it is recommended to install it on a linux operating system. To install all the necessary packages for using TAPsolver, enter the following commands in your terminal:

conda create -n fenicsproject -c conda-forge fenics

source activate fenicsproject

conda install -c conda-forge dolfin-adjoint

pip install matplotlib

For those new to python, [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/) is recommended. Conda makes it easier to install packages and manage environments than the basic python installation.

catpsolver/v1) consists of five files, including tap_sim.py, func_sim.py, reac_odes.py, vari_form.py and input_file.csv. How these files are related to one another are shown in Figure 1.



When the simulator is run, an output directory will be generated with all the relevent information stored in appropriate directories. 

