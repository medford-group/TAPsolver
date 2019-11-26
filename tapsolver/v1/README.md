
# TAPsolver

The primary features of TAPsolver are:
 * Forward simulations of common TAP experiments, including:
   * State-defining
   * State-altering
   * Multi-site 
   * Pump-probe 
 * Sensitivity Analysis of thin-zone and outlet data
 * Kinetic parameter fitting for reaction mechanisms

<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/expSimExample.png">
</p>

<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/CO.gif">
</p>

# Installation

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install) is the simplest way to install [FEniCS](https://fenicsproject.org/) (on [Ubuntu](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0) systems, especially) and is highly recommended. It is possible to install FEniCS in other ways, but Conda was found to be the quickest.

To properly establish the environment for FEniCS through conda, run the following commands:


	conda create -n fenics-2019.1 -c uvilla -c conda-forge fenics==2019.1.0 matplotlib scipy jupyter

	pip install -e git+https://github.com/hippylib/hippylib@master#egg=hippylib

	conda install -c conda-forge dolfin-adjoint

	conda install -c conda-forge imageio

	conda install -c anaconda pandas

	(make sure conda-forge is in your list of channels by running the following commaned: "conda config --append channels conda-forge")

	conda install -c rmg muq2 

	pip install git+https://bitbucket.org/dolfin-adjoint/pyadjoint.git@master

	conda install xlrd

If you are having trouble installing FEniCS, contact Adam Yonge (ayonge3@gatech.edu) or the [FEniCS developers](https://fenicsproject.org/community/).

# Running TAPsolver

An explanation of the input file and how to properly manipulate it can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file) and a simple tutorial can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples/coAdsorption).

Currently, TAPsolver consists of five files, including the input file. With the five core files downloaded in the working directory and with the fenicsproject environment activated, run the following command:

	python tap_sim.py

This will call the simulation process and store all desired information in the appropriate folders.  

This is the current format of TAPsolver and is subject to change. Steps will be taken to merge this process with other methods of TAP analysis developed by collaborators in the R programming language. 
