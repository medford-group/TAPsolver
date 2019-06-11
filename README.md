# TAPsolver

## Installation

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install) is the simplest way to install [FEniCS](https://fenicsproject.org/) (on [Ubuntu](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0) systems, especially) and is highly recommended. It is possible to install FEniCS in other ways, but Conda was found to be the quickest.

To install FEniCS through conda, run the following the terminal:

	conda create -n fenicsproject -c conda-forge fenics

	source activate fenicsproject

	conda install -c conda-forge dolfin-adjoint

	pip install matplotlib

If you are having trouble installing FEniCS, contact Adam Yonge (ayonge3@gatech.edu) or the [FEniCS developers](https://fenicsproject.org/community/).

## Running TAPsolver

With the five core files downloaded in the working directory and with the fenicsproject environment activated, run the following command:

	python tap_sim.py

This will call the simulation process and store all desired information in the appropriate folders.  

This is the current format of TAPsolver and is subject to change. Steps will be taken to merge this process with other methods of TAP analysis developed by collaborators in the R programming language. 

## Resources

An EXPLANATION OF THE INPUT FILE  and how to properly manipulate it can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/outline/input_file).

A simple TUTORIAL can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples/coAdsorption)

## Example Output

![Example curves](./docs/figures/CO.gif)
![Example curves](./docs/figures/flux_data.png)
