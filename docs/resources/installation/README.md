
NOTE TO ACS WORKSHOP ATTENDEES: Please download the examples folder above which will be used throughout the workshop. If you have any issues installing the program, please email me at ayonge3@gatech.edu

# Mac and Ubuntu Installation

Using [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install) is the simplest way to install [FEniCS](https://fenicsproject.org/) (on [Ubuntu](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0) systems, especially) and is highly recommended. To install the necessary programs for running tapsolver, enter the following lines in the terminal (make sure conda-forge is in your list of channels by running the following commaned: "conda config --append channels conda-forge"):

	
	conda create -n tapsolver python=3.7 fenics jsonpickle -c conda-forge/label/cf202003
	
	pip install --upgrade git+https://github.com/medford-group/TAPsolver.git@mbdoe_publication
	
	pip install --upgrade git+https://github.com/dolfin-adjoint/pyadjoint.git@faster-ufl
	
	pip install notebook

# Windows Installation

From the FEniCS installation page: "To install FEniCS on Windows 10, enable the Windows Subsystem for Linux and install the Ubuntu distribution."

Install the [ubuntu operating system](https://ubuntu.com/wsl) on windows and then follow the instructions presented in the 'Mac and Ubuntu Installation' section above.
