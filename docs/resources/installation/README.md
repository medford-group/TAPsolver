
* [Overview](https://github.com/medford-group/TAPsolver/tree/master)
* [Interface options](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/interfaceOptions)
* [Examples](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples/coAdsorption)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion)
* [Documentation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file)

# Mac and Ubuntu Installation

[Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install) is the simplest way to install [FEniCS](https://fenicsproject.org/) (on [Ubuntu](https://tutorials.ubuntu.com/tutorial/tutorial-ubuntu-on-windows#0) systems, especially) and is highly recommended. It is possible to install FEniCS in other ways, but Conda was found to be the quickest.

To install the necessary programs for running tapsolver, enter the following lines in the terminal:

	conda create -n fenics-2019 -c conda-forge fenics

	conda install -c conda-forge dolfin-adjoint

	conda install -c conda-forge imageio

	conda install -c anaconda pandas

	conda install xlrd

If you are having trouble installing FEniCS, contact Adam Yonge (ayonge3@gatech.edu) or the [FEniCS developers](https://fenicsproject.org/community/).

# Windows Installation

From the FEniCS installation page: "To install FEniCS on Windows 10, enable the Windows Subsystem for Linux and install the Ubuntu distribution."

To install the program, install the ubuntu operating system on windows and then follow the instructions presented in the 'Mac and Ubuntu Installation' section above.

# Installation with hippylib and muq2

        conda create -n fenics-2019.1 -c uvilla -c conda-forge fenics==2019.1.0 matplotlib scipy jupyter

        pip install -e git+https://github.com/hippylib/hippylib@master#egg=hippylib

        conda install -c conda-forge dolfin-adjoint

        conda install -c conda-forge imageio

        conda install -c anaconda pandas

        (make sure conda-forge is in your list of channels by running the following commaned: "conda config --append channels conda-forge")

        conda install -c rmg muq2

        pip install git+https://bitbucket.org/dolfin-adjoint/pyadjoint.git@master

        conda install xlrd

