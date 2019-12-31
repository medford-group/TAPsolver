
# TAPsolver

* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Interface options](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/ineterfaceOptions)
* [Examples](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples/coAdsorption)
* [Documentation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion) 

## Overview
TAPsolver is an open-source, python based Temporal Analysis of Products (TAP) simulation and analysis tool. The python packages [FEniCS](https://fenicsproject.org/)  (Finite Element Computational Software) and [Dolfin-Adjoint](http://www.dolfin-adjoint.org/en/latest/) act as the foundation of TAPsolver and allow for the flexible and efficient application of [automatic differentiation](https://towardsdatascience.com/automatic-differentiation-explained-b4ba8e60c2ad). 

Simulation options include:

* State-defining
* State-altering
* Multi-site 
* Pump-probe 

Processing options include:

* Sensitivity analysis
* Parameter optimization
* Moment analysis
* Y-procedure / G-procedure

By generalizing TAP simulation and analysis options and providing roboust documentation, TAPsolver should offer the TAP community a centralized location for method development and application, as well as a consistent tool for processing experimental data. 


<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/expSimExample.png">
</p>

<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/CO.gif">
</p>

# Running TAPsolver

An explanation of the input file and how to properly manipulate it can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file) and a simple tutorial can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples/coAdsorption).

Currently, TAPsolver consists of five files, including the input file. With the five core files downloaded in the working directory and with the fenicsproject environment activated, run the following command:

	python tap_sim.py

This will call the simulation process and store all desired information in the appropriate folders.  

This is the current format of TAPsolver and is subject to change. Steps will be taken to merge this process with other methods of TAP analysis developed by collaborators in the R programming language. 
