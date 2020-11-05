
* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Documentation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/Documentation)
* [Examples](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion) 

# TAPsolver
TAPsolver is an open-source, python based Temporal Analysis of Products (TAP) simulation and analysis tool. The python packages [FEniCS](https://fenicsproject.org/)  (Finite Element Computational Software) and [Dolfin-Adjoint](http://www.dolfin-adjoint.org/en/latest/) act as the foundation of TAPsolver and allow for the flexible and efficient application of [automatic differentiation](https://towardsdatascience.com/automatic-differentiation-explained-b4ba8e60c2ad). 

A primary goal of TAPsolver is to make new methods of TAP analysis easily accesable upon publication, leading to community driven development and the unification of workflows. If you have a new approach for processing TAP data and would like to implement it in TAPsolver, please read the [Questions & Development section](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion).

**Simulation options include**:

* *State-defining*
* *State-altering*
* *Multi-site* 
* *Pump-probe* 

**Processing options include**:

* *Sensitivity analysis*
* *Parameter optimization*
* *Uncertainty quantification*
* *Moment analysis*
* *Y-procedure / G-procedure*

<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/cleanOptimizationExample.gif">
</p>

TAPsolver was developed through a colloboration between the Georgia Institute of Technology and Idaho National Lab. The publication associated with TAPsolver can be found [here](https://arxiv.org/abs/2008.13584). 

Financial support for the development of TAPsolver was provided by the U.S. Department of Energy (USDOE), Office of Energy Efficiency and Renewable Energy (EERE), Advanced Manufacturing Office Next Generation R\&D Projects under contract no. DE-AC07-05ID14517.
