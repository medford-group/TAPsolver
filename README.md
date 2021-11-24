
<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/tapsolver_logo.jpg">
</p>

TAPsolver is an open-source, python program for the analysis of Temporal Analysis of Products (TAP) reactor experiments. The methods of analysis include:

* *PDE-based simulations, including state-altering, state-defining and pump-probe operating conditions*
* *PDE-constrained optimization of kinetic parameters*
* *Hessian assesment of parameter uncertainties*
* *Ensemble-based evaluation of initial state/reactor configuration uncertainties*
* *Model-based design of experiments and model discrimination (under development)*
* *Efficient kinetic parameter uncertainty propagation to TAP experiments (under development)*
* *Bayesian parameter estimation (under development)*

<p align="center">
  <img src="https://github.com/medford-group/TAPsolver/blob/master/docs/figures/tapsolver_concept.gif">
</p>

The python packages [FEniCS](https://fenicsproject.org/)  (Finite Element Computational Software) and [Dolfin-Adjoint](http://www.dolfin-adjoint.org/en/latest/) act as the foundation of TAPsolver and allow for the flexible and efficient application of [algorithmic differentiation](https://towardsdatascience.com/automatic-differentiation-explained-b4ba8e60c2ad). TAPsolver was also closely developed alongside [TAPSAP](https://github.com/IdahoLabResearch/tapsap), which is used for processing and statistically analyzing raw TAP experimental data. 


A primary goal of TAPsolver is to make new methods of TAP analysis easily accesable upon publication, leading to community driven development and the unification of workflows. If you have a new approach for processing TAP data and would like to implement it in TAPsolver, please read the [Questions & Development section](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion).

* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Documentation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/Documentation)
* [Examples](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion)

### Collaborators and Funding

TAPsolver was developed through a colloboration between the Georgia Institute of Technology and Idaho National Lab. The publication associated with TAPsolver can be found [here](https://arxiv.org/abs/2008.13584). 

Financial support for the development of TAPsolver was provided by the U.S. Department of Energy (USDOE), Office of Energy Efficiency and Renewable Energy (EERE), Advanced Manufacturing Office Next Generation R\&D Projects under contract no. DE-AC07-05ID14517.

### Other Software
Idaho National Laboratory is a cutting edge research facility which is a constantly producing high quality research and software. Feel free to take a look at our other software and scientific offerings at:

[Primary Technology Offerings Page](https://www.inl.gov/inl-initiatives/technology-deployment)

[Supported Open Source Software](https://github.com/idaholab)

[Raw Experiment Open Source Software](https://github.com/IdahoLabResearch)

[Unsupported Open Source Software](https://github.com/IdahoLabCuttingBoard)

### License

Copyright 2020 Battelle Energy Alliance, LLC

Licensed under the GPL-2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  https://opensource.org/licenses/GPL-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


Licensing
-----
This software is licensed under the terms you may find in the file named "LICENSE" in this directory.


Developers
-----
By contributing to this software project, you are agreeing to the following terms and conditions for your contributions:

You agree your contributions are submitted under the GPL-2 license. You represent you are authorized to make the contributions and grant the license. If your employer has rights to intellectual property that includes your contributions, you represent that you have received permission to make contributions and grant the required license on behalf of that employer.
