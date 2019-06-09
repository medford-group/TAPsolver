## File structure

The current version of [TAPsolver](https://github.com/medford-group/TAPsolver/tree/master/tapsolver/v1) consists of five files, including tap_sim.py, func_sim.py, reac_odes.py, vari_form.py and input_file.csv. How these files are related to one another are shown in Figure 1.



When the simulator is run, an output directory will be generated with all the relevent information stored in appropriate directories.

Two forms of TAPsolver currently exist. One for FEniCS version 2017 and the other for FEniCS version 2018.

## User Input File Structure

A csv file is used as the input to the simulator. Four main componenents of TAPSolver input are described and include "Reactor Information", "Feed and Surface Composition", "Data Storage Options", and "Reaction Information". The terms in each of these components are discussed below.

### Reactor Information

#### Pulse Duration

Options: 0 to inf (float)

Desciption: The amount of time (in seconds) each pulse is simulated over (typically how long it takes for outlet flux curves to reach zero).

#### Time Steps

Options: 0 to inf (int)

Desciption: The number of time steps you want to occur during each pulse simulation. This will determine how smooth (and accurate) your simulation will be. The more complex the simulation, the more time steps will be required (due to the stiffness of the solution). A value of 500 is recommended as a minimum for relatively simple simulations.

#### Reactor Length

Options: 0 to inf (float)

Desciption: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

#### Mesh Size

Options: 0 to inf (int)

Desciption: How refined the simulation will be spatially (the length of the elements). Similar to time steps, a larger mesh size is typically required for challenging (stiff) simulations.

#### Catalyst Fraction

Options: 0 to 1 (float)

Desciption: This is the fraction of the reactor occupied by the catalyst zone. For TAP reactors, the fraction will typically stay below 0.05 ([to maintain the thin zone assumption](https://www.sciencedirect.com/science/article/pii/S000925099800534X)).

#### Reactor Radius

Options: 0 to inf (float)

Desciption: The radius of the tap reactor. It is common to keep the radius [significantly smaller](https://www.math.wustl.edu/~feres/ReactionDiffusionWallFerGreg.pdf) than the length of the reactor to maintain nearly unidirectional diffusion.

#### Void Fraction Inert

Options: 0 to 1 (float)

Desciption: How much open space is available in the [inert zone](https://en.wikipedia.org/wiki/Porosity) 

#### Void Fraction Catalyst

Options: 0 to inf (float)

Desciption: See 'Void Fraction Inert'. At times, the catalyst zone is narrow
enough to assume that the void fraction will not vary between the inert and catalyst region.

#### Reactor Temperature

Options: 0 to inf (float)

Desciption: Establish a constant TAP reactor temperature, which is acceptable due to the [underlying theory of TAP](https://pubs.rsc.org/en/content/articlehtml/2017/cy/c7cy00678k). Should use K as unit.

#### Reference Diffusion Inert

Options: 0 to inf (float)

Desciption: The known diffusion coefficient in the material of some gas, which can be scaled with Graham's Law (a common method in the TAP community) to find the diffusion coefficient of other gasses without having to measure them experimentally..

#### Reference Diffusion Catalyst

Options:

Desciption: See 'Reference Diffusion Inert'. At times, the catalyst zone is narrow enough to assume that the diffusion will not vary between the inert and catalyst region (similar to 'Void Fraction Catalyst'). 

#### Reference Temperature

Options: 0 to inf (float)

Desciption: The temperature at which the reference diffusion coefficient was measured. Allows the diffusion coefficient to be [scaled by temperature](https://en.wikipedia.org/wiki/Knudsen_diffusion), too

#### Reference Mass

Options: 0 to inf (float, though commonly an int for amu) 

Desciption: 

#### Output Folder Name

Options:

Desciption:

#### Experimental Data Folder

Options:

Desciption:

#### Noise

Options:

Desciption:

#### Reactor Type

Options:

Desciption:

#### Theta

Options:

Desciption:

#### Solver Method

### Feed and Surface Composition

Options:

Desciption:

#### Number of Reactants

Options:

Desciption:

#### Number of Pulses 

Options:

Desciption:

#### Reference Pulse Size

Options:

Desciption:

#### Pulse Ratio

Options:

Desciption:

#### Pulse Time

Options:

Desciption:

#### Mass List

Options:

Desciption:

#### Initial Surface Composition

Options:

Desciption:

#### Number of Active Sites

Options:

Desciption:

#### Number of Inerts

Options:

Desciption:

### Data Storage Options

#### Store Outlet Flux

Options:

Desciption:

#### Store Graph

Options:

Desciption:

#### Display Experimental Data 

Options:

Desciption:

#### Display Graph

Options:

Desciption:

#### Sensitivity Analysis

Options:

Desciption:

#### Fit Parameters

Options:

Desciption:

#### Optimization Method

Options:

Desciption:

#### Objective Points

Options:

Desciption:

#### RRM Analysis 

Options:

Desciption:

#### MKM Analysis 

Options:

Desciption:

#### Petal Plots

Options:

Desciption:

### Reaction Information



















