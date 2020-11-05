
* [Overview](https://github.com/medford-group/TAPsolver/tree/master)
* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Examples](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion)

## File structure

The current version of [TAPsolver](https://github.com/medford-group/TAPsolver/tree/master/tapsolver/v1) consists of five files, including tap_sim.py, func_sim.py, reac_odes.py, vari_form.py and input_file.csv. 

These separate scripts are actively being bundled into a single package that can be installed directly in anaconda.

When the simulator is run, an output directory will be generated with all the relevent information stored in appropriate directories, depending on what information is desired.

## User Input File Structure

A [csv file](https://github.com/medford-group/TAPsolver/blob/master/tapsolver/v1/input_file.csv) is used as the input to the simulator. Four main componenents of TAPSolver input are described and include "Reactor Information", "Feed and Surface Composition", "Data Storage Options", and "Reaction Information". The terms in each of these components are discussed below.

### Reactor Information

#### Zone Length

Options: 0 to inf (float)

Description: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

#### Zone Void

Options: 0 to inf (float)

Description: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

#### Reactor Radius

Options: 0 to inf (float)

Description: The radius of the tap reactor. It is common to keep the radius [significantly smaller](https://www.math.wustl.edu/~feres/ReactionDiffusionWallFerGreg.pdf) than the length of the reactor to maintain nearly unidirectional diffusion.

#### Reactor Temperature

Options: 0 to inf (float)

Description: Establish a constant TAP reactor temperature, which is acceptable due to the [underlying theory of TAP](https://pubs.rsc.org/en/content/articlehtml/2017/cy/c7cy00678k). Should use K as unit.

#### Mesh Size

Options: 0 to inf (int)

Description: How refined the simulation will be spatially (the length of the elements). Similar to time steps, a larger mesh size is typically required for challenging (stiff) simulations.

#### Output Folder Name

Options: example or path/to/example (excluding ./ at start of path)

Description: This identifies the name of the directory and its location relative to the core TAPsolver scripts. It is important to come up with directory names with each simulation to avoid deleting previous results.

Examples:

save in current directory: example

save in nested directory: stored_data/example

save in previous directory: ../example

If no experimental data is being used, then write None.

#### Experimental Data Folder

Options: ./example_data or ./path/to/data (including ./ at start of path)

Description: When making comparisons directly to experimental data or when fitting parameters, specify the location of these files with this parameter.

#### Reference Diffusion Inert

Options: 0 to inf (float)

Description: The known diffusion coefficient in the material of some gas, which can be scaled with Graham's Law (a common method in the TAP community) to find the diffusion coefficient of other gasses without having to measure them experimentally..

#### Reference Diffusion Catalyst

Options:

Description: See 'Reference Diffusion Inert'. At times, the catalyst zone is narrow enough to assume that the diffusion will not vary between the inert and catalyst region (similar to 'Void Fraction Catalyst'). 

#### Reference Temperature

Options: 0 to inf (float)

Description: The temperature at which the reference diffusion coefficient was measured. Allows the diffusion coefficient to be [scaled by temperature](https://en.wikipedia.org/wiki/Knudsen_diffusion), too

#### Reference Mass

Options: 0 to inf (float, though commonly an int for amu) 

Description: Mass of the gas with known diffusion coefficients (see 'Reference Diffusion' Sections)

### Feed and Surface Composition

#### Intensity

Options: (0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),...

Description: Allows for pump-probe style simulations with a similar input to the 'Pulse Ratio'

Example: 0.0,0.03,0.0,0.0

#### Time

Options: (0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),...

Description: Allows for pump-probe style simulations with a similar input to the 'Pulse Ratio'

Example: 0.0,0.03,0.0,0.0

#### Mass

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Description: Mass of each monitored gas species (in amu).

Example: 28,16,44,40

#### Initial Concentration

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Desciption: The initial concentration of the surface species at each mesh element in the thin catalyst zone.

### Reaction Information

Description: Used to specify the elementary reactions involved in the microkinetic model. 

