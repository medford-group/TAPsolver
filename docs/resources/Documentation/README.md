
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

## Reactor Information

### Pulse Duration

Options: 0 to inf (float)

Description: The amount of time (in seconds) each pulse is simulated over (typically how long it takes for outlet flux curves to reach zero).

### Time Steps

Options: 0 to inf (int)

Description: The number of time steps you want to occur during each pulse simulation. This will determine how smooth (and accurate) your simulation will be. The more complex the simulation, the more time steps will be required (due to the stiffness of the solution). A value of 500 is recommended as a minimum for relatively simple simulations.

### Reactor Length

Options: 0 to inf (float)

Description: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

### Mesh Size

Options: 0 to inf (int)

Description: How refined the simulation will be spatially (the length of the elements). Similar to time steps, a larger mesh size is typically required for challenging (stiff) simulations.

### Catalyst Fraction

Options: 0 to 1 (float)

Description: This is the fraction of the reactor occupied by the catalyst zone. For TAP reactors, the fraction will typically stay below 0.05 ([to maintain the thin zone assumption](https://www.sciencedirect.com/science/article/pii/S000925099800534X)).

### Reactor Radius

Options: 0 to inf (float)

Description: The radius of the tap reactor. It is common to keep the radius [significantly smaller](https://www.math.wustl.edu/~feres/ReactionDiffusionWallFerGreg.pdf) than the length of the reactor to maintain nearly unidirectional diffusion.

### Catalyst Location

Options: 0 to 1 (float)

Description: Where to locate the center of the catalyst zone in the TAP reactor (typically at 0.5, the halfway point, but is also often placed somewhere between 0.4 and 0.6).

### Void Fraction Inert

Options: 0 to 1 (float)

Description: How much open space is available in the [inert zone](https://en.wikipedia.org/wiki/Porosity) 

### Void Fraction Catalyst

Options: 0 to inf (float)

Description: See 'Void Fraction Inert'. At times, the catalyst zone is narrow
enough to assume that the void fraction will not vary between the inert and catalyst region.

### Reactor Temperature

Options: 0 to inf (float)

Description: Establish a constant TAP reactor temperature, which is acceptable due to the [underlying theory of TAP](https://pubs.rsc.org/en/content/articlehtml/2017/cy/c7cy00678k). Should use K as unit.

### Reference Diffusion Inert

Options: 0 to inf (float)

Description: The known diffusion coefficient in the material of some gas, which can be scaled with Graham's Law (a common method in the TAP community) to find the diffusion coefficient of other gasses without having to measure them experimentally..

### Reference Diffusion Catalyst

Options:

Description: See 'Reference Diffusion Inert'. At times, the catalyst zone is narrow enough to assume that the diffusion will not vary between the inert and catalyst region (similar to 'Void Fraction Catalyst'). 

### Reference Temperature

Options: 0 to inf (float)

Description: The temperature at which the reference diffusion coefficient was measured. Allows the diffusion coefficient to be [scaled by temperature](https://en.wikipedia.org/wiki/Knudsen_diffusion), too

### Reference Mass

Options: 0 to inf (float, though commonly an int for amu) 

Description: Mass of the gas with known diffusion coefficients (see 'Reference Diffusion' Sections)

### Output Folder Name

Options: example or path/to/example (excluding ./ at start of path)

Description: This identifies the name of the directory and its location relative to the core TAPsolver scripts. It is important to come up with directory names with each simulation to avoid deleting previous results.

Examples: 

save in current directory: example

save in nested directory: stored_data/example

save in previous directory: ../example

If no experimental data is being used, then write None.

### Experimental Data Folder

Options: ./example_data or ./path/to/data (including ./ at start of path)

Description: When making comparisons directly to experimental data or when fitting parameters, specify the location of these files with this parameter.

### Noise

Options: 'TRUE' or 'FALSE'

Description: Include simple tap noise in the [Flux Data](https://biblio.ugent.be/publication/8525766/file/8525861.pdf)

### Reactor Type

Options: 'tap' or 'tap-diffusion'

Description: Specify what type of reactor is being considered (including the transport). Currently need further development before 'tap-diffusion' can be properly applied.

### Knudsen Test

Options: 'TRUE' or 'FALSE'

Description: Used to test if the inert data is within the knudsen regime. Should return approximately 0.31 for all of the inerts.

## Feed and Surface Composition

#### Number of Reactants

Options: 1 to inf (int)

Description: How many gas species are being pulsed in and directly interact with the material(excluding any inerts).

#### Number of Pulses 

Options: 1 to inf (int)

Description: The total number of pulses included in the simulation, allowing for state defining or state altering experiments.

#### Reference Pulse Size

Options: 0 to inf (float)

Description: A reference pulse intensity (nanomoles) for all other species being fed to be based off (see 'Pulse Ratio').

#### Pulse Ratio

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Description: Relative to the 'Reference Pulse Size', what is the intensity of each gas species being fed to the reactor.

Example: 1,0.5,1.5,0,1

#### Pulse Time

Options: (0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),...

Description: Allows for pump-probe style simulations with a similar input to the 'Pulse Ratio'

Example: 0.0,0.03,0.0,0.0

#### Mass List

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Description: Mass of each monitored gas species (in amu).

Example: 28,16,44,40

#### Advection Value

Options: 0 to inf (float)

Description: The value of the advection term (some velocity that is helping carry to species).

#### Initial Surface Composition

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Desciption: The initial concentration of the surface species at each mesh element in the thin catalyst zone.

#### Number of Active Sites

Options: 1 to 4 (int)

Desciption: The total number of unique active sites included in the simulation.

#### Number of Inerts

Options: 1 to inf (int)

Desciption: The total number of inerts being pulsed in the reactor.

### Data Storage Options

#### Store Outlet Flux

Options: 'TRUE' or 'FALSE'

Desciption: Tells the script whether to save the outlet flux data or not

#### Store Graph

Options: 'TRUE' or 'FALSE'

Desciption: Tells the script whether to save the outlet flux graph or not.

#### Display Experimental Data 

Options: 'TRUE' or 'FALSE'

Desciption: Tells the script whether or not to show the experimental data on the same graph as the outlet flux.

#### Display Graph

Options: 'TRUE' or 'FALSE'

Desciption: Tells the script whether or not to display the outlet flux graph after the simulation

#### Sensitivity Analysis

Options: 'TRUE' or 'FALSE'

Desciption: Tells the script whether or not a sensitivity analysis on the reaction mechanism should be performed.

#### Fit Parameters

Options: 'TRUE' or 'FALSE'

Desciption: Tells the script whether or not parameter fitting should be performed on the provided experimental data.

#### RRM Analysis 

Options: 'TRUE' or 'FALSE'

Desciption: (Subject to change)

#### Thin-Zone Analysis 

Options: 'TRUE' or 'FALSE'

Desciption: Store thin zone concentrations and rates.

#### Display Objective Points

Options: 'TRUE' or 'FALSE'

Description: Display the points used for optimization on the output graph.

#### Fitting Gif

Options: 'TRUE' or 'FALSE'

Description: Generate a gif showing the outlet curves for each step in the optimization routine.

#### Uncertainty Quantification

Options: 'TRUE' or 'FALSE'

Description: Used to calculate the uncertainty in parameters found during optimization.

### Reaction Information

Description: Used to specify the elementary reactions involved in the microkinetic model. 


