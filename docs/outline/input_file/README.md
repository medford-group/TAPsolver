## File structure

The current version of [TAPsolver](https://github.com/medford-group/TAPsolver/tree/master/tapsolver/v1) consists of five files, including tap_sim.py, func_sim.py, reac_odes.py, vari_form.py and input_file.csv. How these files are related to one another are shown in Figure 1.



When the simulator is run, an output directory will be generated with all the relevent information stored in appropriate directories.

Two forms of TAPsolver currently exist. One for FEniCS version 2017 and the other for FEniCS version 2018.

## User Input File Structure

A [csv file](https://github.com/medford-group/TAPsolver/blob/master/tapsolver/v1/input_file.csv) is used as the input to the simulator. Four main componenents of TAPSolver input are described and include "Reactor Information", "Feed and Surface Composition", "Data Storage Options", and "Reaction Information". The terms in each of these components are discussed below.

## Reactor Information

### Pulse Duration

Options: 0 to inf (float)

Desciption: The amount of time (in seconds) each pulse is simulated over (typically how long it takes for outlet flux curves to reach zero).

### Time Steps

Options: 0 to inf (int)

Desciption: The number of time steps you want to occur during each pulse simulation. This will determine how smooth (and accurate) your simulation will be. The more complex the simulation, the more time steps will be required (due to the stiffness of the solution). A value of 500 is recommended as a minimum for relatively simple simulations.

### Reactor Length

Options: 0 to inf (float)

Desciption: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

### Mesh Size

Options: 0 to inf (int)

Desciption: How refined the simulation will be spatially (the length of the elements). Similar to time steps, a larger mesh size is typically required for challenging (stiff) simulations.

### Catalyst Fraction

Options: 0 to 1 (float)

Desciption: This is the fraction of the reactor occupied by the catalyst zone. For TAP reactors, the fraction will typically stay below 0.05 ([to maintain the thin zone assumption](https://www.sciencedirect.com/science/article/pii/S000925099800534X)).

### Reactor Radius

Options: 0 to inf (float)

Desciption: The radius of the tap reactor. It is common to keep the radius [significantly smaller](https://www.math.wustl.edu/~feres/ReactionDiffusionWallFerGreg.pdf) than the length of the reactor to maintain nearly unidirectional diffusion.

### Void Fraction Inert

Options: 0 to 1 (float)

Desciption: How much open space is available in the [inert zone](https://en.wikipedia.org/wiki/Porosity) 

### Void Fraction Catalyst

Options: 0 to inf (float)

Desciption: See 'Void Fraction Inert'. At times, the catalyst zone is narrow
enough to assume that the void fraction will not vary between the inert and catalyst region.

### Reactor Temperature

Options: 0 to inf (float)

Desciption: Establish a constant TAP reactor temperature, which is acceptable due to the [underlying theory of TAP](https://pubs.rsc.org/en/content/articlehtml/2017/cy/c7cy00678k). Should use K as unit.

### Reference Diffusion Inert

Options: 0 to inf (float)

Desciption: The known diffusion coefficient in the material of some gas, which can be scaled with Graham's Law (a common method in the TAP community) to find the diffusion coefficient of other gasses without having to measure them experimentally..

### Reference Diffusion Catalyst

Options:

Desciption: See 'Reference Diffusion Inert'. At times, the catalyst zone is narrow enough to assume that the diffusion will not vary between the inert and catalyst region (similar to 'Void Fraction Catalyst'). 

### Reference Temperature

Options: 0 to inf (float)

Desciption: The temperature at which the reference diffusion coefficient was measured. Allows the diffusion coefficient to be [scaled by temperature](https://en.wikipedia.org/wiki/Knudsen_diffusion), too

### Reference Mass

Options: 0 to inf (float, though commonly an int for amu) 

Desciption: Mass of the gas with known diffusion coefficients (see 'Reference Diffusion' Sections)

### Output Folder Name

Options: example or path/to/example (excluding ./ at start of path)

Desciption: This identifies the name of the directory and its location relative to the core TAPsolver scripts. It is important to come up with directory names with each simulation to avoid deleting previous results.

Examples: 

save in current directory: example

save in nested directory: stored_data/example

save in previous directory: ../example

### Experimental Data Folder

Options: ./example_data or ./path/to/data (including ./ at start of path)

Desciption: When making comparisons directly to experimental data or when fitting parameters, specify the location of these files with this parameter.

### Noise

Options: 'TRUE' or 'FALSE'

Desciption: Include simple tap noise in the [Flux Data](https://biblio.ugent.be/publication/8525766/file/8525861.pdf)

### Reactor Type

Options: 'tap' or 'tap-diffusion'

Desciption: Specify what type of reactor is being considered (including the transport). Currently need further development before 'tap-diffusion' can be properly applied.

### Theta

Options: '0', '0.5' or '1'

Desciption: Specify the step type fenics uses to solve 

### Solver Method

Options: - - - -

Desciption: Will include options for implicit and explicit time stepping.

## Feed and Surface Composition

#### Number of Reactants

Options: 1 to inf (int)

Desciption: How many gas species are being pulsed in and directly interact with the material(excluding any inerts).

#### Number of Pulses 

Options: 1 to inf (int)

Desciption: The total number of pulses included in the simulation, allowing for state defining or state altering experiments.

#### Reference Pulse Size

Options: 0 to inf (float)

Desciption: A reference pulse intensity (nanomoles) for all other species being fed to be based off (see 'Pulse Ratio').

#### Pulse Ratio

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Desciption: Relative to the 'Reference Pulse Size', what is the intensity of each gas species being fed to the reactor.

Example: 1,0.5,1.5,0,1

#### Pulse Time

Options: (0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),...

Desciption: Allows for pump-probe style simulations with a similar input to the 'Pulse Ratio'

Example: 0.0,0.03,0.0,0.0

#### Mass List

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Desciption: Mass of each monitored gas species (in amu).

Example: 28,16,44,40

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

#### Optimization Method

Options: 'Newton-CG', 'BFGS', 'SLSQP', 'CG', 'basinhopping', 'COBYLA', 'TNC'

Desciption: See the [scipy documentation](https://docs.scipy.org/doc/scipy/reference/optimize.html) for more details on each of these methods. Recommendations will be included in future versions of the documentation.

#### Objective Points

Options: - - - - (subject to change)

Desciption: The number and distribution of the points that will be used to fit the curves in the 'optimization method'. 

#### RRM Analysis 

Options: 'TRUE' or 'FALSE'

Desciption: (Subject to change)

#### MKM Analysis 

Options: 'TRUE' or 'FALSE'

Desciption: (Subject to change) Primarily used to store thin zone concentrations and rates.

#### Petal Plots

Options: 'TRUE' or 'FALSE'

Desciption: (Subject to change) Generate simple petal plots showing the rates relative to different concentrations.

### Reaction Information

Description: Used to specify the elementary reactions involved in the microkinetic model. 

