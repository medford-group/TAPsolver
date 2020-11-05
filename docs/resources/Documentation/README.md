
* [Overview](https://github.com/medford-group/TAPsolver/tree/master)
* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Examples](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/examples)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion)

#User Input File Structure

A [csv file](https://github.com/medford-group/TAPsolver/blob/master/tapsolver/v1/input_file.csv) is used as the input to the simulator. Four main componenents of TAPSolver input are described and include "Reactor Information", "Feed and Surface Composition", "Data Storage Options", and "Reaction Information". The terms in each of these components are discussed below.

## Reactor Information

### Zone Length

Options: 0 to inf (float)

Description: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

### Zone Void

Options: 0 to inf (float)

Description: How long the reactor is (typically cm), including both inert zones and the thin catalyst zone.

### Reactor Radius

Options: 0 to inf (float)

Description: The radius of the tap reactor. It is common to keep the radius [significantly smaller](https://www.math.wustl.edu/~feres/ReactionDiffusionWallFerGreg.pdf) than the length of the reactor to maintain nearly unidirectional diffusion.

### Reactor Temperature

Options: 0 to inf (float)

Description: Establish a constant TAP reactor temperature, which is acceptable due to the [underlying theory of TAP](https://pubs.rsc.org/en/content/articlehtml/2017/cy/c7cy00678k). Should use K as unit.

### Mesh Size

Options: 0 to inf (int)

Description: How refined the simulation will be spatially (the length of the elements). Similar to time steps, a larger mesh size is typically required for challenging (stiff) simulations.

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

## Feed and Surface Composition

### Intensity

Options: (0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),...

Description: Allows for pump-probe style simulations with a similar input to the 'Pulse Ratio'

Example: 0.0,0.03,0.0,0.0

### Time

Options: (0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),(0 to 'Pulse Duration' (float)),...

Description: Allows for pump-probe style simulations with a similar input to the 'Pulse Ratio'

Example: 0.0,0.03,0.0,0.0

### Mass

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Description: Mass of each monitored gas species (in amu).

Example: 28,16,44,40

### Initial Concentration

Options: (0 to inf (float)),(0 to inf (float)),(0 to inf (float)),...

Desciption: The initial concentration of the surface species at each mesh element in the thin catalyst zone.

## Reaction Information

Description: Used to specify the elementary reactions involved in the microkinetic model. 

# Processing Functions

## Forward Problem (Simulation)

### run_tapsolver(sim_time,add_noise=False,pulse_num = 1,catalyst_data=False,sim_file = './sim_file.csv')

sim_time: The time (in seconds) you would like the simulation to run for. Values can be any positive float or interger.

add_noise: Should white noise be included in the simulation? (True or False)

pulse_num: How many pulses should be simulated? (default = 1, any int > 1)

catalyst_data: Would you like information related to the catalyst zone (surface/gas concentrations, rates) to be stored? (True or False)

sim_file: What previously defined input file would you like to use to run your simulation (see notes above).

## Inverse Problem (Optimization)

### fit_tap(sim_time,optim = 'L-BFGS-B',sim_file = './sim_file.csv')

sim_time: The time (in seconds) you would like the simulation to run for. Values can be any positive float or interger.

optim: What optimization routine would you like to use? Default is 'L-BFGS-B'

sim_file: What previously defined input file would you like to use to run your simulation (see notes above).

### run_sensitivity(sim_time,sens_type=None,sim_file = './sim_file.csv')

sim_time: The time (in seconds) you would like the simulation to run for. Values can be any positive float or interger.

sens_type: Run a total sensitivity analysis (used during optimization) or a transient sensitivity analysis (the influence each kinetic parameter has on the outlet flux over time). Options are 'total' or 'trans'.

sim_file: What previously defined input file would you like to use to run your analysis (see notes above).

### run_uncertainty(sim_time,sim_file = './sim_file.csv')

sim_time: The time (in seconds) you would like the simulation to run for. Values can be any positive float or interger.

sim_file: What previously defined input file would you like to use to run your analysis (see notes above).

## Visualization 

### flux_graph(sim_file = './sim_file.csv',pulse=None,disp_exper=False,disp_analytic=False,disp_objective=False,show_graph=True,store_graph=False,output_name='./flux.png')

sim_file: What previously defined input file would you like to use to run your analysis (see notes above).

pulse: What pulses would you like to visualized? Can be an int or list of ints representing each pulse.

disp_exper: Display the referenced experimental data? (Default False)

disp_analytic: Display the analytical solution (pure diffsion of a gas species)? Default False.

disp_objective: Display the objective points from the experimental data used to fit parameters? (Default False)

show_graph: Display the graph after it is generated? (Default True)

store_graph: Store the graph after it is generated? (Default False)

output_name: Name of graph (if stored). (Default = './flux.png')


### pulse_gif(sim_file = './sim_file.csv',output_name = './output.gif')

sim_file: What previously defined input file would you like to use to run your analysis (see notes above).

output_name: Name of the generated pulse gif. (Default = './output.gif')


### fitting_gif(sim_time,sim_file = './sim_file.csv',x_scale='',y_scale='',outputName='./flux.png')

sim_time: The time (in seconds) you would like the simulation to run for. Values can be any positive float or interger.

sim_file: What previously defined input file would you like to use to run your analysis (see notes above).

x_scale: How would you like the x_axis (time) to be presented? Options are '' and 'log'. 

y_scale: How would you like the y_axis (outlet flux) to be presented? Options are '' and 'normalized'.

outputName: What would you like the name of the generated gif to be? Default is './flux.png'
