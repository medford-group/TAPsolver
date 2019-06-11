# Carbon Moxide Adsorption TAP Simulation

How to generate synthetic data for CO adsorption on a catalyst surface is outlined below. Adjustments to multiple pulses and multiple active sites are also discussed. The first two sections should be straight forward, but the later two could be confusing.

## Defining the Reactor

Basic reactor parameters must be defined by the user, including the 'reactor length', 'radius', 'catalyst zone fraction', 'inert/catalyst void fraction' and 'reactor temperature'. The user must also define parameters that control the precision of the simulator, including time steps and mesh size. Values of 500 for each should be acceptable for simple reactions (like CO adsorption) and reasonable pulse durations.

The diffusion for each species is determined through scaling to some reference species. If it is known, the user can enter the diffusion information in 'Reference Temperature', 'Reference Diffusion Inert', 'Reference Diffusion Catalyst' and 'Reference Temperature rows'.

Each time the simulation is run, a folder will be generated to store the results. The name of this folder (and its location) can be specified by 'Output Folder Name'. Similarly, if you plan on processing some experimental data (or previously generated synthetic data), you can specify its location by the 'Experimental Data Folder' variable. 

![Reactor Information](../../../figures/Reac_information.png)

*Section of the input_file.csv file that details the reactor setup and some simulation details*


## Data Storage and Simulation Options

It is possible to specify what processes are desired during a particular simulation and what should be stored after these processes have occured. This can be done in the 'Data Storage Options' section, where the user can enter 'TRUE' or 'FALSE' for each of the options. For example, if the user would like to perform a sensitivity analysis, but try to fit parameters, they would set 'Sensitivity Analysis' to 'TRUE' and 'Fit Parameters' to 'FALSE'. Details for the remaining options are listed [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file). 

![Reactor Information](../../../figures/data_storage.png)
*Section of the input_file.csv file that details data storage options*

## Feed & Surface Composition / Reaction Information 

This section is where things can be a little more confusing. In particular, listing the elementary reactions in the micro-kinetic model and specifying inputs like 'Pulse Ratio', 'Pulse Time', 'Mass List' and 'Initial Surface Composition'.

### Entering Elementary Reactions

Elementary Reactions (E.R.) follow a general trend:

- Species_1 + Species_2 -> Species_3 + Species_4

Where Species # can be any gas or surface species desired. There are two options for the arrow: either forward or reversible

- ->    : Irreversible E.R. 
- <->   : Reversible E.R.

The two cells directly to the right of each elementary reaction in the input_file.csv file allow the user to specify the forward and backward rate constants, with the backward rate constan always being listed after the forward rate constant. 

An example the reaction information input section is presented below.

![Reactor Information](../../../figures/reac_info.png)

It should be noted that this structure must be maintained for every elementary reaction. So no replacements can be made for the '+' or '->'.

### Gas and Surface Details

![Reactor Information](../../../figures/Feed_surf.png)

## Multi-pulse

## Multi-site


## Pump-Probe


## Additional Reactions


## Sensitivity Analysis


## Parameter Fitting



