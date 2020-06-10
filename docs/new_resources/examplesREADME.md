
* [Overview](https://github.com/medford-group/TAPsolver/tree/master)
* [Installation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/installation)
* [Interface options](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/interfaceOptions)
* [Documentation](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file)
* [Questions & Development](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/questionsDiscussion)

# Tutorial + Examples

Two examples are presented: one for carbon monoxide adsorption and carbon monoxide oxidation. These give the user an appropriate understanding of TAPsolvers functions and how to properly implement them. It should be easy for the user to implement their own, more complex simulations after following these two examples.

Some details are glossed over, but can be found [here](https://github.com/medford-group/TAPsolver/tree/master/docs/resources/input_file).

## Getting started:

All available functions written in TAPsolver are accesable through the importing of the tapsolver package as follows:

	from tapsolver import *

To run simulations and analyze kinetic data using TAPsolver, the input file must be constructed. The input file is primarily defined through two components: the reactor and the microkinetic model. 

Users can easily generate a standard reactor file with the following command:

	define_reactor(reactor_name='./exampleReactor.csv')

This will provide standard input values which can be manually altered. Reactor files can easily be recycled for future experiments and multiple can be generated before analysis.

Next, a set of elementary reactions must manually be defined in a csv file. The list of elementary reactions are entered in each row of the first column. Examples of how to define the input file can be found here. 

Once both the reactor and microkinetic model are defined, the simulation input file can be constructed by running the following command:

	input_construction(reactor_name='.exampleReactor.csv',reaction_name='./exampleName.csv',new_name='./input_file.csv')

with the parameter 'new name' indicating the input file name. Once this is constructed, the only remaining simulation parameters required are 

If you would like to generate these files in bulk or run using high-throughput computing, a dictionary defining the mechanisms and parameter variations can be written and called using the following command:

	bulk_structure(reactor_setup='./exampleReactor.csv',iterDict=variationDictionaries,baseName='testGeneration')

with iterDict and baseName representing the dictionary of values and directory name for the files to be stored. This dictionary takes in values like mechanisms, pulse intensities, pulse timing and initial surface composition. An example of how to define the dictionary for an oxygen scrambling problem is as follows:

	mechs = {'mech1':'./o2_mech_1.csv','mech2':'./o2_mech_2.csv','mech3':'./o2_mech_3.csv'}

	inten = {'O218':[0.9,1]}

	tim = {'O218':[1e3,7e3,1e4]}

	surf = {'*':[1e3,7e3,1e4]}

	variationDictionaries = {}

	variationDictionaries['mechanisms'] = mechs

	variationDictionaries['intensity'] = inten

	variationDictionaries['time'] = tim

	variationDictionaries['surface'] = surf 

Although this might seem like an unnecessary number of steps to define the input file for a single analysis, I would like to highlight how flexible this truely is. The frustration of working with new code is frequently based on a limited understanding or unnecessarily rigid options. With this input format, the microkinetic model is automatically constructed and filled with user values. The most challenging task is reduced to trying to identify the correct microkinetic model rather than ensuring it is defined correctly. 

## Running a simulation (the forward problem)

Once the input file is defined, it is somewhat straight forward to run a simulation. The only two required inputs are the duration of the simulation (how long are you monitoring the gas) and the input file you're simulating. 

	run_tapsolver(timeFunc=1,inputFile='./input_file.csv')

Additional parameters can be defined in the function:

	pulseNumber = 1 # How many pulses would you like to simulate?

	includeNoise = 'TRUE' # Would you like to include some white noise in the synthetic data?

	store_thin_func = 'FALSE' # Would you like to store the surface compositionover time?

	store_flux_func = 'TRUE' # Would you like to store the outlet flux data?


The results of the simulation should can be found in the output folder name (defined in the input file).


## Analyzing the data (the inverse problem)

There are three primary functions available for the analysis of transient data using TAPsolver: sensitivity analysis, parameter fitting and uncertainty quantification. These can be performed with the following aptly named functions:

	run_sensitivity(timeFunc=0.1, inputFile='./input_file.csv')

	fit_parameters(timeFunc=0.1, inputFile='./input_file.csv')

	run_uncertainty(timeFunc=0.1, inputFile='./input_file.csv')

Two notes should be made. First, the run sensitivity command has two modes of operation: overall (what is used to fit kinetic parameters to the experimental data) and transient (how the sensitivity of each kinetic parameter varies over time). These two options each serve their own useful purpose, but the function will default to the overall sensitivity. Second, the run uncertainty command provides the user with the value of the Hessian for the optimization process. As is noted in the literature (include source for the uq from hessian), an accurate value of the Hessian can accurately define the uncertainty of the fitted kinetic parameters.

## Visualizing results

Last, and potentially the most important component, is how to easily visualize your results. The current visualization options in TAPsolver focus on the simulation and fitting results, although options for visualizing the transient sensitivity have been developed in a non-general format. To view the results of the simulation, run the command:

	fluxgraph(timeFunc=1,input_file='./input_file.csv',pulse=1)

The final value pulse indicates which pulses you would like to visualize and can be defined as an integer representing the pulse number or a list of all the pulses you would like to view together (i.e. [1,5,11,20]). 

You can also generate a gif showing each pulse in your simulation in series with the command:

	pulsingVisualization(input_file='./input_file.csv',fileName='./output.gif')

The convergence of the parameter fitting routine can also be visualized with the command:

	fitting_gif(timeFunc=0.1,input_file='./input_file.csv')
	
