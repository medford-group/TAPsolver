

# Running a TAPsolver simulation (and associated parameters)

Now that the TAPsolver object has been defined, the user can run a simulation using TAPsolver. The barebones function is straightforward and can be called with:

	forward_problem(2,1,new_TAPobject)

where the variables called in the function represent the simulation time (seconds), the number of pulses and the newly defined TAPobject, respectively. Calling the python script where this line is defined will run the TAPsolver simulation associated with the mechanism, reactor and initial conditions. 

There are some additional useful features when calling the forward problem that can help sort data storage. If for some reason the user doesn't wish to store the outlet flux data, they can simply adjust the TAPobject by calling:

	new_TAPobject.store_flux_data = False

Another option that can be useful (but also eat up data space) is clarifying if the catalyst zone data should be stored during/after the simulation. This can be done with:

	new_TAPobject.store_catalyst_data = True

allowing the user to store the gas and adspecies concentrations within the catalyst zone. This can be highly useful for additional method development (i.e. applying KINNs or RRM to TAP data). Last, the user can turn all of the noise for gas and surface data on or off using:

	new_TAPobject.gas_noise = False 

The location of the data storage can also be adjusted with the command:

	new_TAPobject.output_name = 'exp_new'

where exp_new can be replaced with any desired output folder name. 

# Visualizing TAPsolver simulations (and associated parameters)

It is always useful for the user to have the ability to observe their own simulations (as well as there relationship to the experimental data). The barebones function to graph the outlet flux data is:

	flux_graph(new_TAPobject)

This will read the TAPobject and find the associated synthetic data generated during the simulation. Just like the forward_problem function, there are additional variables that can help determine what is and isn't shown during the graphing process. For example, the analytical diffusion for all gas species in the TAP simulation can be visualized with the command:

	new_TAPobject.display_analytical = True

The (real) experimental data can also be visualized by including or excluding a reference folder for it:

	new_TAPobject.data_name = 'experimental_data'

Please go to the next tutorial to learn more about the storing and reading TAPobjects, as well as the data that is generated during the simulation.
