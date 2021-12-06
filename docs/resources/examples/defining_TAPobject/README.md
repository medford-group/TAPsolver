
# Generating a TAPobject

Let's first import tapsolver

	from tapsolver import *

TAP experiments involve a series of gas species diffusing through a packed bed reactor, where they can interact with a reactive material to form adspecies. It is necessary to define the reactor and the molecular species in the reactor before analysis of experimental data can be performed. 

A new reactor in TAPsolver can be defined through:

	new_reactor = reactor()

This automatically defines the dimensional parameters of the reactor. The default settings of these parameters are:

	zone_lengths = {0: 3, 1: 0.06, 2: 3}

	zone_voids = {0: 0.4, 1: 0.4, 2: 0.4}

	reactor_radius = 1.0

where the base units of zone length and reactor radius are cm. These parameter can be adjusted to the reactor of interest with:

	new_reactor.zone_lengths = {0: 2, 1: 0.1, 2: 2}
	
	new_reactor.zone_voids = {0: 0.5, 1: 0.5, 2: 0.5}

	new_reactor.reactor_radius = 0.5

Next, the gas and adspecies can be defined. First, the object where all gas and adspecies will be stored must be specified with:

	new_reactor_species = reactor_species()

Like the reactor, the reactor species class has default values of the parameters the might need to be adjusted, including:

	inert_diffusion = 16

	catalyst_diffusion = 16

	reference_temperature = 385.6

	reference_mass = 40

	temperature = 385.6

These parameters can be adjusted in a similar way to that of the reactor object and are used to specify the transport within the reactor. The next step is to add the individual gasses of interest. This process involves defining a new gas object (different from the reactor species object) and adjust the object variables so they are consitent with the gas. For example, you could define carbon monoxide and add it to the reactor species object through:

	CO = define_gas()

	CO.mass = 28

	CO.intensity = 1

	CO.delay = 0.0

	new_reactor_species.add_gas('CO',CO)

where the mass is in amu, intensity is the total number of molecules (nmol) being introduced to the system, and delay is in seconds (at what time are the gasses being pulsed). The add gas function incorporates the gas into the reactor species object. If the user would like to add an inert gas species, a similar approach is taken:

	argon = define_gas()

	argon.mass = 40

	argon.intensity = 1

	argon.delay = 0.0

	new_reactor_species.add_inert_gas('argon',argon)

With the only difference being the add inert gas function in place of the add gas. The last step is to add any adspecies to the system and follows a similar formula to that of the gasses. 

	s1 = define_adspecies()

	s1.concentration = 0

	new_reactor_species.add_adspecies('CO*',s1)

	s2 = define_adspecies()

	s2.concentration = 30

	new_reactor_species.add_adspecies('*', s2)

Concentration is the density of species in the catalyst zone of the reactor (the initial CO* concentration is 0 nmol/cm3, while the active site is 30 nmol/cm3).

Now that we've defined the reactor and the associated gas and adspecies, we can begin to define the mechanism and kinetic details of the system (if any). Since we've defined CO, CO* and \*, we can just use carbon monoxide adsorption as our example. First, we generate a mechanism object

	new_mechanism = mechanism()

The variables in the mechanism are mostly used internally and there are few benefits for the user to adjust them. The only process necessary for the user is adding individual elementary processes and kinetic parameters (or more generally reaction steps) to the system. This is accomplished through:

	new_mechanism.elementary_processes[0] = elementary_process('CO + * <-> CO*')

where the 0 is the first entry into the micro-kinetic model and 'CO + * <-> CO\*' is the explicit elementary process (adsorption of carbon monoxide on an active site). The additional active sites can be added as:

	new_mechanism.elementary_processes[1] = elementary_process('O2 + 2* <-> 2O*')
	
	new_mechanism.elementary_processes[2] = elementary_process('CO* + O* <-> 2* + CO2')

The kinetic parameters can then be specified with the following commands:

	new_mechanism.elementary_processes[0].forward.k = 1
	
	new_mechanism.elementary_processes[0].backward.k = 1
	
	new_mechanism.elementary_processes[1].forward.k = 1
	
	new_mechanism.elementary_processes[1].backward.k = 1
	
	new_mechanism.elementary_processes[2].forward.k = 1
	
	new_mechanism.elementary_processes[2].backward.k = 1 

Where a forward and backward value for each of the elementary processes is specified (only 1 here). Simulations can be run using the free energy specifications, as well as the pre-exponential/activation energies. An additional tutorial outlining this can be found in the examples. The final step in making your mechanism is constructing the stoichiometric matrix with:

	mechanism_constructor(new_mechanism)

Finally, the TAPobject can be generated.

	TAP_test = TAPobject()

Where the subobject (mechanism, species, reactor) can be defined here:

	TAP_test.mechanism = new_mechanism
	
	TAP_test.reactor_species = new_reactor_species
	
	TAP_test.reactor = testGen1