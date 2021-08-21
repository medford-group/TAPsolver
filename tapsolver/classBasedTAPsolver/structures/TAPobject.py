
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class TAPobject():

	"""
	
	This class acts as a container for the three necessary components to run a TAP simulation: the reactor details, mechanism and initial conditions.
	
	Args:
		
		reactor (reactor object): The reactor specifications.

		mechanism (mechanism object): The elementary processes and their associated kinetics.

		initial_conditions (initial_conditions object): The initial conditions of the TAP reactor/experiment.

	"""


	def __init__(self):
		
		self.reactor = reactor()
		self.mechanism = mechanism()
		self.initial_conditions = initial_conditions()