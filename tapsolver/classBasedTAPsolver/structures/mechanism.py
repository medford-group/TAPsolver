
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class mechanism():

	"""
	
	This class acts as a container for all of the elementary processes (mechanism) and kinetic parameters
	
	Args:
		
		elementary_processes (int): The number of cells in the mesh used to solve the PDEs.

		rate_array (numpy array): The stoichiometric matrix of the microkinetic model.

		reactions (list): The individual elementary processes introduced in the mechanism (also found in elementary_processes).

		reactants (list): The reactive gas and surface species in the mechanism.   

	"""

	def __init__(self):
		
		self.elementary_processes = {}
		self.rate_array = None
		self.reactions = None
		self.reactants = None