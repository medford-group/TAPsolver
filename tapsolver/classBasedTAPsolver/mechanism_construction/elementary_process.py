
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class elementary_process():

	"""
	
	This class provides all kinetic information for an elementary process.
	
	Args:
		
		elementary_processes (int): The number of cells in the mesh used to solve the PDEs.

		rate_array (numpy array): The stoichiometric matrix of the microkinetic model.

		reactions (list): The individual elementary processes introduced in the mechanism (also found in elementary_processes).

		reactants (list): The reactive gas and surface species in the mechanism.   

	"""

	def __init__(self):
		self.processString = ''
		self.forward = {}
		self.backward = {}
