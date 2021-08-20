
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class mechanism():

	"""
	
	This class acts as a container for all of the elementary processes (mechanism) and kinetic parameters
	
	Args:
		
		elementary_processes (int): The number of cells in the mesh used to solve the PDEs.

		catalyst_mesh_density (int): The number of times the number of mesh cells should be doubled in the catalyst zone.

		temperature (float): The temperature of the reactor (in K).

		output_name (str): The name of the output folder for the results of the analysis. 

		inert_diffusion (float): The reference value of the diffusion coefficient in the inert zone of the reactor.

		catalyst_diffusion (float): The reference value of the diffusion coefficient in the catalyst zone of the reactor. 

		reference_temperature (float): The temperature (in K) that the inert_diff and inert_cat values were determined.

		reference_mass (float): The mass (in a.m.u.) that inert_diff and inert_cat values were determined.

		advection (float): The rate of advection (velocity) of the gas species. One advection value for all gasses.


	"""

	def __init__(self):
		
		self.elementaryProcesses = {}
		self.rate_array = None
		self.reactions = None
		self.reactants = None
		self.kinetics = 'k'
