
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from .mechanism import mechanism
#from structures import mechanism

def display_processes(mechanism_data: mechanism):
	
	"""
	
	Function to clearly define all of the parameters 
	
	Args:
		
		elementary_processes (int): The number of cells in the mesh used to solve the PDEs.

		rate_array (numpy array): The stoichiometric matrix of the microkinetic model.

		reactions (list): The individual elementary processes introduced in the mechanism (also found in elementary_processes).

		reactants (list): The reactive gas and surface species in the mechanism.   

	"""

	print('ELEMENTARY PROCESSES')
	
	for k in mechanism_data.processes.keys():
		print('__________________________________')
		print('|'+'Reaction '+str(k)+': '+mechanism_data.elementary_processes[k].processString)
		print('|')
		print('|'+'forward - ')
		print('| k: '+str(mechanism_data.processes[k].f.k))
		print('| Ao: '+str(mechanism_data.processes[k].f.Ao))
		print('| Ea: '+str(mechanism_data.processes[k].f.Ea))
		print('| Ga: '+str(mechanism_data.processes[k].f.Ga))
		print('| dG: '+str(mechanism_data.processes[k].f.dG))
		#print('| link: '+str(mechanism_data.elementary_processes[k].forward.link))
		print('|')
		print('|'+'backward - ')
		print('| k: '+str(mechanism_data.processes[k].b.k))
		print('| Ao: '+str(mechanism_data.processes[k].b.Ao))
		print('| Ea: '+str(mechanism_data.processes[k].b.Ea))
		print('| Ga: '+str(mechanism_data.processes[k].b.Ga))
		print('| dG: '+str(mechanism_data.processes[k].b.dG))
		#print('| link: '+str(mechanism_data.elementary_processes[k].backward.link))
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	print('')
	print('')