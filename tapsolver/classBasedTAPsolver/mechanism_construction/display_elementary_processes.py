
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np
from structures import mechanism

def displayProcesses(mechanism_data: mechanism):
	print('ELEMENTARY PROCESSES')
	
	for k in mechanism_data.elementaryProcesses.keys():
		print('__________________________________')
		print('|'+'Reaction '+str(k)+': '+mechanism_data.elementaryProcesses[k].processString)
		print('|')
		print('|'+'forward - ')
		print('|'+str(mechanism_data.elementaryProcesses[k].forward.k))
		print('|'+str(mechanism_data.elementaryProcesses[k].forward.Ao))
		print('|'+str(mechanism_data.elementaryProcesses[k].forward.Ea))
		print('|'+str(mechanism_data.elementaryProcesses[k].forward.Ga))
		print('|'+str(mechanism_data.elementaryProcesses[k].forward.dG))
		print('|'+str(mechanism_data.elementaryProcesses[k].forward.link))
		print('|')
		print('|'+'backward - ')
		print('|'+str(mechanism_data.elementaryProcesses[k].backward.k))
		print('|'+str(mechanism_data.elementaryProcesses[k].backward.Ao))
		print('|'+str(mechanism_data.elementaryProcesses[k].backward.Ea))
		print('|'+str(mechanism_data.elementaryProcesses[k].backward.Ga))
		print('|'+str(mechanism_data.elementaryProcesses[k].backward.dG))
		print('|'+str(mechanism_data.elementaryProcesses[k].backward.link))
		print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
	print('')
	print('')