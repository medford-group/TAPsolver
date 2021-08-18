
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class experiment():

	"""
	
	This class contains all of the experiment information. It mainly acts as a container for all of the species data / methods as well as a container for meta data.


	Args:
		file_name (str): The name of the read file or experiment.

		collection_time (str): The collection time of the experiment.

		pulse_spacing (float): The end time of the experiment.

		num_samples_per_pulse (int): The number of samples per pulse.

		species_data (dict): A collection of Transient objects.

		reactor (class Reactor): The reactor information.


	"""

	def __init__(self):
		self.data = {}
		self.dataLocation = ''

	def readCSVInput(self, fileName, gasses):
		print('READ IN EXPERIMENTAL DATA')
		print('')
		data = pd.read_csv(fileName,header=None)

		rows_1, cols_1 = np.where(data == 'Reactor_Information')
		rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(data == 'Reaction_Information')

		reactor_info = data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
		self.dataLocation = reactor_info.iloc[7,1]	
		if self.dataLocation == 'none':
			print('')
			print('NO EXPERIMENTAL DATA AVAILABLE/PROVIDED')
			print('')
			print('')
		else:
			for k in gasses.keys():
				print(k)
				try:
					dfnew = pd.read_csv(self.dataLocation+'_folder/flux_data/'+k+'.csv',header=None)
					self.data['time'] = dfnew[:,0]
					self.data[k] = dfnew[:,1]
					print('DATA STORED AND FOUND FOR '+k)
					print('')
				except:
					print('NO EXPERIMENTAL DATA AVAILABLE FOR '+k)
		print('')
		print('')