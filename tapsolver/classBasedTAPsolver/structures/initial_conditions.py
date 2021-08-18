
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import pandas as pd
import numpy as np

class initial_conditions():

	def __init__(self):
		self.gasses = {}
		self.inertGasses = {}
		self.surfaceSpecies = {}

	def showGasses(self):
		print('REACTION/PRODUCT GASSES')
		for k in self.gasses.keys():
			print('_________________')
			print('|'+k)
			print('|'+'Intensity: '+str(self.gasses[k].intensity))
			print('|'+'Delay: '+str(self.gasses[k].delay))
			print('|'+'mass: '+str(self.gasses[k].mass))#['intensity']
			print('~~~~~~~~~~~~~~~~~')
		if len(self.gasses) == 0:
			print('_________________')
			print('| None')
			print('~~~~~~~~~~~~~~~~~')
		print('')
		print('')

		print('INERT GASSES')
		for k in self.inertGasses.keys():
			print('_________________')
			print('|'+k)
			print('|Intensity: '+str(self.inertGasses[k].intensity))
			print('|Delay: '+str(self.inertGasses[k].delay))
			print('|mass: '+str(self.inertGasses[k].mass))#['intensity']
			print('~~~~~~~~~~~~~~~~~')

		if len(self.inertGasses) == 0:
			print('_________________')
			print('| None')
			print('~~~~~~~~~~~~~~~~~')

		print('')
		print('')

	def showSurface(self):
		print('REACTION/PRODUCT GASSES')
		for k in self.surfaceSpecies.keys():
			print('_________________')
			print('|'+k)
			print('|Conecentration: '+str(self.surfaceSpecies[k].concentration))
			print('~~~~~~~~~~~~~~~~~')
		if len(self.surfaceSpecies) == 0:
			print('_________________')
			print('| None')
			print('~~~~~~~~~~~~~~~~~')
		print('')
		print('')

	def readCSVInput(self, fileName):

		data = pd.read_csv(fileName,header=None)
		rows_1, cols_1 = np.where(data == 'Reactor_Information')
		rows_2, cols_2 = np.where(data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(data == 'Reaction_Information')

		feed_surf_info = data.iloc[1+rows_2[0]:rows_4[0]-1,:]

		gas_species = [x for x in feed_surf_info.iloc[0,1:] if type(x) == str]
		surface_species = [x for x in feed_surf_info.iloc[5,1:] if type(x) == str]

		for jnum,j in enumerate(gas_species):


			class specificGas():
				def __init__(self):
					self.mass = 0
					self.intensity = 0
					self.delay = 0


			class specificAdspecies():
				def __init__(self):
					self.concentration = 0

			if j.find('Inert') == 0:
				self.inertGasses[j] = specificGas()
				self.inertgasses[j].intensity = float(feed_surf_info.iloc[1,1+jnum])
				self.inertgasses[j].delay = float(feed_surf_info.iloc[2,1+jnum])
				self.inertgasses[j].mass = float(feed_surf_info.iloc[3,1+jnum])

			else:
				self.gasses[j] = specificGas()
				self.gasses[j].intensity = float(feed_surf_info.iloc[1,1+jnum])
				self.gasses[j].delay = float(feed_surf_info.iloc[2,1+jnum])
				self.gasses[j].mass = float(feed_surf_info.iloc[3,1+jnum])


		for jnum,j in enumerate(surface_species):
			self.surfaceSpecies[j] = specificAdspecies()
			if jnum == len(surface_species)-1:
				self.surfaceSpecies[j].conectration = feed_surf_info.iloc[6,1+jnum]