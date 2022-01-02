
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class experimental_data():

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