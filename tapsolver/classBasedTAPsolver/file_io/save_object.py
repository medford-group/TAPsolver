
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

import jsonpickle
import json
import sys

def save_object(storing_object,file_name):


	"""

	This function simplifies the storing of TAP objects for later use in simulation
	and analysis of TAP experiments and kinetic models. Makes sharing results less 
	complex and is offers an alternative to the traditional specification of the csv
	file.

	Args:

		storing_object (TAP structure): The TAP object/structure/class that the
		user wants to save and use later.

		file_name (str): The directory path and file name to store the json
		structure under. 
	
	Returns:

		Saved, encoded object in user specified directory  

	Implementor:

		Adam Yonge

	Link:
		
		https://jsonpickle.github.io/
	
	"""

	encoded_object = jsonpickle.encode({1: storing_object})
	
	with open(file_name, 'w', encoding='utf-8') as f:
		json.dump(encoded_object, f, ensure_ascii=False, indent=4)
	f.close()
