from .mechanism import mechanism
#from structures import mechanism
import jsonpickle
import json
import sys

def read_mechanism_object(file_name):

	loaded_mechanism = mechanism()

	with open(file_name) as f:
		data = json.loads(f.read())
	f.close()

	sameObject = jsonpickle.decode(data)
	sameObject2 = sameObject["1"]

	loaded_mechanism.rate_array = sameObject2.rate_array
	loaded_mechanism.reactants = sameObject2.reactants 
	loaded_mechanism.reactions = sameObject2.reactions

	for jnum,j in enumerate(sameObject2.elementary_processes):
		loaded_mechanism.elementary_processes[jnum] = sameObject2.elementary_processes[j] 
		
	return loaded_mechanism