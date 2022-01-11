#from structures import exp
import copy
import jsonpickle
import json
import sys

def read_experimental_data_object(file_name):

	with open(file_name) as f:
		data = json.loads(f.read())
	f.close()

	sameObject = jsonpickle.decode(data)
	sameObject2 = sameObject["1"]
	new_object = copy.deepcopy(sameObject2)
	for k in sameObject2:
		for j in sameObject2[k]:
			new_object[k][float(j)] = sameObject2[k][j]
			new_object[k][j].pop()
	
	return new_object