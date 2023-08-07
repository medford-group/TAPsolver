#from structures import exp
import copy
import jsonpickle
import json
import sys

def read_transient_sensitivity(file_name):

	with open(file_name) as f:
		data = json.loads(f.read())
	f.close()

	sameObject = jsonpickle.decode(data)
	sameObject2 = sameObject["1"]

	new_object = {}

	for k in sameObject2:
		new_object[k] = {}
		for j in sameObject2[k]:
			new_object[k][int(j)] = {}
			for z in sameObject2[k][j]:
				new_object[k][int(j)][z] = sameObject2[k][j][z]
	
	return new_object