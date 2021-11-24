from structures import exp
import jsonpickle
import json
import sys

def read_reactor_object(file_name):

	loaded_reactor = reactor()

	with open(file_name) as f:
		data = json.loads(f.read())
	f.close()

	sameObject = jsonpickle.decode(data)
	sameObject2 = sameObject["1"]

	loaded_reactor.zone_lengths = {0:sameObject2.zone_lengths['0'], 1: sameObject2.zone_lengths['1'], 2: sameObject2.zone_lengths['2']}
	loaded_reactor.zone_voids = {0:sameObject2.zone_voids['0'], 1: sameObject2.zone_voids['2'], 2: sameObject2.zone_voids['2']}

	loaded_reactor.reactor_radius = sameObject2.reactor_radius	
	loaded_reactor.temperature = sameObject2.temperature	

	return loaded_reactor