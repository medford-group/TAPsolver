import copy
import sys

def new_experiments(original_TAPobject,new_label,original_gasses=None,original_inert_gasses=None,original_adspecies=None,original_reactants=None):
	
	new_TAPobject = copy.deepcopy(original_TAPobject)
	
	if original_gasses == None:
		original_gasses = original_TAPobject.reactor_species.gasses	
	if original_inert_gasses == None:
		original_inert_gasses = original_TAPobject.reactor_species.inert_gasses
	if original_adspecies == None:
		original_adspecies = original_TAPobject.reactor_species.adspecies
	if original_reactants == None:
		original_reactants = original_TAPobject.mechanism.reactants
	
	if (type(original_TAPobject.reactor_species.temperature) == int) or (type(original_TAPobject.reactor_species.temperature) == float):
		new_TAPobject.reactor_species.temperature = {0:copy.deepcopy(original_TAPobject.reactor_species.temperature),1:copy.deepcopy(original_TAPobject.reactor_species.temperature)}
		newest_temperature = 1
	else:
		newest_temperature = len(original_TAPobject.reactor_species.temperature.keys())
		new_TAPobject.reactor_species.temperature[newest_temperature] = original_TAPobject.reactor_species.temperature[0]
	old_name_list = []
	new_name_list = []
	
	for k in original_gasses:
		old_name_list.append(k)	
		new_name_list.append(k+new_label)	
		temp_gas = copy.deepcopy(original_TAPobject.reactor_species.gasses[k])
		temp_gas.temperature_used = newest_temperature
		new_TAPobject.reactor_species.add_gas(k+new_label,temp_gas) 	
	
	for k in original_inert_gasses:
		temp_inert_gas = copy.deepcopy(original_TAPobject.reactor_species.inert_gasses[k])
		temp_inert_gas.temperature_used = newest_temperature
		new_TAPobject.reactor_species.add_inert_gas(k+new_label,temp_inert_gas)
	
	for k in original_adspecies:
		old_name_list.append(k)
		new_name_list.append(k+new_label)
		temp_adspecies = copy.deepcopy(original_TAPobject.reactor_species.adspecies[k])
		new_TAPobject.reactor_species.add_adspecies(k+new_label,temp_adspecies)
	
	if original_TAPobject.mechanism.reactants != None:
		if original_TAPobject.mechanism.reactants == list:
			temp_reactants = copy.deepcopy(original_TAPobject.mechanism.reactants[0])
			for knum,k in enumerate(old_name_list):
				try:
					index1 = temp_reactants.index(k)
					temp_reactants[index1] = new_name_list[knum]
				except:
					pass
			new_TAPobject.append(new_name_list)
		else:
			new_TAPobject.mechanism.reactants = [copy.deepcopy(original_TAPobject.mechanism.reactants)]
			if newest_temperature >= 2:
				temp_reactants = copy.deepcopy(original_TAPobject.mechanism.reactants[0])
			else:
				temp_reactants = copy.deepcopy(original_TAPobject.mechanism.reactants)
			
			for knum,k in enumerate(old_name_list):
				try:
					index1 = temp_reactants.index(k)
					temp_reactants[index1] = new_name_list[knum]
				except:
					pass
			if newest_temperature >= 2:
				new_TAPobject.mechanism.reactants[0].append(temp_reactants)
				new_TAPobject.mechanism.reactants = new_TAPobject.mechanism.reactants[0]
			else:
				new_TAPobject.mechanism.reactants.append(temp_reactants)
				
	return new_TAPobject