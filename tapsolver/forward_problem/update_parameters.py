from structures import *
import pandas as pd
import copy
import sys

def update_parameters(TAPobject_data_original: TAPobject):
	TAPobject_data = copy.deepcopy(TAPobject_data_original)
	new_addition = pd.read_csv('./'+TAPobject_data.output_name+'/optimization_results.csv')
	value_row = list(new_addition.loc[len(new_addition)-1])
	
	name_row = list(new_addition.columns)
	print(value_row)
	print(name_row)
	for jnum,j in enumerate(name_row):
		#print(j.split('['))
		#sys.exit()
		try:
			elementary_step = int(j.split('[')[1].split(']')[0])
			direction = j.split('].')[1].split('.')[0]
			if 'Ga' in j:
				if direction == 'forward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ga = value_row[jnum]
				elif direction == 'backward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ga = value_row[jnum]
			elif 'dG' in j:
				if direction == 'forward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].forward.dG = value_row[jnum]
				elif direction == 'backward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].backward.dG = value_row[jnum]
			elif 'Ao' in j:		
				if direction == 'forward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ao = value_row[jnum]
				elif direction == 'backward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ao = value_row[jnum]
			elif 'Ea' in j:		
				if direction == 'forward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].forward.Ea = value_row[jnum]
				elif direction == 'backward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].backward.Ea = value_row[jnum]
			elif 'k' in j:		
				if direction == 'forward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].forward.k = value_row[jnum]
				elif direction == 'backward':
					TAPobject_data.mechanism.elementary_processes[elementary_step].backward.k = value_row[jnum]

		except:
			pass
	return TAPobject_data