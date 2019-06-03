import numpy as np
import pandas as pd
import sys
from sklearn.preprocessing import scale
from scipy import stats
from scipy.spatial import distance_matrix
from scipy.cluster.hierarchy import dendrogram, linkage
from reaction_list_odes import display_odes, reac_list_parsing, find_partial, find_double_partial, Langmuir_rates

import time

file_location = "../ross_R_code/"

names_of_molecules = ['N','[HH]','N#N']

###Generate the mechanism###

reactions = Langmuir_rates(names_of_molecules)

############################

###Read in Data#############

data_combinations = []

types_of_things = ['g','s']

experimental_data = {}

desired_column = 364

y_column = {}

for k,i in enumerate(names_of_molecules):
	new_name = i+types_of_things[0]
	data_combinations.append(new_name)
	try:
		temp = pd.read_csv(file_location+i+"GasCon.csv")
		experimental_data[new_name] = temp.ix[:,desired_column]
		temp = pd.read_csv(file_location+i+"ReactionRate.csv")
		y_column['y'] = temp.ix[:,desired_column]
	except OSError as e:
		#print("No gas phase concentration data, using rate for species in concentrations place")
		#print("")
		temp = pd.read_csv(file_location+i+"ReactionRate.csv")
		experimental_data[new_name] = temp.ix[:,desired_column]

for k,i in enumerate(names_of_molecules):
	new_name = i+types_of_things[1]
	data_combinations.append(new_name)
	try:
		temp = pd.read_csv(file_location+i+"SurfaceCon.csv")
		experimental_data[new_name] = temp.ix[:,desired_column]
	except:
		print("Necessary files not found in current folder. May need to gather more data, rename files or find the correct directory")	

for k,i in enumerate(names_of_molecules):
	for j,n in enumerate(names_of_molecules):
		name1 = i+types_of_things[0]
		name2 = n+types_of_things[1]
		name3 = name1+name2
		data_combinations.append(name3)
		experimental_data[name3] = experimental_data[name1]*experimental_data[name2]


for k,i in enumerate(names_of_molecules):
	#for j in range(0,(len(names_of_molecules)-k)):
	for j in range(k,len(names_of_molecules)):
		name1 = i+types_of_things[1]
		name2 = names_of_molecules[j]+types_of_things[1]
		name3 = name1+name2
		data_combinations.append(name3)
		experimental_data[name3] = experimental_data[name1]*experimental_data[name2]

frames = []

for i,key in enumerate(experimental_data):
	frames.append(experimental_data[key])

new_df = pd.concat(frames,axis=1)

############################



###Construct odes###########

rate_array, reactions, reactants = reac_list_parsing(reactions)

############################

#print(rate_array)
#print(reactions)
#print(reactants)
#print(rate_array[0])
#print(display_odes(reactants[0],rate_array,reactions,reactants))

for j in range(0,3):
	#print('d(R_'+ reactants[j] + ')/' + 'd('+reactants[k]+') = '+find_partial(reactants[j],reactants[k],rate_array,reactions,reactants))
	print('d('+reactants[j]+')/dt = '+display_odes(reactants[j],rate_array,reactions,reactants))
print("")
for j in range(0,3):
	for k in range(0,6):
		print('d(R_'+ reactants[j] + ')/' + 'd('+reactants[k]+') = '+find_partial(reactants[j],reactants[k],rate_array,reactions,reactants))
sys.exit()

###perform analysis (ROSS)##

x_scale = scale(new_df)
y_scale = scale(y_column['y'])
x_scale_pd = pd.DataFrame(x_scale)
y_scale_pd = pd.DataFrame(y_scale)

distanceMat = pd.DataFrame(1 - x_scale_pd.corr())
#print(distanceMat)
Z = linkage(distanceMat,'complete')

print(Z)

############################


#Need to establish way to make sure I'm storing the correct y rate
#sys.exit()	


#More code below, do not delete sys.exit()















###Unnecessary Notes#######

#To get just the column names
test = list(table1.columns.values)
#To get the whole data frame
test = list(table1)
print(type(table1))
#ix seems like the simplest way to get this
test = list(table1.ix[:,0])
print(type(test))
test.pop(0)

for k in test:
	print(k)

#Now I have to parse the info based on the gas and surface components

#Get the script to show how many gas species
test2 = []
for k in test:
	if k[-1:] == 'g':
		test2.append(k[:-1])
	else:
		break

print(test2)

#Just to keep track of time
#start_time = time.time()
#print(time.time() - start_time)

###########################
