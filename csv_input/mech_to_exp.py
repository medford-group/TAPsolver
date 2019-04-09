import numpy as np

def reac_list_parsing(reactions):

	reacs = reactions.copy()
	reactants = []

	#need to make this a function

	for k,i in enumerate(reacs):
		reacs[k] = reacs[k].replace('<->','')
		reacs[k] = reacs[k].replace('->','')
		reacs[k] = reacs[k].replace('+','')
		reacs[k] = reacs[k].split()

		for n,j in enumerate(reacs[k].copy()):
			c = j
			if j[0].isdigit():
				j = j[1:]
			if j in reactants:
				reacs[k].remove(c)
			else:
				reactants.append(j)

	insert_location = 0

	for k,i in enumerate(reactants):
		if '*' not in i:
			reactants.insert(insert_location, reactants.pop(k))
			insert_location += 1
		if i == '*':
			reactants.pop(k)

	#Add the surface site to the end of the list of reactants
	reactants.insert(len(reactants), '*')

	#Display all of the reactants and reactions
	#print("")
	#print("All elementary reactions considered")
	#print(reactions)
	#print("")
	#print("Gas phase molecules & surface intermediates")
	#print(reactants)
	#print("")

	#make empty matrix for storing coefficients in reactions
	rate_array = np.zeros((len(reactions),len(reactants)))

	#step through reactions and store stoichiometry
	for k,n in enumerate(reactions):
		reactions[k] = reactions[k].replace('+','')
		try:
			neg,pos = reactions[k].split('<->')
		except ValueError:
			neg,pos = reactions[k].split('->')

		new_neg = neg.split()
		new_pos = pos.split()

		for val in new_neg:
			into_array = -1
			new_val = val
			if val[0].isdigit():
				new_val = val[1:]
				into_array = into_array*int(val[0])
			rate_array[k,reactants.index(new_val)] = into_array

		for val in new_pos:
			into_array = 1
			new_val = val
			if val[0].isdigit():
				new_val = val[1:]
				into_array = into_array*int(val[0])
			rate_array[k,reactants.index(new_val)] = into_array

	return rate_array, reactions, reactants
	#display the generated rate_array
	#print(rate_array)

gasses = ['CO','O2','CO2']

mech_1 =['CO + * <-> CO*','O2 + 2* <-> 2O*','CO + O* <-> CO2*','CO2* <-> CO2 + *'] #Eley - Rideal
mech_2 =['O2 + 2* <-> 2O*','CO + * <-> CO*','CO* + O* <-> CO2*','CO2* <-> CO2 + *'] #Langmuir - Hinshelwood
mech_3 =['CO + * <-> CO*','O2 + 2* <-> 2O*','2CO* <-> C* + CO2*','C* + O* <-> CO*','CO2* <-> CO2 + *'] #CO / CO combinations

mechanisms = [mech_1,mech_2,mech_3]
mech_info = {}
for k in range(0,len(mechanisms)):
	mech_info[k] = {}
for k_num,k in enumerate(mechanisms):

	out_1,out_2,out_3 = reac_list_parsing(k)

	mech_info[k_num]['array'] = out_1
	mech_info[k_num]['reacs'] = out_2
	mech_info[k_num]['species'] = out_3

print(mech_info[0]['species'])
print(mech_info[1]['species'])
print(mech_info[2]['species'])

for k in range(0,len(mechanisms)):
	for j in range(0,len(gasses)):
		
