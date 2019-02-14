from numpy import inf 
import numpy as np
from math import isnan
from math import sqrt
import sympy
import time 
import sys

from imolecule.notebook import generate
from imolecule.format_converter import convert
import networkx.algorithms.isomorphism as iso
from networkx import is_isomorphic
import networkx as nx
#from rxn_network import scissions, recursive_scissions, get_mechanisms
#from data_structures import MolGraph, RxnGraph


###########################PUT IN LATEX FORMAT##############################

#print('\documentclass[a4paper]{article}')

#print('%% Language and font encodings')
#print('\\usepackage[english]{babel}')
#print('\\usepackage[utf8x]{inputenc}')
#print('\\usepackage[T1]{fontenc}')
#print('\\usepackage{tikz}')
#print('\\usepackage[super]{natbib}')
#print('\\usepackage{setspace}')
#print('%% Sets page size and margins')
#print('\\usepackage[a4paper,top=3cm,bottom=2cm,left=3cm,right=3cm,marginparwidth=1.75cm]{geometry}')

#print('%% Useful packages')
#print('\\usepackage{amsmath}')
#print('\\usepackage{graphicx}')
#print('\\usepackage[colorinlistoftodos]{todonotes}')
#print('\\usepackage[colorlinks=true, ')

#print('\\begin{document}')




############################################################################




#####################################################
#I've used this script to break down reaction lists into there components. Once a matrix is constructed, it can be fed to the
#F_constructor.py or fenics_simulator for other use. Originally, the goal was to root the reactivities generated in the RRM 
#back to intrinsic kinetics, but my priorities have shifted to simulation and parameter fitting.



#Test case for simple system of A, B and C with 6 possible reactions in the network
#reactions = ['A + * <-> A*', 'B + * <-> B*', 'C + * <-> 2C*', 'A* + B* <-> C* + 2*', 'A + B* <-> C', 'A + B* <-> C*']

#More specific: ammonia decomposition. Large number of reactants and reactions
#reactions = ['NH3 + * <-> NH3*', 'NH3 + 2* <-> NH2* + H*',     'H2 + * <-> H2*', 'H2 + 2* <-> 2H*',   'N2 + * <-> N2*', 'N2 + 2* <-> 2N*', 'N2H4 + * <-> N2H4*', 'N2H2 + * <-> N2H2*', \
#		'NH2* + * <-> NH* + H*','NH* + * <-> N* + H*',\
#		'N2H4* + * <-> 2N2H2*', 'N2H4* + * <-> N2H3* + H*', 'N2H3* + * <-> N2H2* + H*', 'N2H2* + * <-> N2H* + H*','N2H* + * <-> N2* + H*', 'N2H2* + * <-> 2NH2*', 'NH3* + H* <-> NH2* + H2*']

#ammonia decomposition generated from the reaction network generator
#reactions = ['H2 <-> H2', 'NN <-> NN', 'NHNH <-> NHNH', 'NH2NH2 <-> NH2NH2', 'NH3 <-> NH3', 'H2 <-> H + H', 'NN <-> N + N', 'NHNH <-> NH + NH', 'NHNH <-> NNH + H', 'NH2NH2 <-> NH2 + NH2', 'NH2NH2 <-> NH2NH + H', 'NH3 <-> NH2 + H', 'NH <-> N + H', 'NNH <-> NH + N', 'NNH <-> NN + H', 'NH2 <-> NH + H', 'NH2NH <-> NH2N + H', 'NH2NH <-> NHNH + H', 'NH2NH <-> NH + NH2', 'NH2N <-> NH2 + N', 'NH2N <-> NHN + H']

#Copy of reactions for ease in generating complete list of reactants

def variational_list_parsing(reactions):

	reacs = reactions.copy()
	reactants = []
	rev_irr = []

	#need to make this a function
	print(reacs)
	for k,i in enumerate(reacs):
		if '<->' in reacs[k]:
			rev_irr.append(1)
			reacs[k] = reacs[k].replace('<->','')
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
		else:
			rev_irr.append(0)
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

	#make empty matrix for storing coefficients in reactions
	rate_array = np.zeros((len(reactions),len(reactants)))

	#step through reactions and store stoichiometry
	for k,n in enumerate(reactions):
		reactions[k] = reactions[k].replace('+','')
		if '<->' in reactions[k]:
			neg,pos = reactions[k].split('<->')
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
		else:
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

	return rate_array, reactions, reactants, rev_irr

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

# rate_array, reactions, reactants, rev_irr = reac_list_parsing(reactions_n)

#rate_array, reactions, reactants = reac_list_parsing(reactions)

#function used to make sure that components appear on same side of elementary
#reactions
def test_value(tester1,tester2):
	tester1_truth = False
	if tester1 > 0:
		tester1_truth = True
	tester2_truth = False
	if tester2 > 0:
		tester2_truth = True
	if (tester1_truth == tester2_truth) and tester1!= 0 and tester2!=0:
		return True
	return False

def display_odes(num,rate_array,reactions,reactants):
	if type(num) != str:
		print("failed")
		return
	#store elementary reaction number involved	
	deriv_of_reacs = []
	#store coefficient associated with elementary reaction
	deriv_of_reacs_val = []
	#initiate reaction expression as empty string
	reaction_expression = ''
	for k,j in enumerate(reactions):
		if rate_array[k,reactants.index(num)] != 0:
			deriv_of_reacs.append(k)
			deriv_of_reacs_val.append(rate_array[k,reactants.index(num)])
	for step,k in enumerate(deriv_of_reacs):
		if reaction_expression != '':
			reaction_expression = reaction_expression + ' + ('+str(int(deriv_of_reacs_val[step]))+')*( '
		else:
			reaction_expression = reaction_expression + ' ('+str(int(deriv_of_reacs_val[step]))+')*( '
		rate_constant = "ke"+str(k)
		reaction_expression = reaction_expression + rate_constant
		for v,z in enumerate(reactants):
			if (test_value(-1,rate_array[k,v]) == True):
				#rate_constant = "k"+str(k+1)+"f"
				if abs(rate_array[k,v]) > 1:
					reaction_expression = reaction_expression + '*([' + reactants[v] + ']^'+str(int(abs(rate_array[k,v])))+')'
				else:
					reaction_expression = reaction_expression + '*[' + reactants[v] + ']'
		rate_constant = " - kd"+str(k)
		reaction_expression = reaction_expression + rate_constant
		for v,z in enumerate(reactants):
			if (test_value(1,rate_array[k,v]) == True):
				if abs(rate_array[k,v]) > 1:
					reaction_expression = reaction_expression + '*([' + reactants[v] + ']^'+str(int(abs(rate_array[k,v])))+')'
				else:
					reaction_expression = reaction_expression + '*[' + reactants[v] + ']'
		reaction_expression = reaction_expression + ' )'
	return reaction_expression #'d('+num+')/dt = '+

def find_partial(num,den,rate_array,reactions,reactants):
	#check if function input num and den are strings
	if type(num) != str and type(den) != str:
		print("failed")
		return
	#store elementary reaction number involved	
	deriv_of_reacs = []
	#store coefficient associated with elementary reaction
	deriv_of_reacs_val = []
	#initiate reaction expression as empty string
	reaction_expression = ''
	#find each elementary reaction involved
	for k,j in enumerate(reactions):
		if rate_array[k,reactants.index(num)] != 0:
			deriv_of_reacs.append(k)
			deriv_of_reacs_val.append(rate_array[k,reactants.index(num)])
	#for each 
	for step,k in enumerate(deriv_of_reacs):
		if rate_array[k,reactants.index(den)] != 0:
			if reaction_expression != '':
				reaction_expression = reaction_expression + ' + '
			reaction_expression = reaction_expression +str(int(deriv_of_reacs_val[step]))+'*('
			
			rate_constant = ''
			pass_test = False

			if rate_array[k,reactants.index(den)] < 0:
				rate_constant = "ke"+str(k)
				pass_test = True
			elif rate_array[k,reactants.index(den)] > 0 and '<->' in reactions[k]:  #!!!!
				rate_constant = "kd"+str(k)
				pass_test = True
			reaction_expression = reaction_expression + rate_constant
			if pass_test == True:
				for v,z in enumerate(reactants):	
					if (test_value(rate_array[k,reactants.index(den)],rate_array[k,v]) == True):
						if (reactants.index(den) != v):
							if abs(rate_array[k,v]) > 1:
								reaction_expression = reaction_expression + '*([' + reactants[v] + ']^'+str(int(abs(rate_array[k,v])))+')'
							else:
								reaction_expression = reaction_expression + '*[' + reactants[v] + ']'
						elif abs(int(rate_array[k,v])) == 2:
							reaction_expression = reaction_expression + '*(2*[' + reactants[v] + ']'


			reaction_expression = reaction_expression + ')'
	if reaction_expression == '':
		return '0'
	return(reaction_expression)

#Function used to evaluate the 2nd derivative of a specific combination of values
def find_double_partial(num,den,den2,rate_array,reactions,reactants):
	if type(num) != str and type(den) != str:
		print("failed")
		return
	deriv_of_reacs = []
	deriv_of_reacs_val = []
	reaction_expression = ''
	for k,j in enumerate(reactions):
		if rate_array[k,reactants.index(num)] != 0:
			deriv_of_reacs.append(k)
			deriv_of_reacs_val.append(rate_array[k,reactants.index(num)])
	for step,k in enumerate(deriv_of_reacs):
		if (den == den2) and abs(rate_array[k,reactants.index(den)]) != 2:
			pass
		elif test_value(rate_array[k,reactants.index(den)],rate_array[k,reactants.index(den2)]) == True:
			if reaction_expression != '':
				reaction_expression = reaction_expression + ' + '
			reaction_expression = reaction_expression +str(int(deriv_of_reacs_val[step]))+'*('
			if rate_array[k,reactants.index(den)] < 0:
				rate_constant = "ke"+str(k)
			else:
				rate_constant = "kd"+str(k)
			reaction_expression = reaction_expression + rate_constant
			for v,z in enumerate(reactants):
				#Need to check to see if the appropriate 2nd derivatives are being returned

				if test_value(rate_array[k,reactants.index(den)],rate_array[k,v]) == True:
					if (reactants.index(den) != v) and (reactants.index(den2) != v):
						if abs(rate_array[k,v]) > 1:
							reaction_expression = reaction_expression + '*([' + reactants[v] + ']^'+str(int(abs(rate_array[k,v])))+')'
						else:
							reaction_expression = reaction_expression + '*[' + reactants[v] + ']'
					elif abs(int(rate_array[k,v])) == 2 and (den != den2):
						reaction_expression = reaction_expression + '*(2*[' + reactants[v] + ']'
					else:
						reaction_expression = reaction_expression + '*2'

			reaction_expression = reaction_expression + ')'
	if reaction_expression == '':
		return ' 0'
	return(reaction_expression)



######################################################################################################
######################################################################################################
######################################################################################################
#Code below is from the front_end.py script. Didn't want to simplify code and have pulse simulator no 
#longer work
#Make graph molecule from smiles input and store in dict_of_reactants
def make_molegraph(name_in_smiles, dict_of_reactants):
	molecule_num = len(dict_of_reactants)
	new_name = 'A'+str(molecule_num)
	A = MolGraph()
	sto_test = A.generate(name_in_smiles,'smi')
	dict_of_reactants[new_name] = sto_test
	return dict_of_reactants

#Develope a rate expression based on langmuir hinshelwood kinetics
def Langmuir_rates(name_in_smiles):

	#Need to make the equations for the adsorption and desorption of the molecules.
	#This is where it is done.
	def generate_gas_equations(network_input):
		simp_ads = []
		names_of_inputs = []
		surface_reactions = []
		for j in network_input:
			the_name = j.__altname__()
			new_reac_name = the_name+' + * <-> '+the_name+'*'
			names_of_inputs.append(the_name+'*')
			simp_ads.append(new_reac_name)
		for i in range(0,len(rxn2)):
			surface_reactions.append(str(rxn2[i]))
		return simp_ads, names_of_inputs, surface_reactions

	#Block of code between 143 and 315 is somewhat confusing and unorganized
	list_reacts = []
	dict_reacts = {}
	for i in range(0,len(name_in_smiles)):
		dict_reacts = make_molegraph(name_in_smiles[i], dict_reacts)

	for key in dict_reacts:
		list_reacts.append(dict_reacts[key])

	list_of_reactants = list_reacts

	rxn2 = recursive_scissions(list_of_reactants)

	simp_ads, names_of_inputs, surface_reactions = generate_gas_equations(list_of_reactants)

	diss_reacs =[]
	#print(names_of_inputs)
	for k in surface_reactions:
		temp_parts = k.replace('+','')
		neg,pos = k.split('<->')
		neg = neg.replace(' ','')
		if neg in names_of_inputs:
			new_exp = neg[:-1]+' + 2* <-> '+pos
		else:
			break
		diss_reacs.append(new_exp)

	merged_reactions = simp_ads+diss_reacs+surface_reactions

	#print(merged_reactions)

	return merged_reactions
#sys.exit()

#Ammonia Sythesis
#p_feed = ['N#N','[HH]']
#r_feed = ['N']

#Water-gas shift
#p_feed = ['[HH]','C(=O)=O']
#r_feed = ['[C-]#[O+]','O']

#list_of_reactants = p_feed + r_feed

#reactions = Langmuir_rates(list_of_reactants)

reactions = ['H2 + 2* <-> 2H*','CO + * <-> CO*', 'CO + O* <-> CO2','CO* + O* <-> CO2']

def deriv_and_const(reactions,monitor):

	rate_array, reactions, reactants = reac_list_parsing(reactions)

	#print(reactants)
	#print(reactions)
	#print()
	deriv_dict = {}
	rate_equations = []
	#sys.exit()


	#all_vals = 'NN'
	#all_vals_2 = 'NN'
	#all_vals_3 = all_vals_2
	
	for k,k2 in enumerate(reactants[0:monitor]):
		rate_equations.append('d('+str(k2)+')/dt = '+display_odes(k2,rate_array,reactions,reactants))
		#print('d('+str(k2)+')/dt = '+display_odes(k2,rate_array,reactions,reactants))
		for j,j2 in enumerate(reactants):
			#print(k2)
			#print(j2)
			new_expression = find_partial(k2,j2,rate_array,reactions,reactants)
			deriv_dict[(k,j)] = new_expression
			if new_expression != '0':
				pass
				#print('d(R_'+ k2 + ')/' + 'd('+j2+') = '+new_expression)
				
				#print(" ")

	return deriv_dict, reactants, rate_equations

#print('d('+str(all_vals)+')/dt = '+display_odes(all_vals))
#print("")
#print('d(R('+str(all_vals)+'))/(d('+all_vals_2+')) = '+find_partial(all_vals,all_vals_2))
#print("")
#print('d(R('+str(all_vals)+'))/(d( '+all_vals_2+' | '+all_vals_3+' )) = '+find_double_partial(all_vals,all_vals_2,all_vals_3))
#print("")

#for all_vals in reactants:
	#print('d('+str(all_vals)+')/dt = '+display_odes(all_vals,rate_array,reactions,reactants))
#	print('\\begin{equation}')
#	print('\\frac{\partial '+all_vals+'}{\partial t}} = '+display_odes(all_vals,rate_array,reactions,reactants))
#	print('\\end{equation}')
	#\frac{\partial C_{g}}{\partial x}\vert_{\hat\beta} = 0, \; x = 0
	#\end{equation}
	#print('d('+str(reactants[i])+')/dt = '+display_odes(reactants[i]))
#print('\\end{document}')



#generates each partial derivative for each ode
#print("")
#print('1st derivative')
#for all_vals in reactants:
#	for all_vals_2 in reactants:
#		#pass
#		new_expression = find_partial(all_vals,all_vals_2,rate_array,reactions,reactants)
#		if new_expression != '0':
#			print('d(R_'+ all_vals + ')/' + 'd('+all_vals_2+') = '+new_expression)
#			print(" ")

#steps through every combination of partial derivatives.
#can definetly reduce the necessary steps, but including all
#combinations (even redundant ones) doesn't take much time
#print('2nd derivative')
#for all_vals in reactants:
#	for blah,all_vals_2 in enumerate(reactants[:15]):
#		for all_vals_3 in reactants[blah:15]:

#			#print(all_vals)
#			why_not = find_double_partial(all_vals,all_vals_2,all_vals_3,rate_array,reactions,reactants)
#			if why_not != ' 0':
#				print('d(R_'+ all_vals + ')/' + 'd('+all_vals_2+' | '+all_vals_3+') =' +why_not)
#				print(" ")

#sys.exit()

#Ammonia Sythesis
#p_feed = ['N#N','[HH]']
#r_feed = ['N']

#list_of_reactants = p_feed + r_feed

#reactions = Langmuir_rates(list_of_reactants)

#rate_array, reactions, reactants = reac_list_parsing(reactions)

#print(display_odes('NN*',rate_array,reactions,reactants))



#Code below is just to store repeatedly used rate expressions.
#Checks if the same products/reactants have been used by creating
#a seperate file and storing or checking if it already exists
