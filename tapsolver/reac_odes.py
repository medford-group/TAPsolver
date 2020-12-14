from numpy import inf
import numpy as np
from math import isnan
from math import sqrt
from math import *
import sympy
import time 
import sys
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
import time
import re

#from imolecule.notebook import generate
#from imolecule.format_converter import convert
#import networkx.algorithms.isomorphism as iso
#from networkx import is_isomorphic
#import networkx as nx
#from rxn_network import scissions, recursive_scissions, get_mechanisms
#from data_structures import MolGraph, RxnGraph

def variational_list_parsing(reactions):

	active_site_names = ['*','^','@','#']

	reacs = reactions.copy()
	reactants = []
	rev_irr = []

	#need to make this a function
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
	active_sites = 0
	for k,i in enumerate(reactants):
		if '*' not in i and '^' not in i and '@' not in i and '#' not in i:
			reactants.insert(insert_location, reactants.pop(k))
			insert_location += 1

		if i == '*':
			reactants.pop(k)
			if active_sites == 0:
				active_sites += 1
		if i == '^':
			reactants.pop(k)
			if active_sites == 1:
				active_sites += 1
		if i == '@':
			reactants.pop(k)
			if active_sites == 2:
				active_sites += 1
		if i == '#':
			reactants.pop(k)
			if active_sites == 3:
				active_sites += 1


	####Add the surface site to the end of the list of reactants
	###reactants.insert(len(reactants), '*')

	for kj in range(0,int(active_sites)):
		reactants.insert(len(reactants), active_site_names[kj])		

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
		rate_constant = "kf"+str(k)
		reaction_expression = reaction_expression + rate_constant
		for v,z in enumerate(reactants):
			if (test_value(-1,rate_array[k,v]) == True):
				#rate_constant = "k"+str(k+1)+"f"
				if abs(rate_array[k,v]) > 1:
					reaction_expression = reaction_expression + '*([' + reactants[v] + ']^'+str(int(abs(rate_array[k,v])))+')'
				else:
					reaction_expression = reaction_expression + '*[' + reactants[v] + ']'
		rate_constant = " - kb"+str(k)
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
				rate_constant = "kf"+str(k)
				pass_test = True
			elif rate_array[k,reactants.index(den)] > 0 and '<->' in reactions[k]:  #!!!!
				rate_constant = "kb"+str(k)
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
				rate_constant = "kf"+str(k)
			else:
				rate_constant = "kb"+str(k)
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

#Ammonia Sythesis
#p_feed = ['N#N','[HH]']
#r_feed = ['N']

#Water-gas shift
#p_feed = ['[HH]','C(=O)=O']
#r_feed = ['[C-]#[O+]','O']

#list_of_reactants = p_feed + r_feed

#reactions = Langmuir_rates(list_of_reactants)

#reactions = ['H2 + 2* <-> 2H*','CO + * <-> CO*', 'CO + O* <-> CO2','CO* + O* <-> CO2']

def deriv_and_const(reactions,monitor):

	rate_array, reactions, reactants = reac_list_parsing(reactions)

	#print(reactants)
	#print(reactions)
	#print()
	deriv_dict = {}
	rate_equations = []


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

def R_RRM_func(input_reactants,current_path,data_folder):

	path = current_path+'/'+data_folder+'/RRM_results/'
	try:  
		os.mkdir(path)
	except OSError:  
		print ("Creation of the directory %s failed" % path)
		pass
	else:  
		print ("Successfully created the directory %s " % path)
	
	#for j in r_const:
	#	r_const[j] = Constant(r_const[j])

	#user_data = pd.read_csv('./sim_file.csv',header=None)
	
	#reactor_info = user_data.iloc[2:18,:] 
	#feed_surf_info = user_data.iloc[21:27,:]
	#data_storage = user_data.iloc[29:35,:]
	#reaction_info = user_data.iloc[38:,:]
	

	user_data = pd.read_csv('./sim_file.csv',header=None)

	rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
	rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
	rows_3, cols_3 = np.where(user_data == 'Data_Storage_Options')
	rows_4, cols_4 = np.where(user_data == 'Reaction_Information')

	reactor_info = user_data.iloc[(2+rows_1[0]):(rows_2[0]-1),:] 
	feed_surf_info = user_data.iloc[2+rows_2[0]:rows_3[0]-1,:]
	data_storage = user_data.iloc[2+rows_3[0]:rows_4[0]-1,:]
	reaction_info = user_data.iloc[2+rows_4[0]:,:]

	reactor_kinetics_input = {}

	for k in range(0,len(feed_surf_info.index)):
		try:
			reactor_kinetics_input[feed_surf_info.iloc[k,5]] = float(feed_surf_info.iloc[k,1]) 
		except ValueError:
			reactor_kinetics_input[feed_surf_info.iloc[k,5]] = feed_surf_info.iloc[k,1]



	reactor_info = reactor_info.T
	a = reactor_info.iloc[5]
	reactor_info.iloc[0] = a
	reactor_info = reactor_info.iloc[:2,:]
	
	feed_surf_info = feed_surf_info.T
	a = feed_surf_info.iloc[5]
	feed_surf_info.iloc[0] = a
	feed_surf_info = feed_surf_info.iloc[:2,:]
	
	reactor_info = pd.concat([reactor_info,feed_surf_info], axis=1)
	
	data_storage = data_storage.T
	a = data_storage.iloc[5]
	data_storage.iloc[0] = a
	data_storage = data_storage.iloc[:2,:]
	
	reactor_info = pd.concat([reactor_info,data_storage], axis=1)
	
	reaction_info = reaction_info.T
	a = reaction_info.iloc[5]
	reaction_info.iloc[0] = a
	reaction_info = reaction_info.iloc[:2,:]
	
	reactor_info = pd.concat([reactor_info,reaction_info], axis=1)
	
	reactor_info.to_csv(data_folder+'/RRM_params.csv',header=None)
	
	#print(reactor_info)
	
	command = 'Rscript'
	
	#path2script = './simpleRRM.R'#'./yongeRRM.R'
	path2script = './yongeRRM.R'#'./yongeRRM.R
	
	feed = 1
	reactants = input_reactants
	in_reactants = reactants
	print(reactants)
	reactants = str(reactants)
	reactants = reactants.replace('[','')
	reactants = reactants.replace(']','')
	reactants = reactants.replace("!","")
	reactants = reactants.replace(' ','')
	
	arg1 = current_path+'/'+data_folder+'/'#"/home/adam/research_medford/python_code/tap_code/csv_input_sim/eley_eluc_folder/"
	arg2 = reactants
	arg3 = str(feed)
	
	pass_arg = []
	pass_arg.append(arg1)
	pass_arg.append(arg2)
	pass_arg.append(arg3)
	cmd = [command, path2script,arg1,arg2,arg3] # ... <- the arguments that you need to call to get the script to run properly
	
	subprocess.run(cmd) #Will just save the output data as a file
	

	textstr2 = 'RRM Model:'
	
	for k in range(len(in_reactants)):
		new_model = ''
		new_model += 'R_'+in_reactants[k]+' = '
		new_model += '\u03A8'+'_'+str(0)+' + '
		for j in range(len(in_reactants)):
			new_model += '\u03A8'+'_'+str(j+1)+'*C_'+''+in_reactants[j]+' + '
	
		for j in range(len(in_reactants)):
			new_model += '\u03A8'+'_'+str(len(in_reactants)+j+1)+'*U_'+''+in_reactants[j]+' + '
	
		for j in range(len(in_reactants)):
			for k in range(len(in_reactants)):
				new_model += '\u03A8'+'_'+str((j+2)*len(in_reactants)+k+1)+'*C_'+''+in_reactants[j]+'*U_'+''+in_reactants[k]+' + '
	
		for j in range(len(in_reactants)):
			for k in range(j,len(in_reactants)):
				#print((j-1)*len(in_reactants)+2*len(in_reactants)+len(in_reactants)*len(in_reactants)+k+1)
				new_model += '\u03A8'+'_'+str((j)*len(in_reactants)+2*len(in_reactants)+len(in_reactants)*len(in_reactants)+k+1-j)+'*U_'+''+in_reactants[j]+'*U_'+''+in_reactants[k]
				if j == (len(in_reactants)-1) and k == (len(in_reactants)-1):
					pass
				else:
					new_model += ' + '
		#print(new_model)
		#print()
		
		textstr2 = '\n'.join((textstr2,new_model))
	
	#	print()
	#	print()
	

	#total_u = 0
	#for j_k in range(0,len(in_reactants)):
	#	total_u += len(in_reactants)-j_k
	#total_u += len(in_reactants)*len(in_reactants)+len(in_reactants)*len(in_reactants)
	#print(total_u)
	#sys.exit()
	
	df_time = pd.read_csv(current_path+'/'+data_folder+'/flux_data/Inert.csv',header=None).iloc[:,0]

	rows = in_reactants.copy()
	for k in in_reactants.copy():
		rows.append('U_'+k)

	RRM_data = {}

	for k in in_reactants:
		RRM_data[k] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+k+'_reactivities.csv',header=None)

	y_proc_data = {}
	for k in range(int(reactor_kinetics_input['Number of Pulses'])):
		y_proc_data[k+1] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+str(k+1)+'_y_proc.csv',header=None).iloc[1:,:]
	#new_frame = y_proc_data[0+1].iloc[:,1+0+len(in_reactants)].copy()
	#print(new_frame)
	#sys.exit()
	
	#fig2, ax2 = plt.subplots()
	#for k in range(1,len(y_proc_data[1].columns)-3):
	#	#pd.to_numeric(new_frame,errors='coerce')
	#	ax2.plot(df_time,pd.to_numeric(y_proc_data[1].iloc[:,k]))
	#plt.show()

	#step through the pulses
	for z in range(int(reactor_kinetics_input['Number of Pulses'])):
		f, axarr = plt.subplots(len(in_reactants)*2, len(in_reactants))
		f = plt.gcf()
		f.set_size_inches(3*(len(in_reactants)), 1.5*len(rows))

		try:
			for ax, col in zip(axarr[0], in_reactants):
				#new_in = r'$R_%d$' % col
				ax.set_title(r'$\partial$'+'('+col+')'+' /')
				for ax, row in zip(axarr[:,0], rows):
					ax.set_ylabel('/ '+r'$\partial$'+'('+row+')', rotation=90, size='large')
		except TypeError:
			axarr[0].set_title(cols[0])
			for ax, row in zip(axarr[:], rows):
				ax.set_ylabel('/ '+r'$\partial$'+row, rotation=90, size='large')


		#step through the numerator
		for k in range(len(in_reactants)):
			#step through the gas denominator
			for j in range(len(in_reactants)):
				
				try:
					test = axarr[j,k]
				except IndexError:
					test = axarr[j]

				df_merged = pd.DataFrame(index=y_proc_data[1].index,columns=range(1))#.shape[0]

				df_merged.iloc[:] = float(0)#float(RRM_data[in_reactants[k]].iloc[z+1,j+2])
				
				#RRM_data[k] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+k+'_reactivities.csv',header=None)

				#print('dR_'+in_reactants[k]+' / dC_'+in_reactants[j],end=' = ')
				#print('\u03A8'+'_'+str((j+1)),end=' + ')

				for j_nu in range(0,len(in_reactants)):
					new_frame = y_proc_data[z+1].iloc[:,1+j_nu]#.copy()
					
					new_frame = pd.to_numeric(new_frame,errors='coerce')
					
					new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,2*len(in_reactants)+(k)*len(in_reactants)+2+j_nu])

					df_new = new_frame.to_frame()
					df_new.columns = [0]
					df_merged = df_merged.add(df_new,fill_value=0)
					#df_merged = df_merged + new_frame.to_frame()
					#df_merged[0] = df_merged.iloc[0] + new_frame.to_frame().iloc[0]
					#print(df_merged)
					#sys.exit()
					#df_merged.iloc[0] = df_merged.iloc[0] + new_frame#+ RRM_data[in_reactants[k]].iloc[z+1,2*len(in_reactants)+(k)*len(in_reactants)+2+j_nu]*y_proc_data[z+1].iloc[:,j_nu]

				#sys.exit()
				#print(df_merged)
				#sys.exit()

				#for j_nu in range(0,len(in_reactants)):
				#	print('\u03A8'+'_'+str(2*len(in_reactants)+(k)*len(in_reactants)+1+j_nu)+' * U_'+in_reactants[j_nu],end=' + ')
				test.plot(df_time, df_merged[0])

			#step through the surface denominator
			for j in range(len(in_reactants)):
				try:
					test = axarr[j+len(in_reactants),k]
				except IndexError:
					test = axarr[j+len(in_reactants)]

				df_merged = pd.DataFrame(index=y_proc_data[1].index,columns=range(1))#.shape[0]

				df_merged.iloc[:] = 0#float(RRM_data[in_reactants[k]].iloc[z+1,len(in_reactants)+j+2])



				#aditional terms
				#for step_1 in range(0,len(in_reactants)):

				print()
				# Graphs for dU portion
				#for j_nu in range(0,len(in_reactants)):
				curr_value = k + 2*len(in_reactants) 
				print(curr_value)
				for step_2 in range(0,k):
					new_frame = y_proc_data[z+1].iloc[:,1+step_2+len(in_reactants)].copy()
					#print(new_frame)
					
					new_frame = pd.to_numeric(new_frame,errors='coerce')

					new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1])
					#print(float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1]))

					df_new = new_frame.to_frame()
					df_new.columns = [0]
					df_merged = df_merged.add(df_new,fill_value=0)

					curr_value += (len(in_reactants)-step_2-1)

					df_merged.iloc[0] = df_merged.iloc[0] + new_frame
					print(curr_value)
				#if k > 1:
				#	sys.exit()
				curr_value = k + 2*len(in_reactants) + len(in_reactants)*len(in_reactants) + 1 #step_1
				for step_3 in range(k,len(in_reactants)):#curr_value+
					new_frame = y_proc_data[z+1].iloc[:,1+step_3].copy()
					new_frame = pd.to_numeric(new_frame,errors='coerce')
					#new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1])
					curr_value += 1
					df_new = new_frame.to_frame()
					df_new.columns = [0]
					df_merged = df_merged.add(df_new,fill_value=0)
					#print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_3],end=' + ')
					df_merged.iloc[0] = df_merged.iloc[0] + new_frame
				#sys.exit()
				#print(df_time)
				#print(df_merged[0])
				#print(type(df_time))
				#print(type(df_merged[0]))
				#sys.exit()
				test.plot(df_time, df_merged)
				
				
				####Generate the derivative equations
				#print('dR_'+in_reactants[k]+'/dU_'+in_reactants[j],end=' = ')
				#print('\u03A8'+'_'+str(len(in_reactants)+(j+1)),end=' + ')
				#curr_value = k + 2*len(in_reactants) + len(in_reactants)*len(in_reactants) + 1 #step_1
				#for j_nu in range(0,len(in_reactants)):
				#	print('\u03A8'+'_'+str(2*len(in_reactants)+(j_nu)*len(in_reactants)+1+j)+' * C_'+in_reactants[j_nu],end=' + ')


				#for step_2 in range(0,k):
				#	print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_2],end=' + ')
					
			
				#for step_3 in range(k,len(in_reactants)):#curr_value+
				#	print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_3],end=' + ')
				#	curr_value += 1
				#print()
				#print()


		props = dict(boxstyle='round', facecolor='white', alpha=0.4)
		plt.gcf().text(0.5, 0.05, textstr2, fontsize=8, bbox=props,ha='center')
		f.tight_layout()	
		f.subplots_adjust(bottom=.2,left=0.3,top=0.9,right=0.7)
		f.suptitle("RRM Jacobian", fontsize=14)
		#f.text(0.15, 0.6, 'Derivative Denominator', va='center', rotation='vertical', fontsize=14)
		
		plt.savefig(current_path+'/'+data_folder+'/graphs/RRM_p'+str(z+1)+'.png')#current_path+'/'+data_folder+
		plt.show()

def petal_plots_RRM(folder_location):

	print('good pass')

def MKM_graphs(kin_params,reactions_in,fold_loc,graph_display):
	#print(kin_params)
	kf0 = 1
	gasses = 2
	
	folder_location = './'+fold_loc+'/'#./new_test_folder/'
	
	reactions = reactions_in#['CO + * <-> CO*','CO + O* <-> CO2','CO* + O* <-> CO2']
	reactions_input = reactions.copy()
	deriv_dict, reactants, rate_equations = deriv_and_const(reactions,gasses)
	
	
	deriv = len(reactants)

	textstr = 'Elementary Reactions:'
	for k_elem in reactions_input:
		textstr = '\n'.join((textstr,k_elem))
	
	textstr2 = 'Rate Expressions:'
	for k_reac in rate_equations:
		textstr2 = '\n'.join((textstr2,k_reac))
	
	cols = reactants[:gasses]
	
	rows = reactants#['CO','CO2','CO*','O*','*']
	
	df_time = pd.read_csv(folder_location+'*'+'.csv',header=None).iloc[:,0]
	
	#example_array = [['-1*(ke0*[*]*[CO*]) + -1*(ke1*[O*])','1*(ke1*[O*])'],['-1*(kd1)','0'],['0','1*(ke2*[O*])'],['-1*(ke1*[CO])','1*(ke1*[CO]) + 1*(ke2*[CO*])'],['-1*(ke0*[CO])','0']]
	
	f, axarr = plt.subplots(deriv, gasses)
	
	f = plt.gcf()
	f.set_size_inches(3*len(cols), 1.5*len(rows))
	#fig.savefig('test2png.png', dpi=100)
	
	
	try:
		for ax, col in zip(axarr[0], cols):
		    ax.set_title(col)
		    for ax, row in zip(axarr[:,0], rows):
	    		ax.set_ylabel(row, rotation=90, size='large')
	except TypeError:
		axarr[0].set_title(cols[0])
		for ax, row in zip(axarr[:], rows):
			ax.set_ylabel(row, rotation=90, size='large')
	
	
	for k_der in range(0,deriv):
		for j_spe in range(0,gasses):
			try:
				test = axarr[k_der,j_spe]
			except IndexError:
				test = axarr[k_der]
			test.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			try:
				new = deriv_dict[(j_spe,k_der)].split('+')#example_array[k_der][j_spe].split('+')
			except NameError:
				new = deriv_dict[(j_spe,k_der)]#example_array[k_der][j_spe]
			df_merged = pd.DataFrame(index=df_time.index,columns=range(1))
			df_merged.iloc[:] = 0
			for z in new:
				if '(' in z:
					result = (z.split('('))[1].split(')')[0]
					#result = re.search('%s(.*)%s' % ('(', ')'), z)
					result = result.split('*[')
					for jk_num,jk in enumerate(result.copy()):
						if "]" in jk:
							result[jk_num] = jk.replace("]",'')
					#print(result)
					if len(result) > 1:
						new_df = pd.read_csv(folder_location+'/'+result[1]+'.csv',header=None).iloc[:,1]
						if len(result) > 2:
							for k_cool in range(2,len(result)):
								#print(result[k_cool])
								df_new = pd.read_csv(folder_location+'/'+result[k_cool]+'.csv',header=None).iloc[:,1]
								new_df = new_df*df_new
						if '-' in z:
							new_df *= -kin_params[result[0].replace('k', 'K')]#ke0#str
						else:
							new_df *= kin_params[result[0].replace('k', 'K')]#kin_params[result[0]]#ke0
					elif (z.split('('))[1].split(')')[0] != '':
						new_df = pd.DataFrame(index=df_time.index,columns=range(1))
						if '-' in z:
							new_df.iloc[:] = -kin_params[result[0].replace('k', 'K')]#kin_params[result[0]]
						else:
							new_df.iloc[:] = kin_params[result[0].replace('k', 'K')]#kin_params[result[0]]#ke0
					else:
						new_df.iloc[:] = 0
				temp = df_merged 
				df_merged[0] = temp[0] + new_df
				
				#print(df_merged)
					#test.plot(df_time, df_merged)
			#print(df_merged)
			#sys.exit()
			props = dict(boxstyle='round', facecolor='white', alpha=0.0)
			#axarr[k_der, j_spe].set_xlabel(example_array[k_der][j_spe])
			if deriv_dict[(j_spe,k_der)] != '0':#example_array[k_der][j_spe] != '0':
				result = (deriv_dict[(j_spe,k_der)].split('('))[1].split(')')[0]
				if result != '':
					pass
					#if '+' in deriv_dict[(j_spe,k_der)]:#example_array[k_der][j_spe]:
					#	new_label = deriv_dict[(j_spe,k_der)].split('+')#example_array[k_der][j_spe].split('+')
					#	#test.text(0.5, 0.77, new_label[0],str('+')+new_label[1], transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')
					#	test.text(0.5, -0.3, '\n'.join((new_label[0],str('+')+new_label[1])), transform=test.transAxes, fontsize=7,verticalalignment='bottom', bbox=props,ha='center')
					#else:	
					#	test.text(0.5, -0.3, deriv_dict[(j_spe,k_der)], transform=test.transAxes, fontsize=7,verticalalignment='bottom', bbox=props,ha='center')#, bbox=props
					#	#test.text(0.5, 0.77, '\n'.join(('Deriv:',deriv_dict[(j_spe,k_der)])), transform=test.transAxes, fontsize=8,verticalalignment='top', bbox=props,ha='center')
			else:
				df_merged = pd.DataFrame(index=df_time.index,columns=range(1))
				df_merged.iloc[:] = 0
			
			test.plot(df_time, df_merged)

	props = dict(boxstyle='round', facecolor='white', alpha=0.4)
	plt.gcf().text(0.75, 0.85, textstr, fontsize=8, bbox=props,ha='left')
	plt.gcf().text(0.5, 0.05, textstr2, fontsize=8, bbox=props,ha='center')
	f.tight_layout()	
	f.subplots_adjust(bottom=.2,left=0.3,top=0.9,right=0.7)
	f.suptitle("MKM Jacobian", fontsize=14)
	plt.savefig(fold_loc+'/../graphs/Jacob_MKM.png')
	
	if graph_display == 'true':
		plt.show()
	else:
		plt.clf()
		plt.close()
	
				#try:
				#	result = re.search('[(.*)]',z)
				#	print(result.group(0))
				#	print('Test')
				#except TypeError:
				#	print('Failed')
				#	pass
	
				#if '-' in new:
				#	print('negative')
				#else:
				#	print('positive')

def petal_plots_exp(folder_location):

	print('Petal Plot Pass Confirmed')



def jacobian_visual(kin_params,reactions_in,fold_loc,graph_display,input_reactants,current_path,data_folder,in_reactants,number_of_pulses):
	
	df_time = pd.read_csv(current_path+'/'+data_folder+'/flux_data/Inert.csv',header=None).iloc[:,0]

	rows = in_reactants.copy()
	for k in in_reactants.copy():
		rows.append('U_'+k)

	RRM_data = {}

	for k in in_reactants:
		RRM_data[k] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+k+'_reactivities.csv',header=None)

	y_proc_data = {}
	for k in range(int(number_of_pulses)):
		y_proc_data[k+1] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+str(k+1)+'_y_proc.csv',header=None).iloc[1:,:]
	


	kf0 = 1
	gasses = 2
	
	folder_location = './'+fold_loc+'/'#./new_test_folder/'
	
	reactions = reactions_in#['CO + * <-> CO*','CO + O* <-> CO2','CO* + O* <-> CO2']
	reactions_input = reactions.copy()
	deriv_dict, reactants, rate_equations = deriv_and_const(reactions,gasses)
	
	
	deriv = len(reactants)

	textstr = 'Elementary Reactions:'
	for k_elem in reactions_input:
		textstr = '\n'.join((textstr,k_elem))
	
	textstr2 = 'Rate Expressions:'
	for k_reac in rate_equations:
		textstr2 = '\n'.join((textstr2,k_reac))
	
	cols = reactants[:gasses]
	
	rows = reactants#['CO','CO2','CO*','O*','*']
	
	df_time = pd.read_csv(folder_location+'*'+'.csv',header=None).iloc[:,0]
	
	#example_array = [['-1*(ke0*[*]*[CO*]) + -1*(ke1*[O*])','1*(ke1*[O*])'],['-1*(kd1)','0'],['0','1*(ke2*[O*])'],['-1*(ke1*[CO])','1*(ke1*[CO]) + 1*(ke2*[CO*])'],['-1*(ke0*[CO])','0']]
	
	f, axarr = plt.subplots(deriv, gasses)
	
	f = plt.gcf()
	f.set_size_inches(3*len(cols), 1.5*len(rows))
	#fig.savefig('test2png.png', dpi=100)
	
	
	try:
		for ax, col in zip(axarr[0], cols):
		    ax.set_title(col)
		    for ax, row in zip(axarr[:,0], rows):
	    		ax.set_ylabel(row, rotation=90, size='large')
	except TypeError:
		axarr[0].set_title(cols[0])
		for ax, row in zip(axarr[:], rows):
			ax.set_ylabel(row, rotation=90, size='large')
	
	
	for k_der in range(0,deriv):
		for j_spe in range(0,gasses):
			try:
				test = axarr[k_der,j_spe]
			except IndexError:
				test = axarr[k_der]
			test.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
			try:
				new = deriv_dict[(j_spe,k_der)].split('+')#example_array[k_der][j_spe].split('+')
			except NameError:
				new = deriv_dict[(j_spe,k_der)]#example_array[k_der][j_spe]
			df_merged = pd.DataFrame(index=df_time.index,columns=range(1))
			df_merged.iloc[:] = 0
			for z in new:
				if '(' in z:
					result = (z.split('('))[1].split(')')[0]
					#result = re.search('%s(.*)%s' % ('(', ')'), z)
					result = result.split('*[')
					for jk_num,jk in enumerate(result.copy()):
						if "]" in jk:
							result[jk_num] = jk.replace("]",'')
					#print(result)
					if len(result) > 1:
						new_df = pd.read_csv(folder_location+'/'+result[1]+'.csv',header=None).iloc[:,1]
						if len(result) > 2:
							for k_cool in range(2,len(result)):
								#print(result[k_cool])
								df_new = pd.read_csv(folder_location+'/'+result[k_cool]+'.csv',header=None).iloc[:,1]
								new_df = new_df*df_new
						if '-' in z:
							new_df *= -kin_params[result[0].replace('k', 'k')]#ke0#str
						else:
							new_df *= kin_params[result[0].replace('k', 'k')]#kin_params[result[0]]#ke0
					elif (z.split('('))[1].split(')')[0] != '':
						new_df = pd.DataFrame(index=df_time.index,columns=range(1))
						if '-' in z:
							new_df.iloc[:] = -kin_params[result[0].replace('k', 'k')]#kin_params[result[0]]
						else:
							new_df.iloc[:] = kin_params[result[0].replace('k', 'k')]#kin_params[result[0]]#ke0
					else:
						new_df.iloc[:] = 0
				temp = df_merged 
				df_merged[0] = temp[0] + new_df
				
				#print(df_merged)
					#test.plot(df_time, df_merged)
			#print(df_merged)
			#sys.exit()
			props = dict(boxstyle='round', facecolor='white', alpha=0.0)
			#axarr[k_der, j_spe].set_xlabel(example_array[k_der][j_spe])
			if deriv_dict[(j_spe,k_der)] != '0':#example_array[k_der][j_spe] != '0':
				result = (deriv_dict[(j_spe,k_der)].split('('))[1].split(')')[0]
				if result != '':
					pass
					#if '+' in deriv_dict[(j_spe,k_der)]:#example_array[k_der][j_spe]:
					#	new_label = deriv_dict[(j_spe,k_der)].split('+')#example_array[k_der][j_spe].split('+')
					#	#test.text(0.5, 0.77, new_label[0],str('+')+new_label[1], transform=test.transAxes, fontsize=8,verticalalignment='top',bbox=props,ha='center')
					#	test.text(0.5, -0.3, '\n'.join((new_label[0],str('+')+new_label[1])), transform=test.transAxes, fontsize=7,verticalalignment='bottom', bbox=props,ha='center')
					#else:	
					#	test.text(0.5, -0.3, deriv_dict[(j_spe,k_der)], transform=test.transAxes, fontsize=7,verticalalignment='bottom', bbox=props,ha='center')#, bbox=props
					#	#test.text(0.5, 0.77, '\n'.join(('Deriv:',deriv_dict[(j_spe,k_der)])), transform=test.transAxes, fontsize=8,verticalalignment='top', bbox=props,ha='center')
			else:
				df_merged = pd.DataFrame(index=df_time.index,columns=range(1))
				df_merged.iloc[:] = 0


			test.plot(df_time, df_merged)

			if k_der < 4:
				if k_der < 2:

					df_merged = pd.DataFrame(index=y_proc_data[1].index,columns=range(1))#.shape[0]

					df_merged.iloc[:] = float(0)#float(RRM_data[in_reactants[k]].iloc[z+1,j+2])
				
					#RRM_data[k] = pd.read_csv(current_path+'/'+data_folder+'/RRM_results/'+k+'_reactivities.csv',header=None)

					#print('dR_'+in_reactants[k]+' / dC_'+in_reactants[j],end=' = ')
					#print('\u03A8'+'_'+str((j+1)),end=' + ')

					for j_nu in range(0,len(in_reactants)):
						new_frame = y_proc_data[1].iloc[:,1+j_nu]#.copy()
						
						new_frame = pd.to_numeric(new_frame,errors='coerce')
						
						new_frame *= float(RRM_data[in_reactants[k_der]].iloc[1,2*len(in_reactants)+(k_der)*len(in_reactants)+2+j_nu])

						df_new = new_frame.to_frame()
						df_new.columns = [0]
						df_merged = df_merged.add(df_new,fill_value=0)
						#df_merged = df_merged + new_frame.to_frame()
						#df_merged[0] = df_merged.iloc[0] + new_frame.to_frame().iloc[0]
						#print(df_merged)
						#sys.exit()
						#df_merged.iloc[0] = df_merged.iloc[0] + new_frame#+ RRM_data[in_reactants[k]].iloc[z+1,2*len(in_reactants)+(k)*len(in_reactants)+2+j_nu]*y_proc_data[z+1].iloc[:,j_nu]

					#sys.exit()
					#print(df_merged)
					#sys.exit()

					#for j_nu in range(0,len(in_reactants)):
					#	print('\u03A8'+'_'+str(2*len(in_reactants)+(k)*len(in_reactants)+1+j_nu)+' * U_'+in_reactants[j_nu],end=' + ')
					test.plot(df_time, df_merged[0])
				elif k_der == 100:
			#step through the surface denominator
					for j in range(len(in_reactants)):

						df_merged = pd.DataFrame(index=y_proc_data[1].index,columns=range(1))#.shape[0]

						df_merged.iloc[:] = 0#float(RRM_data[in_reactants[k]].iloc[z+1,len(in_reactants)+j+2])



						#aditional terms
						#for step_1 in range(0,len(in_reactants)):

						print()
						# Graphs for dU portion
						#for j_nu in range(0,len(in_reactants)):
						curr_value = k_der + 2*len(in_reactants)
						for step_2 in range(0,k_der):
							new_frame = y_proc_data[1].iloc[:,1+step_2+len(in_reactants)].copy()
							#print(new_frame)
					
							new_frame = pd.to_numeric(new_frame,errors='coerce')
							#print(k_der)
							#sys.exit()
							new_frame *= float(RRM_data[in_reactants[k_der]].iloc[1,curr_value+1])
							#print(float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1]))

							df_new = new_frame.to_frame()
							df_new.columns = [0]
							df_merged = df_merged.add(df_new,fill_value=0)

							curr_value += (len(in_reactants)-step_2-1)

							df_merged.iloc[0] = df_merged.iloc[0] + new_frame
						#if k > 1:
						#	sys.exit()
						curr_value = k_der + 2*len(in_reactants) + len(in_reactants)*len(in_reactants) + 1 #step_1
						
						for step_3 in range(k,len(in_reactants)):#curr_value+
							new_frame = y_proc_data[1].iloc[:,1+step_3].copy()
							new_frame = pd.to_numeric(new_frame,errors='coerce')
							#new_frame *= float(RRM_data[in_reactants[k]].iloc[z+1,curr_value+1])
							curr_value += 1
							df_new = new_frame.to_frame()
							df_new.columns = [0]
							df_merged = df_merged.add(df_new,fill_value=0)
							#print('\u03A8'+'_'+str(curr_value)+' * U_'+in_reactants[step_3],end=' + ')
							df_merged.iloc[0] = df_merged.iloc[0] + new_frame
				#sys.exit()
						#print(df_time)
						#print(df_merged[0])
						#print(type(df_time))
						#print(type(df_merged[0]))
						#sys.exit()

					test.plot(df_time, df_merged[0])






	props = dict(boxstyle='round', facecolor='white', alpha=0.4)
	plt.gcf().text(0.75, 0.85, textstr, fontsize=8, bbox=props,ha='left')
	plt.gcf().text(0.5, 0.05, textstr2, fontsize=8, bbox=props,ha='center')
	f.tight_layout()	
	f.subplots_adjust(bottom=.2,left=0.3,top=0.9,right=0.7)
	f.suptitle("MKM Jacobian", fontsize=14)
	plt.show()
