import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random
import sys
import time
from reaction_list_odes import reac_list_parsing

def make_f_equation(reactions_n,reactants_number,Inert_only=False):

	F = ''
	rate_array, reactions, reactants, rev_irr = reac_list_parsing(reactions_n)
	#print(reactions)
	#print(reactants)
	#print(rate_array)

	gas_num = 0

	gas_num = len(reactants)
	molecules_in_gas_phase = 0

	for k in reactants:
		if '*' in k:
			break
		else:
			molecules_in_gas_phase += 1 
	
	gas_molecules = reactants[:gas_num]
	#print(reactants_number)
	#print(gas_num)
	#sys.exit()
	tail = len(reactants)
	#print("")
	for k in range(0,molecules_in_gas_phase):
		
		if k < (molecules_in_gas_phase-1):
			F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1) + " 
			#print("((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(1) + " )
			#print("")
		else:
			F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)" 
			#print("((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(1)" )
			#print("")
	for k in range(molecules_in_gas_phase,len(reactants)):
		F = F + " + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + dk*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)"
		#print(" + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + dk*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)")
		#print("")
	for k,z in enumerate(rate_array):
		
		neg = []
		val_neg = []
		pos = []
		val_pos = []
		for j,v in enumerate(z):
			if v < 0:
				neg.append(j)
				val_neg.append(v)
			elif v > 0:
				pos.append(j)
				val_pos.append(v)
			
		together = neg+pos

		new_neg = 'Ke'+str(k)
		for j,v in enumerate(neg):
			new_neg = new_neg+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"#
		new_pos = 'Kd'+str(k)
		for j,v in enumerate(pos):
			new_pos = new_pos+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"#
		
		for j,v in enumerate(together):
			if j < len(neg):
				#print("neg")
				if rev_irr[k] == 1:
					F = F+"- dk*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
					#print("- dk*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
				else:
					F = F+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
					#print("+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
			else:
				#print("pos")
				if rev_irr[k] == 1:
					F = F+"+ dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
					#print("+ dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
				else:
					F = F+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
					#print("- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")

	#Inert equation generator
	#print(len(reactants))
	#sys.exit()
	F = F + "+ ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)+ dk * Constant((np.sum(r_param))) * dot(w, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0) + ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']) )*v_d['v_"+str(len(reactants)+1)+"']*dx(1) + dk * Constant((np.sum(r_param))) * dot(w, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(1)"
	#print("+ Constant(eb[0]) * ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)+ dk * Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param))) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(0) + Constant(eb[1]) * ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']) )*v_d['v_"+str(len(reactants)+1)+"']*dx(1) + dk * Constant(D["+str(molecules_in_gas_phase)+"][1]/(np.sum(r_param))) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(1)")
	#print("")
	#F =  F+ "+ ((u_d['u_"+str(len(reactants))+"'] - u_nd['u_n"+str(len(reactants))+"']) )*v_d['v_"+str(len(reactants))+"']*dx(0)  + dk * Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param)))*dot(grad(u_d['u_"+str(len(reactants))+"']), grad(v_d['v_"+str(len(reactants))+"']))*dx(0) + ((u_d['u_"+str(len(reactants))+"'] - u_nd['u_n"+str(len(reactants))+"']) )*v_d['v_"+str(len(reactants))+"']*dx(1)  + dk*Constant(D["+str(molecules_in_gas_phase)+"][1]/(np.sum(r_param)))*dot(grad(u_d['u_"+str(len(reactants))+"']), grad(v_d['v_"+str(len(reactants))+"']))*dx(1)" 
	#sys.exit()
	element = '['
	for k in range(0,len(reactants)+1):
		if (k+1) < len(reactants)+1:
			element = element + 'P1,'
		else:
			element = element + 'P1]'
	necessary_dictionary = {'F': F,'gas_num': gas_num, 'surf_num': (len(reactants)-gas_num-1),'element': element,'gas_molecules': gas_molecules,'reactants_number':reactants_number,'reactants':reactants,'molecules_in_gas_phase':molecules_in_gas_phase}
	
	#print(F)
	#sys.exit()
	return necessary_dictionary

#test_new_again = make_f_equation(reactions_test)

#print(test_new_again['reactants'])