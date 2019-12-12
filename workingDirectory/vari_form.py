import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
import random
import sys
import time
from reac_odes import deriv_and_const,variational_list_parsing

def make_f_equation(reactions_n,reactants_number,reactor_type,active_sites,number_of_inerts,advection,Inert_only=False):

	active_site_names = ['*','^','@','#']

	if Inert_only == True:
		F = "Constant(eb[0]) * ((u - u_n))*v*dx(0) + dk * ( Dout[0]/Constant(np.sum(r_param)**2) ) * dot(grad(u), grad(v))*dx(0) + Constant(eb[1]) * ((u - u_n) )*v*dx(1) + dk * ( Din[0]/Constant(np.sum(r_param)**2) ) * dot(grad(u), grad(v))*dx(1)"
		necessary_dictionary = {'F': F,'gas_num': 1, 'surf_num': 0}
		return necessary_dictionary

	F = ''
	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions_n,active_sites)
	
	#sys.exit()
	#print(rate_array)
	#sys.exit()
	#print(reactants)
	#Removed the removal of the *
		
	#####If I need to include the number of active sites or not, alter the following code####
	#reactants = reactants[:-1]
	#print(type(rate_array))
	#temp_arr = rate_array.copy()
	#rate_array = np.delete(rate_array,-1,axis=1)

	#new_rate_array = []
	#for k,i in enumerate(temp_arr):
		#print(i)
	#	new_rate_array.append(rate_array[k][:-1])

	#########################################################################################

	gas_num = 0

	gas_num = len(reactants)
	molecules_in_gas_phase = 0

	for k in reactants:
		if '*' in k:
			break
		else:
			molecules_in_gas_phase += 1 
	
	gas_molecules = reactants[:gas_num]
	
	gas_num = len(reactants)+int(number_of_inerts)

	#print(reactants_number)
	#print(gas_num)
	#sys.exit()
	tail = len(reactants)

	#rate_array, reactions, reactants, rev_irr

	#print("")

	#time.sleep(4)
	
	
	if reactor_type == 'tap':
		for k in range(0,molecules_in_gas_phase):
			if k < (molecules_in_gas_phase-1):
				F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) ) *dot(grad(theta*u_d['u_"+str(k+1)+"']+(1-theta)*u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)     + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)       + dk*(Din["+str(k)+"]/(eb[1]*np.sum(r_param)**2) )*dot(grad(theta*u_d['u_"+str(k+1)+"'] + (1-theta)*u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1) + " 
				#F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) ) *dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)       + dk*Constant(1-theta)*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)      + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)       + dk*Constant(theta)*(Din["+str(k)+"]/(eb[1]*np.sum(r_param)**2) )*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)      + dk*Constant(1-theta)*(Din["+str(k)+"]/(eb[1]*np.sum(r_param)**2) )*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1) + " 
				# current # F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) ) *dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)       + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)       + dk*Constant(theta)*(Din["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)     + " 
				
				if advection.lower() == 'true':
					F = F+ "dk*Constant(theta)*Constant((np.sum(r_param)))*dot(advTerm, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(advTerm, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)           + "
				#F = F+"Constant(eb[0])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)       + dk*Constant(1-theta)*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)      + Constant(eb[1])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)       + dk*Constant(theta)*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)      + dk*Constant(1-theta)*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1) + " 
				#print("Constant(eb[0])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0) + Constant(eb[1])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1) + " )
				#print("")
			else:
				F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)         + dk*Constant(1-theta)*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)     + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)      + dk*Constant(theta)*(Din["+str(k)+"]/(eb[1]*np.sum(r_param)**2) )*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)      + dk*Constant(1-theta)*(Din["+str(k)+"]/(eb[1]*np.sum(r_param)**2) )*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)" 
				# current # F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)         + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)      + dk*Constant(theta)*(Din["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)  " 
				
				if advection.lower() == 'true':
					F = F+ " + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(advTerm, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)    + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(advTerm, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)  "

				#F = F+"Constant(eb[0])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + dk*Constant(theta)*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)         + dk*Constant(1-theta)*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)     + Constant(eb[1])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)      + dk*Constant(theta)*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)      + dk*Constant(1-theta)*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_nd['u_n"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)" 
				#print("Constant(eb[0])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0) + Constant(eb[1])*((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)" )
				#print("")
	elif reactor_type == 't_pfr':
		for k in range(0,molecules_in_gas_phase):
			if k < (molecules_in_gas_phase-1):
				F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)       + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)          + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)       + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)  + " 
				#print("((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(1) + " )
				#print("")
			else:
				F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)       + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)      + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)           + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)" 
				#print("((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(1)" )
				#print("")
	elif reactor_type == 't_pfr_diff':
		for k in range(0,molecules_in_gas_phase):
			if k < (molecules_in_gas_phase-1):
				F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)        + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)         + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)               + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)        + dk*Constant(theta)*Constant(Dout["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)          + dk*Constant(theta)*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1)  +" 
				#print("((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(1) + " )
				#print("")
			else:
				F = F+"((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)        + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)        + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)             + dk*Constant(theta)*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(1)     + dk*Constant(theta)*Constant(D["+str(k)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(0)       + dk*Constant(theta)*Constant(D["+str(k)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(k+1)+"']), grad(v_d['v_"+str(k+1)+"']))*dx(1) " 
			#print("((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)  + dk*Constant((np.sum(r_param)))*dot(w, grad(u_d['u_"+str(k+1)+"']))*dx(1)" )
			#print("")


	if reactor_type != 'batch':
		for k in range(molecules_in_gas_phase,len(reactants)):
			F = F + " + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)"
			#print(" + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']))*v_d['v_"+str(k+1)+"']*dx(0) + ((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx(1)")
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
	
			new_neg = 'r_const["kf'+str(k)+'"]'
			for j,v in enumerate(neg):
				new_neg = new_neg+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"#
			new_pos = 'r_const["kb'+str(k)+'"]'
			for j,v in enumerate(pos):
				new_pos = new_pos+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"#

			for j,v in enumerate(together):
				if j < len(neg):
					#print("neg")
					if rev_irr[k] == 1:
						F = F+"- dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("- dk*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
					else:
						F = F+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
				else:
					#print("pos")
					if rev_irr[k] == 1:
						F = F+"+ dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("+ dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
					else:
						F = F+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")

#######################################################

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
	
			new_neg = 'r_const["kf'+str(k)+'"]'
			for j,v in enumerate(neg):
				new_neg = new_neg+"*(u_nd['u_n"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"#
			new_pos = 'r_const["kb'+str(k)+'"]'
			for j,v in enumerate(pos):
				new_pos = new_pos+"*(u_nd['u_n"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"#

			for j,v in enumerate(together):
				if j < len(neg):
					#print("neg")
					if rev_irr[k] == 1:
						F = F+"- dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("- dk*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
					else:
						F = F+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
				else:
					#print("pos")
					if rev_irr[k] == 1:
						F = F+"+ dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("+ dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
					else:
						F = F+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))"
						#print("- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
						#print("")
	
#######################################################


		if reactor_type == 'tap':
			for k in range(0,int(number_of_inerts)):
				F = F + "+  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']))*v_d['v_"+str(len(reactants)+1+k)+"']*dx(0)        + dk *( Dout["+str(molecules_in_gas_phase+k)+"]/Constant(eb[0]*np.sum(r_param)**2) ) * dot(grad( theta*u_d['u_"+str(len(reactants)+1+k)+"'] + (1-theta)*u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)     +  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']) )*v_d['v_"+str(len(reactants)+1+k)+"']*dx(1)             + dk * (Din["+str(molecules_in_gas_phase+k)+"]  /Constant(eb[1]*np.sum(r_param)**2) ) * dot(grad( theta*u_d['u_"+str(len(reactants)+1+k)+"'] + (1-theta)*u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)"
				
				#F = F + "+  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']))*v_d['v_"+str(len(reactants)+1+k)+"']*dx(0)        + dk *Constant(theta)*( Dout["+str(molecules_in_gas_phase+k)+"]/Constant(eb[0]*np.sum(r_param)**2) ) * dot(grad(u_d['u_"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)          + dk *Constant(1-theta)*(Dout["+str(molecules_in_gas_phase+k)+"] /Constant(eb[0]*np.sum(r_param)**2) )  * dot(grad(u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)          +  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']) )*v_d['v_"+str(len(reactants)+1+k)+"']*dx(1)             + dk *Constant(theta)* (Din["+str(molecules_in_gas_phase+k)+"]  /Constant(eb[1]*np.sum(r_param)**2) ) * dot(grad(u_d['u_"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)          + dk *Constant(1-theta)*(Din["+str(molecules_in_gas_phase+k)+"]  /Constant(eb[1]*np.sum(r_param)**2) ) * dot(grad(u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)"
				# current # F = F + "+  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']))*v_d['v_"+str(len(reactants)+1+k)+"']*dx(0)        + Constant(1) *dk *Constant(theta)*( D["+str(molecules_in_gas_phase+k)+"][0] /Constant(eb[0]*np.sum(r_param)**2) ) * dot(grad(u_d['u_"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)       +  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']) )*v_d['v_"+str(len(reactants)+1+k)+"']*dx(1)             + Constant(1) *dk *Constant(theta)* (Din["+str(molecules_in_gas_phase+k)+"]  /Constant(eb[0]*np.sum(r_param)**2) ) * dot(grad(u_d['u_"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)      "
				
				#print("+  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']))*v_d['v_"+str(len(reactants)+1+k)+"']*dx(0)        + Constant(2*eb[0]) *dk *Constant(theta)*( D["+str(molecules_in_gas_phase+k)+"][0] /Constant(np.sum(r_param)**2) ) * dot(grad(u_d['u_"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)          + Constant(2*eb[0]) * dk *Constant(1-theta)*(D["+str(molecules_in_gas_phase+k)+"][0] /Constant(np.sum(r_param)**2) )  * dot(grad(u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(0)          +  ((u_d['u_"+str(len(reactants)+1+k)+"'] - u_nd['u_n"+str(len(reactants)+1+k)+"']) )*v_d['v_"+str(len(reactants)+1+k)+"']*dx(1)             + Constant(2*eb[1]) *dk *Constant(theta)* (D["+str(molecules_in_gas_phase+k)+"][1]  /Constant(np.sum(r_param)**2) ) * dot(grad(u_d['u_"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)          + Constant(2*eb[1]) *dk *Constant(1-theta)*(D["+str(molecules_in_gas_phase+k)+"][1]  /Constant(np.sum(r_param)**2) ) * dot(grad(u_nd['u_n"+str(len(reactants)+1+k)+"']), grad(v_d['v_"+str(len(reactants)+1+k)+"']))*dx(1)")
				if advection.lower() == 'true':
					F = F+ "+ dk *Constant(theta)* Constant((np.sum(r_param))) * dot(advTerm, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)        + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(advTerm, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(1) "


			#F = F + "+ Constant(eb[0]) * ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)        + dk *Constant(theta)* Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(0)          + dk *Constant(1-theta)* Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param)**2)) * dot(grad(u_nd['u_n"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(0)          + Constant(eb[1]) * ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']) )*v_d['v_"+str(len(reactants)+1)+"']*dx(1)             + dk *Constant(theta)* Constant(D["+str(molecules_in_gas_phase)+"][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(1)          + dk *Constant(1-theta)* Constant(D["+str(molecules_in_gas_phase)+"][1]/(np.sum(r_param)**2)) * dot(grad(u_nd['u_n"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(1)"
		elif reactor_type == 't_pfr':
			F = F + "+ ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)          + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(w, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)           + ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']) )*v_d['v_"+str(len(reactants)+1)+"']*dx(1)          + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(w, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(1)   "
		elif reactor_type == 't_pfr_diff':
			F = F + "+ ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)          + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(w, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)        + ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']) )*v_d['v_"+str(len(reactants)+1)+"']*dx(1)           + dk *Constant(theta)* Constant((np.sum(r_param))) * dot(w, grad(u_d['u_"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(1)   + dk *Constant(theta)* Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(0) "		
		#print("+ Constant(eb[0]) * ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']))*v_d['v_"+str(len(reactants)+1)+"']*dx(0)+ dk * Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param)**2)) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(0) + Constant(eb[1]) * ((u_d['u_"+str(len(reactants)+1)+"'] - u_nd['u_n"+str(len(reactants)+1)+"']) )*v_d['v_"+str(len(reactants)+1)+"']*dx(1) + dk * Constant(D["+str(molecules_in_gas_phase)+"][1]/(np.sum(r_param)**2)) * dot(grad(u_d['u_"+str(len(reactants)+1)+"']), grad(v_d['v_"+str(len(reactants)+1)+"']))*dx(1)")
		#print("")
		#F =  F+ "+ Constant(eb[0])*((u_d['u_"+str(len(reactants))+"'] - u_nd['u_n"+str(len(reactants))+"']) )*v_d['v_"+str(len(reactants))+"']*dx(0)  + dk * Constant(D["+str(molecules_in_gas_phase)+"][0]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(len(reactants))+"']), grad(v_d['v_"+str(len(reactants))+"']))*dx(0) + Constant(eb[1])*((u_d['u_"+str(len(reactants))+"'] - u_nd['u_n"+str(len(reactants))+"']) )*v_d['v_"+str(len(reactants))+"']*dx(1)  + dk*Constant(D["+str(molecules_in_gas_phase)+"][1]/(np.sum(r_param)**2))*dot(grad(u_d['u_"+str(len(reactants))+"']), grad(v_d['v_"+str(len(reactants))+"']))*dx(1)" 
		#sys.exit()
	
###################
	element = '['
	for k in range(0,len(reactants)+int(number_of_inerts)):
		if (k+1) < len(reactants)+int(number_of_inerts):
			element = element + 'P1,'
		else:
			element = element + 'P1]'

	necessary_dictionary = {'F': F,'gas_num': gas_num, 'surf_num': (len(reactants)-gas_num-1),'element': element,'gas_molecules': gas_molecules,'reactants_number':reactants_number,'reactants':reactants,'molecules_in_gas_phase':molecules_in_gas_phase}
	#print(F)
	#sys.exit()
	return necessary_dictionary, rate_array, rev_irr

#test_new_again = make_f_equation(reactions_test)

#print(test_new_again['reactants'])

def rateEqs(rate_array,rev_irr):

###########################################################################################################
	# Calls the array and 
	rateStrings = {}
	for z in range(0, rate_array.shape[1]):
		rateStrings[z] = ''

	for k,z in enumerate(rate_array):
		
		F = ''


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

		new_neg = 'kVals["kf'+str(k)+'"]'
		for j,v in enumerate(neg):
			#new_neg = new_neg+"*(cat_data['conVtime_"+str(v)+"']**"+str(abs(val_neg[j]))+")"#
			new_neg = new_neg+"*np.power(cat_data['conVtime_"+str(v)+"'],"+str(abs(val_neg[j]))+")"#
		new_pos = 'kVals["kb'+str(k)+'"]'
		for j,v in enumerate(pos):
			#new_pos = new_pos+"*(cat_data['conVtime_"+str(v)+"']**"+str(abs(val_pos[j]))+")"#
			new_pos = new_pos+"*np.power(cat_data['conVtime_"+str(v)+"'],"+str(abs(val_pos[j]))+")"#
		for j,v in enumerate(together):
			F = rateStrings[v]#''
			if j < len(neg):
				#print("neg")
				if rev_irr[k] == 1:
					F = F+"-"+str(abs(val_neg[j]))+" * "+new_pos+" + "+str(abs(val_neg[j]))+"* "+new_neg
					#print("- dk*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
				else:
					F = F+"+"+str(abs(val_neg[j]))+"* "+new_neg
					#print("+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
			else:
				#print("pos")
				if rev_irr[k] == 1:
					F = F+"+"+str(abs(val_pos[j-len(neg)]))+" * "+new_pos+" - "+str(abs(val_pos[j-len(neg)]))+"* "+new_neg
					#print("+ dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
				else:
					F = F+"-"+str(abs(val_pos[j-len(neg)]))+" * "+new_neg
					#print("- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
			
			rateStrings[v] = F

	return rateStrings

def rrmEqs(rate_array,rev_irr,domain):

###########################################################################################################
	# Calls the array and 
	rateStrings = {}
	for z in range(0, rate_array.shape[1]):
		rateStrings[z] = ''

	for k,z in enumerate(rate_array):
		
		F = ''

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

		new_neg = 'r_const["kf'+str(k)+'"]'
		for j,v in enumerate(neg):
			#new_neg = new_neg+"*(cat_data['conVtime_"+str(v)+"']**"+str(abs(val_neg[j]))+")"#
			new_neg = new_neg+"*(u["+str(v)+"]**"+str(abs(val_neg[j]))+")"#
		new_pos = 'r_const["kb'+str(k)+'"]'
		for j,v in enumerate(pos):
			#new_pos = new_pos+"*(cat_data['conVtime_"+str(v)+"']**"+str(abs(val_pos[j]))+")"#
			new_pos = new_pos+"*(u["+str(v)+"]**"+str(abs(val_pos[j]))+")"#
		
		for j,v in enumerate(together):
			F = F+rateStrings[v]
			if j < len(neg):
				#print("neg")
				if rev_irr[k] == 1:
					F = F+"-"+str(abs(val_neg[j]))+" * "+new_pos+" + "+str(abs(val_neg[j]))+"* "+new_neg
					#print("- dk*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
				else:
					F = F+"+"+str(abs(val_neg[j]))+"* "+new_neg
					#print("+ dk*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
			else:
				#print("pos")
				if rev_irr[k] == 1:
					F = F+"+"+str(abs(val_pos[j-len(neg)]))+" * "+new_pos+" - "+str(abs(val_pos[j-len(neg)]))+"* "+new_neg
					#print("+ dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx(1))"+"- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
				else:
					F = F+"-"+str(abs(val_pos[j-len(neg)]))+" * "+new_neg
					#print("- dk*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx(1))")
					#print("")
		
			rateStrings[v] = F

	for allEquations in rateStrings:
		rateStrings[allEquations] = 'assemble( ( inner(' + rateStrings[allEquations] + ', Sw3) )* '+domain+')'

	return rateStrings
