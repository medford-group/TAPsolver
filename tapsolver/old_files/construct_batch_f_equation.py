def make_batch_equation(reactions_n,reactants_number,arrForward,arrBackward,gForward,temp_change=False):

	tempVariation = temp_change

	active_site_names = ['*','^','@','#']

	F = ''
	rate_array, reactions, reactants, rev_irr = variational_list_parsing(reactions_n)

	gas_num = 0

	gas_num = len(reactants)
	molecules_in_gas_phase = 0

	for k in reactants:
		if '*' in k:
			break
		else:
			molecules_in_gas_phase += 1 
	
	gas_molecules = reactants[:gas_num]
	
	gas_num = len(reactants)

	tail = len(reactants)

	for k in range(0,len(reactants)):
		if k != len(reactants)-1:
			F = F + "((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx + "
		else:
			F = F + "((u_d['u_"+str(k+1)+"'] - u_nd['u_n"+str(k+1)+"']) )*v_d['v_"+str(k+1)+"']*dx "

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
	
		if k in gForward:
			new_neg = '(kbt*constantTemp/hb)*exp(-r_Ga_in["Ga'+str(k)+'"])'

		elif k in arrForward:
			new_neg = 'r_Ao["Aof'+str(k)+'"]*exp(-r_Ea["Eaf'+str(k)+'"]/(Rgas*constantTemp))'
				
		else:
			new_neg = 'r_const["kf'+str(k)+'"]'
				
		for j,v in enumerate(neg):
			new_neg = new_neg+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"#
			
		if k in gForward:
			new_pos = '(kbt*constantTemp/hb)*exp((-(r_Ga_in["Ga'+str(k)+'"]-r_dG_in["dG'+str(k)+'"])))'
			
		elif k in arrBackward:
			new_pos = 'r_Ao["Aob'+str(k)+'"]*exp(-r_Ea["Eab'+str(k)+'"]/(Rgas*constantTemp))'	
		else:
			new_pos = 'r_const["kb'+str(k)+'"]'
				
		for j,v in enumerate(pos):
			new_pos = new_pos+"*(u_d['u_"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"#

		for j,v in enumerate(together):
			if j < len(neg):
				
				if rev_irr[k] == 1:
					F = F+"- dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx)"+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
					
				else:
					F = F+"+ dk*Constant(theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
					
			else:
				if rev_irr[k] == 1:
					F = F+"+ dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx)"+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
					
				else:
					F = F+"- dk*Constant(theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"

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

			if k in gForward:
				new_neg = '(kbt*constantTemp/hb)*exp(-r_Ga_in["Ga'+str(k)+'"])'

			elif k in arrForward:
				new_neg = 'r_Ao["Aof'+str(k)+'"]*exp(-r_Ea["Eaf'+str(k)+'"]/(Rgas*constantTemp))'
				
			else:
				new_neg = 'r_const["kf'+str(k)+'"]'
				
			for j,v in enumerate(neg):
				new_neg = new_neg+"*(u_nd['u_n"+str(v+1)+"']**"+str(abs(val_neg[j]))+")"
			
			if k in gForward:
				new_pos = '(kbt*constantTemp/hb)*exp((-(r_Ga_in["Ga'+str(k)+'"]-r_dG_in["dG'+str(k)+'"])))'
		
			elif k in arrBackward:
				new_pos = 'r_Ao["Aob'+str(k)+'"]*exp(-r_Ea["Eab'+str(k)+'"]/(Rgas*constantTemp))'
				
			else:
				new_pos = 'r_const["kb'+str(k)+'"]'
				
			for j,v in enumerate(pos):
				new_pos = new_pos+"*(u_nd['u_n"+str(v+1)+"']**"+str(abs(val_pos[j]))+")"

			for j,v in enumerate(together):
				if j < len(neg):
					if rev_irr[k] == 1:
						F = F+"- dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx)"+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
						
					else:
						F = F+"+ dk*Constant(1-theta)*("+str(abs(val_neg[j]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
						
				else:
					
					if rev_irr[k] == 1:
						F = F+"+ dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_pos+"*v_d['v_"+str(v+1)+"']*dx)"+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
						
					else:
						F = F+"- dk*Constant(1-theta)*("+str(abs(val_pos[j-len(neg)]))+"* "+new_neg+"*v_d['v_"+str(v+1)+"']*dx)"
						
	element = '['
	for k in range(0,len(reactants)):
		if (k+1) < len(reactants):
			element = element + 'P1,'
		else:
			element = element + 'P1]'

	necessary_dictionary = {'F': F,'gas_num': gas_num, 'surf_num': (len(reactants)-gas_num-1),'element': element,'gas_molecules': gas_molecules,'reactants_number':reactants_number,'reactants':reactants,'molecules_in_gas_phase':molecules_in_gas_phase}
	
	return necessary_dictionary, rate_array, rev_irr