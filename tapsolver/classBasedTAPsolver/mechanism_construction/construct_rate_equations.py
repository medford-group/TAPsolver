def rateEqs(rate_array,rev_irr,gForward,arrForward,arrBackward):

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


		if k in gForward:#,Ga_in,dG_in
			new_neg = '(kbt*constantTemp/hb)*exp(-Ga_in["Ga'+str(k)+'"])'

		elif k in arrForward:
			new_neg = 'r_Ao["Aof'+str(k)+'"]*exp(-r_Ea["Eaf'+str(k)+'"]/(Rgas*constantTemp))'
				
		else:
			new_neg = 'constants_input["kf'+str(k)+'"]'
		
		for j,v in enumerate(neg):
			new_neg = new_neg+"*np.power(cat_dataRate['convtime_"+str(v)+"'],"+str(abs(val_neg[j]))+")"#

		if k in gForward:
			new_pos = '(kbt*constantTemp/hb)*exp((-(Ga_in["Ga'+str(k)+'"]-dG_in["dG'+str(k)+'"])))'
		
		elif k in arrBackward:
			new_pos = 'r_Ao["Aob'+str(k)+'"]*exp(-r_Ea["Eab'+str(k)+'"]/(Rgas*constantTemp))'
				
		else:
			new_pos = 'constants_input["kb'+str(k)+'"]'
		
		for j,v in enumerate(pos):
			new_pos = new_pos+"*np.power(cat_dataRate['convtime_"+str(v)+"'],"+str(abs(val_pos[j]))+")"#
		for j,v in enumerate(together):
			F = rateStrings[v]#''
			if j < len(neg):
				#print("neg")
				if rev_irr[k] == 1:
					F = F+"-"+str(abs(val_neg[j]))+" * "+new_pos+" + "+str(abs(val_neg[j]))+"* "+new_neg
				else:
					F = F+"+"+str(abs(val_neg[j]))+"* "+new_neg
			else:
				if rev_irr[k] == 1:
					F = F+"+"+str(abs(val_pos[j-len(neg)]))+" * "+new_pos+" - "+str(abs(val_pos[j-len(neg)]))+"* "+new_neg
				else:
					F = F+"-"+str(abs(val_pos[j-len(neg)]))+" * "+new_neg
			
			rateStrings[v] = F

	return rateStrings