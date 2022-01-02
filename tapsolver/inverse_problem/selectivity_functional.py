					if selectivity == True:
						w_new = Expression('1',degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)
	
						try:
							if legend_label[2] != 'Inert':
								selectDenom = 0
								for zap in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
									selectDenom += u_n[zap]*to_flux[zap]
								jfunc_2 += assemble(inner(u_n[2]*to_flux[2]/selectDenom,w3)*dP(1))							
							else:
								pass
	
						except UnboundLocalError:
							if legend_label[2] != 'Inert':
								selectDenom = 0
								for zap in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
									selectDenom += u_n[zap]*to_flux[zap]
								jfunc_2 = assemble(inner(u_n[2]*to_flux[2]/selectDenom,w3)*dP(1))
							else:
								pass
						