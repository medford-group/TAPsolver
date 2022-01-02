					if conversion == True:
											
						w_new = Expression('1',degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)
	
						try:
							if legend_label[2] != 'Inert':
								jfunc_2 += assemble(inner(dt*u_n[2]*to_flux[2],w3)*dP(1))							
							else:
								pass
	
						except UnboundLocalError:
							if legend_label[2] != 'Inert':
								jfunc_2 = assemble(inner(dt*u_n[2]*to_flux[2],w3)*dP(1))
							else:
								pass