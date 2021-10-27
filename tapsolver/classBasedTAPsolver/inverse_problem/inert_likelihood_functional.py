				if reac_input['Fit Inert'].lower() == 'true':
					
					for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
						if objSpecies[k_fitting] == '1':
							if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:
								c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
								slope = (-c_exp)/(1/mesh_size)
								intercept = c_exp - ((1-(1/mesh_size))*slope)
								w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)
								w_new2 = interpolate(w_new,V_du)
								w3 = project(w_new2,V_du)
	
								try:
									jfunc_2 += assemble(inner(u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3,u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3)*dP(1))								
									
								except UnboundLocalError:
									jfunc_2 = assemble(inner(u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3,u_n[necessary_values['surf_num']+k_fitting-2]*to_flux[k_fitting] - w3)*dP(1))
				