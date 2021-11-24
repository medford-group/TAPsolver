if round(t,6) in output_fitting[legend_label[k_fitting]]['times']:

	c_exp = output_fitting[legend_label[k_fitting]]['values'][output_fitting[legend_label[k_fitting]]['times'].index(round(t,6))]
	slope = (-c_exp)/(1/mesh_size)
	intercept = c_exp - ((1-(1/mesh_size))*slope)
	w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)
	w_new2 = interpolate(w_new,V_du)
	w3 = project(w_new2,V_du)		
	
	try:
		if legend_label[k_fitting] != 'Inert':
			if sigma == None:
				jfunc_2 += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))																	
			else:
				jfunc_2 += assemble(inner(u_n[k_fitting]*to_flux[k_fitting] - w3,u_n[k_fitting]*to_flux[k_fitting] - w3)*dP(1))/((sigma)**2)
			else:
				pass