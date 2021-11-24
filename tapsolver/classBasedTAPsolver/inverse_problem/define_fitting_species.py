def curveFitting(pulse_time,original_TAPobject):
	#def curveFitting(species_list,sim_steps,folder,timeTot,points,objSpecies):
	frequency = 1
	"""Define the objective function for optimizing kinetic parameters"""
	time_steps = pulse_time*1000
	syn_time_step = pulse_time/time_steps

	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1
		
	user_data = original_TAPobject.experimental_data


	curve_fitting = {}
	#self.inert_gasses_objective
	for k_new in original_TAPobject.gasses_objective: #  and original_TAPobject.inert_gasses_objective
		
		def find_experimental_point(n,exp_step):
			""" Find an appropriate intensity point for the fitting process """

			approx_exp_n = n*(syn_time_step)/exp_step
			if approx_exp_n != n:
				high = math.ceil(approx_exp_n)
				low = int(approx_exp_n)
				return interp(user_data[k_new]['0'][high-1],user_data[k_new][1][low-1],user_data[k_new][0][high-1],user_data[k_new][0][low-1],n*(syn_time_step))

			else:
				user_data[k_new][1][n]
				return user_data[k_new][1][n]

		def exp_point_to_syn_point(n_exp,exp_step):
			"""Align an experimental data point with the associated (or nearest synthetic point)"""
			approx_syn_n = n_exp*exp_step/(syn_time_step)
				
			if int(approx_syn_n) > 0:
				return find_experimental_point(int(approx_syn_n),exp_step)
			else:
				return find_experimental_point(math.ceil(approx_syn_n),exp_step) 

		time_step = []
		times = []
		values = []

		for j in range(0,len(user_data['time']['0'])):
			user_data['time']['0'][j] = round(user_data['time']['0'][j],5) 
		
		exp_time_step = user_data['time']['0'][len(user_data['time']['0'])-2] - user_data['time']['0'][len(user_data['time']['0'])-3]
		
		fitStartTime = 0
		frequency = 1

		for k in range(fitStartTime,int(time_steps),frequency):
			time_step.append(k)
			times.append(k*(syn_time_step))
			#current_list_value = user_data['time']['0'].index(round(k*(syn_time_step),5))
			#print(current_list_value)
			if round(exp_time_step,5) == syn_time_step:
				values.append(user_data[k_new]['0'][k])
			else:
				values.append(find_experimental_point(k,exp_time_step))#k*(syn_time_step)
			#t+=	
			#print(find_experimental_point(k*(syn_time_step),exp_time_step))

		data = {}
		data['time_step'] = time_step
		data['times'] = times
		data['values'] = values

		curve_fitting[k_new] = data

	return curve_fitting
