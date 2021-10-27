def curveFitting(species_list,sim_steps,folder,timeTot,points,objSpecies):
	frequency = 1
	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = timeTot/sim_steps

	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1
		
	user_data = {}

	if species_list[0].find('Inert') == 0:
				
		for klabel,k in enumerate(objSpecies):

			if objSpecies[klabel] == '1':
				try:
					fitStartValue = False
					user_data[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
				except:
					fitStartValue = True
					user_data[species_list[klabel]] = pd.read_csv(folder+'/flux_data/'+species_list[klabel]+'.csv',header=None)
	else:
		species_list = species_list[:len(species_list)]#'./experimental_data/'
		for k in range(0,len(species_list)):
			if objSpecies[k] == '1':
				try:
					fitStartValue = False
					user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
				except:
					fitStartValue = True
					user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	curve_fitting = {}
	exp_data = user_data
	
	for k_newNum, k_new in enumerate(species_list):

		if objSpecies[k_newNum] == '1':
			def find_experimental_point(n,exp_step):
				""" Find an appropriate intensity point for the fitting process """

				approx_exp_n = n*(syn_time_step)/exp_step
				if approx_exp_n != n:
					high = math.ceil(approx_exp_n)
					low = int(approx_exp_n)
					return interp(user_data[k_new][1][high-1],user_data[k_new][1][low-1],user_data[k_new][0][high-1],user_data[k_new][0][low-1],n*(syn_time_step))

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

			exp_time_step = user_data[k_new][0][user_data[k_new].shape[0]-2] - user_data[k_new][0][user_data[k_new].shape[0]-3]
			near_start = round(user_data[k_new].iloc[2,0],6)/(timeTot/sim_steps)
			if fitStartValue == True:
				fitStartTime = 0
				frequency = 1
			else: 
				fitStartTime = 2
				frequency = 1

			for k in range(fitStartTime,int(sim_steps),frequency):
				time_step.append(k)
				times.append(k*(syn_time_step))
				#print(exp_time_step)

				if round(exp_time_step,5) == syn_time_step:
					values.append(user_data[k_new][1][k])
				else:
					values.append(find_experimental_point(k,exp_time_step))#k*(syn_time_step)
				
				#print(find_experimental_point(k*(syn_time_step),exp_time_step))

			data = {}
			data['time_step'] = time_step
			data['times'] = times
			data['values'] = values

			curve_fitting[k_new] = data

	return curve_fitting