def pointFitting(species_list,sim_steps,folder,time,points,objSpecies):

	"""Define the objective function for optimizing kinetic parameters"""

	syn_time_step = time/sim_steps
	
	def interp(y2,y1,x2,x1,xn):
		"""Simple linear interpolation function"""
		return ((y2 - y1)/(x2 - x1))*(xn - x1) + y1

	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'

	for k in range(0,len(species_list)):
		if objSpecies[k] == '1':
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
					return interp(user_data[k_new][1][high],user_data[k_new][1][low],user_data[k_new][0][high],user_data[k_new][0][low],n*(syn_time_step))

				else:
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
			exp_time_step = user_data[k_new][0][1] - user_data[k_new][0][0]

			near_start = round(user_data[k_new].iloc[30,0],6)/(time/sim_steps)
		
			peak_loc = user_data[k_new].iloc[user_data[k_new][1].idxmax()]
			near_peak = peak_loc[0]/(time/sim_steps)
			peak2 = user_data[k_new].loc[user_data[k_new][0] == peak_loc[0]].index

			test3 = int(round((peak2[0]+1)/2,0))
			mid_loc = user_data[k_new].iloc[test3,:]
			near_mid = mid_loc[0]/(time/sim_steps)

			time_step.append(int(near_mid))
			times.append(int(near_mid)*(syn_time_step))

			values.append(find_experimental_point(int(near_mid),exp_time_step))
			
			if points > 1:
				time_step.append(int(near_peak))
				times.append(int(near_peak)*(syn_time_step))
				values.append(find_experimental_point(int(near_peak),exp_time_step))

			if points > 2:
				value_test = 0.9*peak_loc[1]
				sort_3 = user_data[k_new].iloc[(user_data[k_new][1]-value_test).abs().argsort()[:]]

				for k in range(0,sort_3.index[0]):
					if sort_3.iloc[k,0] > peak_loc[0]:
						thr_point = sort_3.iloc[k]
						break
					else:
						pass

				near_3 = thr_point[0]/(time/sim_steps)

				time_step.append(int(near_3))
				times.append(int(near_3)*(syn_time_step))
				values.append(find_experimental_point(int(near_3),exp_time_step))
			
			if points > 3:
				value_test = 0.75*peak_loc[1]
				sort_4 = user_data[k_new].iloc[(user_data[k_new][1]-value_test).abs().argsort()[:]]

				for k in range(0,sort_4.index[0]):
					if sort_4.iloc[k,0] > peak_loc[0]:
						four_point = sort_4.iloc[k]
						break
					else:
						pass

				near_4 = four_point[0]/(time/sim_steps)

				time_step.append(int(near_4))
				times.append(int(near_4)*(syn_time_step))
				values.append(find_experimental_point(int(near_4),exp_time_step))
				

			data = {}
			data['time_step'] = time_step
			data['times'] = times
			data['values'] = values

			curve_fitting[k_new] = data
			
	return curve_fitting