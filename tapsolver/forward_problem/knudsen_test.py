def knudsenTest(species_list,sim_steps,folder,time,points,intensity,fraction):

	user_data = {}
	species_list = species_list[:len(species_list)]#'./experimental_data/'

	for k in range(0,len(species_list)):
		user_data[species_list[k]] = pd.read_csv(folder+'/flux_data/'+species_list[k]+'.csv',header=None)
	print()
	print('Knudsen Regime Fingerprint Test (should be ~ 0.31):')
	curve_fitting = {}
	exp_data = user_data
	for k_num_new, k_new in enumerate(species_list):
		peak_loc = user_data[k_new].iloc[user_data[k_new][1].idxmax()]
		near_peak = peak_loc[0]/(time/sim_steps)
		peak2 = user_data[k_new].loc[user_data[k_new][0] == peak_loc[0]].index

		print(k_new+': '+str(peak_loc[0]*peak_loc[1]/(intensity*float(fraction[k_num_new]))))
		print()
