
def call_sim(input_file, kForward, kBackward):
	kForward = kForward.strip('[]').split(',')
	kForward = [float(i) for i in kForward]
	kBackward = kBackward.strip('[]').split(',')
	kBackward = [float(i) for i in kBackward]
	reactor_kinetics_input,kinetic_parameters,kin_in = read_input(input_file)

	for i in np.arange(start = 0, stop = len(kForward)):
		kin_in['kf' + str(i)] = kForward[i]
		if kBackward[i] == 0.0:
			continue
		else:
			kin_in["kb" + str(i)] = kBackward[i]
	kinetic_parameters = kin_in

	print(kinetic_parameters)

	# kin_in['kf0'] = kForward
	# kinetic_parameters['kf0'] = kForward
	# kin_in['kb0'] = kBackward
	# kinetic_parameters['kb0'] =  kBackward
	# # need to determine the gas names wrt to input kinetics rather than <-> or ->
	# if len(kForward) == 1:
	# 	if kBackward[0] == 0:
	# 		reactor_kinetics_input['reactions_test'] = ['CO + * -> CO*']
	# 	else:
	# 		reactor_kinetics_input['reactions_test'] = ['CO + * <-> CO*']
	# elif len(kForward) == 2:
	# 	reactor_kinetics_input['reactions_test'] = ['CO + * <-> CO*', 'CO* + O* -> CO2 + 2*']
	# else: 
	# 	reactor_kinetics_input['reactions_test'] = ['CO + * <-> CO*', 'CO* + O* -> CO2 + 2*']
	if reactor_kinetics_input['Sensitivity Analysis'].lower() == 'true' or reactor_kinetics_input['RRM Analysis'].lower() == 'true' or reactor_kinetics_input['Uncertainty Quantification'].lower() == 'true':
		for parameters in kinetic_parameters:
			reactor_kinetics_input,kinetic_parameters,kin_in = read_input()
			reactor_kinetics_input['Sensitivity Parameter'] = parameters
			reactor_kinetics_input['Display Graph'] = 'FALSE'
			if reactor_kinetics_input['Fit Parameters'].lower() == 'true' or reactor_kinetics_input['Fit Inert'].lower() == 'true':
				print('')
				print('')

				print('Running the Sensitivity/RRM Analysis and Parameter Fitting methods simultaniously is not possible due to conflicts between the tangent linear and adjoint methods.')
				print('Please run again with one of these methods excluded.')
				sys.exit()
			print('')
			print('Processing '+parameters)
			graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)
	else:
		graph_data, legend_label,in_reactants = tap_simulation_function(reactor_kinetics_input,kinetic_parameters)

	# inertValues = pd.read_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/Inert-1.csv')
	# # assume only one reactant for the time being
	# reactantValue = pd.read_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'NH3.csv')
	# productValue = pd.read_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'N.csv')
	# productValue2 = pd.read_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'H.csv')
	# tempInert = inertValues.to_numpy()[:,1]
	# gasConNH3 = reactantValue
	# rateNH3 = inertValues
	# rateH = productValue2
	# rateN = productValue
	# for i in np.arange(1,reactantValue.shape[1]):
	# 	tempReactant = reactantValue.to_numpy()[:,1]
	# 	tempProduct1 = productValue.to_numpy()[:,1]
	# 	tempProduct2 = productValue2.to_numpy()[:,1]
	# 	tempData = yProcedure(tempReactant, 40, isProduct = False, inertPulse = tempInert, inertMass = 40, timeStep = 0.00025, inertZone1 = 2.4375, catalystZone = 0.125, smoothing = 2)
	# 	gasConNH3[i] = tempData['gasConcentration']
	# 	rateNH3[i] = tempData['reactionRate']
	# 	tempData = yProcedure(tempProduct1, 40, isProduct = True, inertPulse = tempInert, inertMass = 40, timeStep = 0.00025, inertZone1 = 2.4375, catalystZone = 0.125, smoothing = 2)
	# 	rateN[i] = tempData['reactionRate']
	# 	tempData = yProcedure(tempProduct2, 40, isProduct = True, inertPulse = tempInert, inertMass = 40, timeStep = 0.00025, inertZone1 = 2.4375, catalystZone = 0.125, smoothing = 2)
	# 	rateH[i] = tempData['reactionRate']

	# rateNH3.to_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'rateNH3.csv')
	# gasConNH3.to_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'gasConNH3.csv')
	# rateH.to_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'rateH.csv')
	# rateN.to_csv(reactor_kinetics_input['Output Folder Name']+'_folder/flux_data/' + 'rateN.csv')
	


call_sim(sys.argv[1], sys.argv[2], sys.argv[3])