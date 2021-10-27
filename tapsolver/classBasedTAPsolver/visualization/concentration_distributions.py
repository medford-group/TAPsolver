def concDistPlot(input_file = './input_file.csv',dataType='surf',output_name='./cat.gif',pulse=1):
	
	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm = 'new')
	reac_input = reactor_kinetics_input
	fit_temperature = False
	reac_input['Advection'] = 'FALSE'
	reac_input['Reactor Type'] = 'tap'
	reac_input['Scale Output'] = 'FALSE'
	r_links = reactor_kinetics_input['linked parameters'].copy()
	print(r_links)
	linkForward = reactor_kinetics_input['link forward'].copy()
	linkBackard = reactor_kinetics_input['link backward'].copy()
	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,reac_input['linked names'],linkForward,linkBackard,fit_temperature)
	#fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),reac_input['Scale Output'])

	reac_input['Advection'] = 'FALSE'
	reac_input['Reactor Type'] = 'tap'
	reac_input['Scale Output'] = 'FALSE'
	fit_temperature = False

	try:
		dataLocation = reac_input['Output Folder Name']+'_folder/thin_data/'
		dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)
	except:
		dataLocation = reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(pulse)+'/'
		dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)

	dataDictionary = {}

	for knum,k in enumerate(necessary_values['reactants']):
		dataDictionary[k] = pd.read_csv(dataLocation+k+'.csv',header=None)

	dfColumns = pd.read_csv(dataLocation+necessary_values['reactants'][0]+'.csv',header=None)

	fig,ax = plt.subplots()
	for knum,k in enumerate(necessary_values['reactants']):
		if dataType == 'surf':
			if knum > necessary_values['molecules_in_gas_phase']:
				ax.plot(dataDictionary[k].iloc[int(len(dfColumns.columns)/2),:],label = k)
		elif dataType == 'gas':
			if knum < necessary_values['molecules_in_gas_phase']:
				ax.plot(dataDictionary[k].iloc[dataDictionary[k].iloc[:,int(len(dfColumns.columns)/2)],:],label = k)
		else:
			if necessary_values['reactants'] in dataType:
				ax.plot(dataDictionary[k].iloc[dataDictionary[k].iloc[:,int(len(dfColumns.columns)/2)],:],label = k)

	plt.show()
