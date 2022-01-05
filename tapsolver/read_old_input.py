def readInput(sim_file,inputForm = 'old'):
	
	"""
	Convert the input file into dictionaries for TAPsolver to use
	"""
	if inputForm == 'old':
	
		user_data = pd.read_csv(sim_file,header=None)
		
		rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
		rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(user_data == 'Reaction_Information')
	
		reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
		feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
		reaction_info = user_data.iloc[1+rows_4[0]:,:]
	
		reactor_kinetics_input = {}
		
		for k in range(0,len(reactor_info.index)):
			try:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
	
	
		for k in range(0,len(feed_surf_info.index)):
			try:
				reactor_kinetics_input[feed_surf_info.iloc[k,0]] = float(feed_surf_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[feed_surf_info.iloc[k,0]] = feed_surf_info.iloc[k,1]
	
	
		reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']*(reactor_kinetics_input['Catalyst Location'] - 0.5*(reactor_kinetics_input['Catalyst Fraction']))#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
		reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length'] - (reactor_kinetics_input['len_cat'] + reactor_kinetics_input['len_inert_1'])#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	
		reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()

	else:
		
		user_data = pd.read_csv(sim_file,header=None)
		
		rows_1, cols_1 = np.where(user_data == 'Reactor_Information')
		rows_2, cols_2 = np.where(user_data == 'Feed_&_Surface_Composition')
		rows_4, cols_4 = np.where(user_data == 'Reaction_Information')

		linkedKinetics = False
		if user_data[0].str.contains('Linked Kinetics').any():
			linkedKinetics = True
			rows_5, cols_5 = np.where(user_data == 'Linked Kinetics')
		thermoConstraints = False
		if user_data[0].str.contains('Thermodynamic Constraints').any() and linkedKinetics == True:
			thermoConstraints = True
			rows_6, cols_6 = np.where(user_data == 'Thermodynamic Constraints')    		
		elif user_data[0].str.contains('Thermodynamic Constraints').any():
			thermoConstraints = True
			rows_5, cols_5 = np.where(user_data == 'Thermodynamic Constraints')    		
		reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:]
		feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
	#	data_storage = user_data.iloc[1+rows_3[0]:rows_4[0]-1,:]
		if thermoConstraints == False and linkedKinetics == False:
			reaction_info = user_data.iloc[1+rows_4[0]:,:]
		elif linkedKinetics == True and thermoConstraints == False:
			reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			linked_kinetics = user_data.iloc[1+rows_5[0]:,:]
		elif linkedKinetics == False and thermoConstraints == True:		
			reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			thermo_constraints = user_data.iloc[1+rows_5[0]:,:]
		elif linkedKinetics == True and thermoConstraints == True:
			reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
			linked_kinetics = user_data.iloc[(1+rows_5[0]):(rows_6[0]-1),:]
			thermo_constraints = user_data.iloc[1+rows_6[0]:,:]		


		#thermoConstraints = False
		#if user_data[0].str.contains('Thermodynamic Constraints').any():
		#	thermoConstraints = True
		#	rows_5, cols_5 = np.where(user_data == 'Thermodynamic Constraints')    		
		#
		#reactor_info = user_data.iloc[(1+rows_1[0]):(rows_2[0]-1),:] 
		#feed_surf_info = user_data.iloc[1+rows_2[0]:rows_4[0]-1,:]
		#
		#if thermoConstraints == False:
		#	reaction_info = user_data.iloc[1+rows_4[0]:,:]
		#else:
		#	reaction_info = user_data.iloc[(1+rows_4[0]):(rows_5[0]-1),:]
		#	thermo_constraints = user_data.iloc[1+rows_5[0]:,:]

		reactor_kinetics_input = {}
		
		number_of_gasses = 0
		number_of_surface = 0
		for z in feed_surf_info.iloc[0,1:]:
			if type(z) == (str):
				number_of_gasses += 1 
		for z in feed_surf_info.iloc[5,1:]:
			if type(z) == (str):
				number_of_surface += 1 

		gas_species = feed_surf_info.iloc[0,1:number_of_gasses+1]
		surface_species = feed_surf_info.iloc[6,1:number_of_surface+1]
		
		reactor_kinetics_input['Number of Inerts'] = 0
		reactor_kinetics_input['Pulse Size'] = ''
		reactor_kinetics_input['Pulse Time'] = ''
		reactor_kinetics_input['Mass List'] = ''
		reactor_kinetics_input['Number of Reactants'] = 0

		for jnum,j in enumerate(gas_species):
			
			if j.find('Inert') == 0:
				reactor_kinetics_input['Number of Inerts'] += 1
			if jnum == len(gas_species)-1:

				reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])
				reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])
				reactor_kinetics_input['Mass List'] = reactor_kinetics_input['Mass List']+str(feed_surf_info.iloc[3,1+jnum])
				
			else:
				reactor_kinetics_input['Pulse Size'] = reactor_kinetics_input['Pulse Size']+str(feed_surf_info.iloc[1,1+jnum])+','
				reactor_kinetics_input['Pulse Time'] = reactor_kinetics_input['Pulse Time']+str(feed_surf_info.iloc[2,1+jnum])+','
				reactor_kinetics_input['Mass List'] = reactor_kinetics_input['Mass List']+str(feed_surf_info.iloc[3,1+jnum])+','
			reactor_kinetics_input['Number of Reactants'] += 1
		reactor_kinetics_input['Initial Surface Composition'] = ''

		for jnum,j in enumerate(surface_species):

			if jnum == len(surface_species)-1:
				reactor_kinetics_input['Initial Surface Composition'] = reactor_kinetics_input['Initial Surface Composition']+str(feed_surf_info.iloc[6,1+jnum])
			else:
				reactor_kinetics_input['Initial Surface Composition'] = reactor_kinetics_input['Initial Surface Composition']+str(feed_surf_info.iloc[6,1+jnum])+','

		for k in range(0,len(reactor_info.index)):
			try:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = float(reactor_info.iloc[k,1])
			except ValueError:
				reactor_kinetics_input[reactor_info.iloc[k,0]] = reactor_info.iloc[k,1]
		
		reactor_kinetics_input['Reactor Length'] =  float(reactor_info.iloc[0,1]) + float(reactor_info.iloc[0,2]) + float(reactor_info.iloc[0,3])		
		
		reactor_kinetics_input['Catalyst Fraction'] = float(reactor_info.iloc[0,2])/reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['Catalyst Location'] = (float(reactor_info.iloc[0,1])+(float(reactor_info.iloc[0,2])/2))/reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['Void Fraction Inert'] = float(reactor_info.iloc[1,1])
		reactor_kinetics_input['Void Fraction Catalyst'] = float(reactor_info.iloc[1,2])

		reactor_kinetics_input['len_inert_1'] = reactor_kinetics_input['Reactor Length']*(reactor_kinetics_input['Catalyst Location'] - 0.5*(reactor_kinetics_input['Catalyst Fraction']))#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
		reactor_kinetics_input['len_cat'] = (reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length'] 
		reactor_kinetics_input['len_inert_2'] =  reactor_kinetics_input['Reactor Length'] - (reactor_kinetics_input['len_cat'] + reactor_kinetics_input['len_inert_1'])#reactor_kinetics_input['Reactor Length']/2 -  0.5*(reactor_kinetics_input['Catalyst Fraction'])*reactor_kinetics_input['Reactor Length']
	
		reactor_kinetics_input['reactions_test'] = reaction_info.iloc[:,0].tolist()


	
	kinetic_parameters = {}
	Ao = {}
	Ea = {}
	Ga = {}
	dG = {}
	link = {}

	fittingParametersList = []

	gForward = []
	arrForward = []
	arrBackward = []
	linkForward = []
	linkBackward = []

	for j in range(0,len(reaction_info.index)):
		if reaction_info.iloc[j,1].find("#") > 0:
			Anew, Eanew = reaction_info.iloc[j,1].split("#")
			if Anew.find("!") < 0:
				#!?
				if Anew.find("{") == 0:
					link['Ga'+str(j)] = float(Anew)
					fittingParametersList.append('Ga'+str(j))
				else:
					Ga['Ga'+str(j)] = float(Anew)
					kinetic_parameters['Ga'+str(j)] = float(Anew)
					fittingParametersList.append('Ga'+str(j))
				
			else:
				#!?
				if Anew.find("{") == 0:
					fittingParametersList.append(Anew)
				else:
					Ga['Ga'+str(j)] = float(Anew[:-1])
			if Eanew.find("!") < 0:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					dG['dG'+str(j)] = float(Eanew)
					fittingParametersList.append('dG'+str(j))
			
			else:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					dG['dG'+str(j)] = float(Eanew[:-1])
			gForward.append(j)
		elif reaction_info.iloc[j,1].find("$") > 0:
			Anew, Eanew = reaction_info.iloc[j,1].split("$")
			if Anew.find("!") < 0:
				#!?
				if Anew.find("{") == 0:
					fittingParametersList.append(Anew)
				else:	
					Ao['Aof'+str(j)] = float(Anew)
					fittingParametersList.append('Aof'+str(j))
			
			else:
				#!?
				if Anew.find("{") == 0:
					fittingParametersList.append(Anew)
				else:
					Ao['Aof'+str(j)] = float(Anew[:-1])
			if Eanew.find("!") < 0:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					Ea['Eaf'+str(j)] = float(Eanew)
					fittingParametersList.append('Eaf'+str(j))
			
			else:
				#!?
				if Eanew.find("{") == 0:
					fittingParametersList.append(Eanew)
				else:
					Ea['Eaf'+str(j)] = float(Eanew[:-1])
			arrForward.append(j)
		else:
			if reaction_info.iloc[j,1].find("!") < 0:
				if reaction_info.iloc[j,1].find("{") == 0:
					link['kf'+str(j)] = reaction_info.iloc[j,1][1:-1]
					if reaction_info.iloc[j,1] not in fittingParametersList:
						fittingParametersList.append(reaction_info.iloc[j,1])
					linkForward.append(j)
				else:
					kinetic_parameters['kf'+str(j)] = float(reaction_info.iloc[j,1])
					fittingParametersList.append('kf'+str(j))
			
			else:
				if reaction_info.iloc[j,1].find("{") == 0:
					link['kf'+str(j)] = reaction_info.iloc[j,1][1:-1]
					if reaction_info.iloc[j,1] not in fittingParametersList:
						fittingParametersList.append(reaction_info.iloc[j,1])
					linkForward.append(j)
				else:
					new_value = float(reaction_info.iloc[j,1])
					kinetic_parameters['kf'+str(j)] = new_value#float(reaction_info.iloc[j,1])
		if str(reaction_info.iloc[j,2]) != 'nan':
			if str(reaction_info.iloc[j,2]).find("$") > 0:
				Anew, Eanew = str(reaction_info.iloc[j,2]).split("$")
				if Anew.find("!") < 0:
					#!?
					if Anew.find("{") == 0:
						fittingParametersList.append(Anew)
					else:
						Ao['Aob'+str(j)] = float(Anew)
						fittingParametersList.append('Aob'+str(j))
					
				else:
					#!?
					if Anew.find("{") == 0:
						fittingParametersList.append(Anew)
					else:
						Ao['Aob'+str(j)] = float(Anew[:-1])						
				if Eanew.find("!") < 0:
					#!?
					if Eanew.find("{") == 0:
						fittingParametersList.append(Eanew)
					else:
						Ea['Eab'+str(j)] = float(Eanew)
						fittingParametersList.append('Eab'+str(j))
					
				else:
					#!?
					if Eanew.find("{") == 0:
						fittingParametersList.append(Eanew)
					else:
						Ea['Eab'+str(j)] = float(Eanew[:-1])
				arrBackward.append(j)
			else:
				if str(reaction_info.iloc[j,2]).find("!") < 0:
				
					if reaction_info.iloc[j,2].find("{") == 0:
						link['kb'+str(j)] = reaction_info.iloc[j,2][1:-1]
						if reaction_info.iloc[j,2] not in fittingParametersList:
							fittingParametersList.append(reaction_info.iloc[j,2])
						linkBackward.append(j)
					else:
						kinetic_parameters['kb'+str(j)] = float(reaction_info.iloc[j,2])
						fittingParametersList.append('kb'+str(j))
				else:
					if reaction_info.iloc[j,2].find("{") == 0:
						link['kb'+str(j)] = reaction_info.iloc[j,2][1:-2]
						if reaction_info.iloc[j,2] not in fittingParametersList:
							fittingParametersList.append(reaction_info.iloc[j,2])
						linkBackward.append(j)
					else:
						new_value = float(reaction_info.iloc[j,2][:-1])
						kinetic_parameters['kb'+str(j)] = new_value
		else:
			pass


	kin_in = kinetic_parameters.copy()
	Ao_in = Ao.copy()
	Ea_in = Ea.copy()
	Ga_in = Ga.copy()
	dG_in = dG.copy()

	linkedParameters = {}
	if linkedKinetics == True:
		for j in range(0,len(linked_kinetics.index)):
			linkedParameters[linked_kinetics.iloc[j,0]] = float(linked_kinetics.iloc[j,1])
	reactor_kinetics_input['linked parameters'] = linkedParameters
	reactor_kinetics_input['linked names'] = link
	reactor_kinetics_input['link forward'] = linkForward
	reactor_kinetics_input['link backward'] = linkBackward

	thermo_equations = []
	thermo_values = []
	if thermoConstraints == True:
		for j in range(0,len(thermo_constraints.index)):
			thermo_equations.append(thermo_constraints.iloc[j,0])
			try:
				if type(float(thermo_constraints.iloc[j,1])) == float and math.isnan(float(thermo_constraints.iloc[j,1])) == False:
					thermo_values.append(thermo_constraints.iloc[j,1])
				else:
					thermo_values.append('none')
			except:
				print(type(thermo_constraints.iloc[j,1]))
				pass

	reactor_kinetics_input['thermo equations'] = thermo_equations
	reactor_kinetics_input['thermo values'] = thermo_values

	return reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,fittingParametersList,arrForward,arrBackward
