def flux_graph(input_file = './input_file.csv',pulse=None,dispExper=False,dispAnalytic=False,dispObjective=False,show_graph=True,store_graph=False,output_name='./flux.png'):

	timeFunc = 0.4

	reactor_kinetics_input,kinetic_parameters,kin_in,Ao_in,Ea_in,Ga_in,dG_in,gForward,kin_fit,arrForward,arrBackward = readInput(input_file,inputForm = 'new')

	reac_input = reactor_kinetics_input
	reac_input['Pulse Duration'] = timeFunc
	reac_input['Infinite Inert'] = 'FALSE'
	reac_input['Display Experimental Data'] = 'False'
	reac_input['Display Objective Points'] = 'FALSE'
	reac_input['Objective Points'] = 'all'

	if dispExper != False:
		reac_input['Display Experimental Data'] = 'TRUE'
		reac_input['Display Objective Points'] = 'FALSE'
		reac_input['Objective Points'] = 'all'

	if dispObjective != False:
		reac_input['Display Objective Points'] = 'TRUE'
		reac_input['Objective Points'] = 'all'
#
	if dispAnalytic != False:
		reac_input['Infinite Inert'] = 'TRUE'


	r_links = reactor_kinetics_input['linked parameters'].copy()
	print(r_links)
	linkForward = reactor_kinetics_input['link forward'].copy()
	linkBackard = reactor_kinetics_input['link backward'].copy()

	reac_input['Advection'] = 'FALSE'
	reac_input['Time Steps'] = reac_input['Pulse Duration']*1000
	dt = 0.001
	reac_input['Reference Pulse Size'] = 1

	fit_temperature = False
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,reac_input['linked names'],linkForward,linkBackard,fit_temperature)

	reactor_kinetics_input['Reactor Type'] = 'tap'
	fig2,ax2,legend_label,header,colors = establishOutletGraph(reac_input['Reactor Type'],necessary_values['molecules_in_gas_phase'],necessary_values['reactants'],int(reac_input['Number of Inerts']),'FALSE')

	monitored_gas = necessary_values['molecules_in_gas_phase']
	reac_input['Scale Output'] = 'FALSE'

	reac_input['Objective Species'] = '1,1,1,0,0'

	if '/' in str(reac_input['Pulse Size']):
		list_of_feed = reac_input['Pulse Size'].split('/')
		list_of_time = reac_input['Pulse Time'].split('/')
		pulse_variation = len(list_of_feed)
		pulse_time = len(list_of_time)
		reactant_feed = []
		reactant_time = []
		for k in range(0,len(reac_input['Pulse Size'].split('/'))):
			reactant_feed.append(list_of_feed[k].split(','))
			reactant_time.append(list_of_time[k].split(','))
	else:
		pulse_variation = 1
		reactant_feed = []
		reactant_time = []
		
		try:
			reactant_feed.append(reac_input['Pulse Size'].split(','))
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
		except AttributeError:
			reactant_feed.append(reac_input['Pulse Size'])
			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
	
	species_pulse_list = reactant_feed[0]
	current_reac_set = 1
	if type(reactant_time[current_reac_set-1]) == list:
		species_time = reactant_time[current_reac_set-1]
	else:
		species_time = reactant_time

	if str(reac_input['Objective Species']).find(',') != -1:
		objSpecies = list(reac_input['Objective Species'].split(','))
	else:
		print(reac_input['Objective Species'])
		objSpecies = [str(int(reac_input['Objective Species']))]

	if reac_input['Experimental Data Folder'].lower() != 'none':
		output_fitting = curveFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)

	r_param, dx_r, dx2_r, frac_length, cat_location = establishMesh(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['Mesh Size'])
	cat_location = 1 - reac_input['Catalyst Location']
	dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])
	eb = np.array((reac_input['Void Fraction Inert'],reac_input['Void Fraction Catalyst'],reac_input['Void Fraction Inert']))
	ca = (reac_input['Reactor Radius']**2)*3.14159 

	# Define: Pulse
	point_volume = dx_r * ca * eb[0]
	Inert_pulse_conc = reac_input['Reference Pulse Size']/(point_volume)
		
	# Define the diffusion constants
	ref_rate = np.append(reac_input['Reference Diffusion Inert'],reac_input['Reference Diffusion Catalyst'])
	ref_rate = np.append(ref_rate,ref_rate[0])  	
	D = np.empty((len(reac_input['Mass List'].split(',')),3))
	
	def diff_func(ref_mass,ref_T,mol_mass,ref_r):
		return ref_r*(mp.sqrt(ref_mass*reac_input['Reactor Temperature'])/mp.sqrt(ref_T*mol_mass))
	Dout = []
	Din = []
	constantTemp = reac_input['Reactor Temperature']
	constantTemp = Constant(constantTemp)
	testTemp = reac_input['Reactor Temperature']
	testTemp = Constant(testTemp)

	if fit_temperature == True:
		controls = []
		new = Control(constantTemp)
		controls.append(new)

	compMass_list = reac_input['Mass List'].split(',')

	for k,j in enumerate(reac_input['Mass List'].split(',')):
		for k_2,j_2 in enumerate(ref_rate):
			D[k,k_2] = diff_func(reac_input['Reference Mass'],reac_input['Reference Temperature'],float(j),j_2) ###??? should this (np.sum(r_param)**2) be here? This puts it in dimensional form!
						
			if k_2 == 0:
				Dout.append(Constant(D[k,k_2]))
			if k_2 == 1:
				Din.append(Constant(D[k,k_2]))

	# Experimental Plot
#

	if reac_input['Experimental Data Folder'].lower() != 'none':
		
		for k,j in enumerate(legend_label):
			if j != 'timing': #and j != 'conVtime_1' and j != 'conVtime_2' and j != 'conVtime_3' and j != 'conVtime_0' 
				dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/flux_data/'+legend_label[k]+'.csv',header=None)
				if reac_input['Display Experimental Data'].lower() == 'true':
					ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
			else:
				pass

		if reac_input['Display Objective Points'].lower() == 'true':
			for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
				print(legend_label[k_fitting])
				if objSpecies[k_fitting] == '1':
					ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting])
#	

	# Inert Plot
	if reac_input['Infinite Inert'].lower() == 'true':
		# For reactant / product species
		analyticalTiming = np.arange(0, reac_input['Time Steps']*dt, dt).tolist()
		for kjc in range(0,monitored_gas):
			outlet = []
			
			if reac_input['Scale Output'].lower() == 'true':
				factor = 1
			else:
				factor = float(species_pulse_list[kjc])*reac_input['Reference Pulse Size']
			
			for k in range(0,int(reac_input['Time Steps'])+1):
				analyticalValue = 0
				if k*dt - species_time[kjc] > 0:
					for n in range(0,50):
						analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[kjc])*(D[kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
				else: 
					analyticalValue = -1
				if analyticalValue < 0:
					outlet.append(0)
				else:
					outlet.append(factor*(D[kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
			
			# Store analytical solution data
			np.savetxt('./analyticalCO19000.csv', outlet, delimiter=",")
			try:
				ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
			except ValueError:
				outlet = outlet[:-1]
				ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
		
		# For inert species
		for kjc in range(0,int(reac_input['Number of Inerts'])):
			outlet = []
			outlet.append(0)
		
			zAnalytical = 1
			while zAnalytical*dt < species_time[monitored_gas+kjc]:
				outlet.append(0)
				zAnalytical+=1	
			outlet.append(0)
			if reac_input['Scale Output'].lower() == 'true':
				factor = 1
			else:
				factor = float(species_pulse_list[monitored_gas+kjc])*reac_input['Reference Pulse Size']
			for k in range(zAnalytical+1,int(reac_input['Time Steps'])+1):
				analyticalValue = 0
				for n in range(0,50):
					analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[monitored_gas+kjc])*(D[monitored_gas+kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
				outlet.append(factor*(D[monitored_gas+kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
			try:
				ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical '+legend_label[monitored_gas+kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
			except ValueError:
				outlet = outlet[:-1]
				ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical '+legend_label[monitored_gas+kjc], alpha=0.7)
				#ax2.plot(graph_data['timing'],outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)

	for knum,k in enumerate(necessary_values['reactants']):
		if knum < necessary_values['molecules_in_gas_phase']:

			dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/'+k+'.csv',header=None)
			graphLimit = dfTemp[0][len(dfTemp[0]) - 1]
			
			if type(pulse) == int:
				ax2.plot(dfTemp[0],dfTemp[pulse],color=colors[knum], ls = '--',label=legend_label[knum], alpha=0.7)
			
			elif type(pulse) == list:
				legendInclude = True
				for j in pulse:
					if legendInclude == True:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label=legend_label[knum], alpha=0.7)
						legendInclude = False
					else:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label='_nolegend_', alpha=0.7)
			else:
				legendInclude = True
				for j in range(1, len(dfTemp.columns)):
					if legendInclude == True:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label=legend_label[knum], alpha=0.7)
						legendInclude = False
					else:
						ax2.plot(dfTemp[0],dfTemp[j],color=colors[knum], ls = '--', label='_nolegend_', alpha=0.7)
	ax2.set_xlim(0,graphLimit)
	for k in range(0,int(reac_input['Number of Inerts'])):
		
		dfTemp = pd.read_csv(reac_input['Output Folder Name']+'_folder/flux_data/Inert-'+str(k+1)+'.csv',header=None)
		if type(pulse) == int:
			ax2.plot(dfTemp[0],dfTemp[pulse],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label=legend_label[necessary_values['molecules_in_gas_phase']+k], alpha=0.7)
			
		elif type(pulse) == list:
			legendInclude = True
			for j in pulse:
				if legendInclude == True:
					#print(legend_label[necessary_values['molecules_in_gas_phase']+k-1])
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label=legend_label[necessary_values['molecules_in_gas_phase']+k], alpha=0.7)
					legendInclude = False
				else:
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label='_nolegend_', alpha=0.7)

		else:
			legendInclude = True
			for j in range(1, len(dfTemp.columns)):
				if legendInclude == True:
					#print(legend_label[necessary_values['molecules_in_gas_phase']+k-1])
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label=legend_label[necessary_values['molecules_in_gas_phase']+k], alpha=0.7)
					legendInclude = False
				else:
					ax2.plot(dfTemp[0],dfTemp[j],color=colors[necessary_values['molecules_in_gas_phase']+k], ls = '--', label='_nolegend_', alpha=0.7)

	ax2.legend(title="Gas Species",loc = 1)

	if store_graph == True:
		plt.savefig(output_name)
	
if show_graph == True:
		plt.show()

			if reac_input['Infinite Inert'].lower() == 'true':
				if k_pulse == 0:
	
					reac_time['Time Steps'] = 3000
					analyticalTiming = np.arange(0, reac_input['Time Steps']*dt, dt).tolist()
					for kjc in range(0,monitored_gas):
						outlet = []
						
						if reac_input['Scale Output'].lower() == 'true':
							factor = 1
						else:
							factor = float(species_pulse_list[kjc])*reac_input['Reference Pulse Size']
						
						for k in range(0,int(reac_input['Time Steps'])+1):
							analyticalValue = 0
							if k*dt - species_time[kjc] > 0:
								for n in range(0,50):
									analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[kjc])*(D[kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
							else: 
								analyticalValue = -1
							if analyticalValue < 0:
								outlet.append(0)
							else:
								outlet.append(factor*(D[kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
						
						np.savetxt('./analyticalCO19000.csv', outlet, delimiter=",")
	
						try:
							ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)

						except ValueError:
							outlet = outlet[:-1]
							ax2.plot(analyticalTiming,outlet,color=colors[kjc],label='Analytical '+legend_label[kjc], alpha=0.7)
					
					for kjc in range(0,int(reac_input['Number of Inerts'])):
						outlet = []
						outlet.append(0)
					
						zAnalytical = 1
	
						while zAnalytical*dt < species_time[monitored_gas+kjc]:
							outlet.append(0)
							zAnalytical+=1	
						outlet.append(0)
	
						if reac_input['Scale Output'].lower() == 'true':
							factor = 1
						else:
							factor = float(species_pulse_list[monitored_gas+kjc])*reac_input['Reference Pulse Size']
	
						for k in range(zAnalytical+1,int(reac_input['Time Steps'])+1):
							analyticalValue = 0
							for n in range(0,50):
								analyticalValue += ((-1)**n)*(2*n+1)*np.exp((-(n+0.5)**2)*(3.14159**2)*((k*dt - species_time[monitored_gas+kjc])*(D[monitored_gas+kjc][0]/(eb[0]*(np.sum(r_param)**2)))))
							outlet.append(factor*(D[monitored_gas+kjc][0]*3.14159/(eb[0]*(np.sum(r_param)**2)))*analyticalValue)
	
						try:
							ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)
						except ValueError:
							outlet = outlet[:-1]
							ax2.plot(analyticalTiming,outlet,color=colors[monitored_gas+kjc],label='Analytical'+legend_label[monitored_gas+kjc], alpha=0.7)




			if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Display Experimental Data'].lower() == 'true':
	
				if reac_input['Display Objective Points'].lower() == 'true' or reac_input['Fit Parameters'].lower() == 'true':
					if k_pulse == 0:
						for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
							if objSpecies[k_fitting] == '1':
								ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting])
			
				if reac_input['Display Objective Points'].lower() == 'true':# and reac_input['Fit Inert'].lower() == 'true':
					for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
						if objSpecies[k_fitting] == '1':
							ax2.scatter(output_fitting[legend_label[k_fitting]]['times'],output_fitting[legend_label[k_fitting]]['values'],marker='^',color=colors[k_fitting],label='Fitting'+legend_label[k_fitting], alpha=0.3)
	
	
				for k,j in enumerate(graph_data):
					if j != 'timing':
	
						if k_pulse > 0:
							pass
						else:
							dfNew = pd.read_csv(reac_input['Experimental Data Folder']+'/flux_data/'+legend_label[k]+'.csv',header=None)

							ax2.scatter(dfNew[0][:],dfNew[1][:],color=colors[k],label='exp '+legend_label[k], alpha=0.2) #, ls = '--'
					else:
						pass
	