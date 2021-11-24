
import numpy as np

def flux_graph(TAPobject_data: TAPobject,data_sets):

	time_steps = 3*1000
	dt = 3/time_steps
	
	fig2, ax2 = plt.subplots()

	ax2.set_xlabel('$Time\ (s)$', fontsize = 14)


	colors = ['b','orange','g','r','k','y','c','m','brown','darkgreen','goldenrod','lavender','lime']

	if TAPobject_data.scaled_graph == True:
		ax2.set_ylabel('$Outlet\ Flow\ (1/s)$', fontsize = 14)
	else:
		ax2.set_ylabel('$Outlet\ Flow\ (nmol/s)$', fontsize = 14)
	
	legend_label = []
	
	if reac_input['Experimental Data Folder'].lower() != 'none':
		output_fitting = curveFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
		
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

	# Inert Plot
	if TAPobject_data.display_analytical == True:
		# For reactant / product species

		analyticalTiming = np.arange(0, time_steps*dt, dt).tolist()
		for kjc in TAPobject_data.reactor_species.gasses:
			outlet = []
		
			if reac_input['Scale Output'].lower() == 'true':
				factor = 1
			else:
				factor = TAPobject_data.reactor_species.gasses[kjc].intensity*TAPobject_data.reactor_species.reference_pulse_size
			
			for k in range(0,time_steps+1):
				analyticalValue = 0
				if k*dt - TAPobject_data.reactor_species.gasses[kjc].delay > 0:
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
	