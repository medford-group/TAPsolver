

			if k_pulse == 0:
				dictionary_of_numpy_data = {}
				for j_species in range(0,monitored_gas+int(reac_input['Number of Inerts'])):
					dictionary_of_numpy_data[legend_label[j_species]] = np.empty(shape=[0, len(graph_data['timing'])])
					new_data = np.asarray(graph_data['timing'])
					dictionary_of_numpy_data[legend_label[j_species]] = np.vstack((dictionary_of_numpy_data[legend_label[j_species]],new_data))
				
				for k_species in range(monitored_gas+1,len(necessary_values['reactants'])+1):
					dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.empty(shape=[0, len(graph_data['timing'])])
					new_data = np.asarray(graph_data['timing'])
					dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]] = np.vstack((dictionary_of_numpy_data[necessary_values['reactants'][k_species-1]],new_data))
	









