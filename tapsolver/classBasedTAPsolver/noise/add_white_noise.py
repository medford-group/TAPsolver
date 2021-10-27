
			sig = 0.1
			beta_2 = 0.00270
			w_2 = 2*3.14159*70
	
			for k,j in enumerate(graph_data):
				if j != 'timing':
					if reac_input['Noise'].lower() == 'true':
						for z in range(0,int(reac_input['Time Steps'])):
							graph_data[j][z] += np.random.normal(0,1)*sig +beta_2*np.cos(w_2*(k*dt))
					if k_pulse > 0:
						ax2.plot(graph_data['timing'],graph_data[j],color=colors[k], ls = '--', alpha=0.7)
					else:
						ax2.plot(graph_data['timing'],graph_data[j],color=colors[k],label=legend_label[k], ls = '--', alpha=0.7)
					new_data = np.asarray(graph_data[j])
					dictionary_of_numpy_data[legend_label[k]] = np.vstack((dictionary_of_numpy_data[legend_label[k]],new_data))
				else:
					pass	
	
			name_list = necessary_values['reactants'][monitored_gas:]