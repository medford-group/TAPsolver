def generateGif(molecules,exp_loc,fit_loc,all_steps,constants,reactions,time_data,xscale,yscale):
	
	"""
	Return a gif showing changes made during the optimization process
	"""

	def add_subplot_axes(ax,rect,axisbg='w'):
		
		"""
		Generates the subplot to help visualize some other quantitiy
		"""

		fig = plt.gcf()
		box = ax.get_position()
		width = box.width
		height = box.height
		inax_position  = ax.transAxes.transform(rect[0:2])
		transFigure = fig.transFigure.inverted()
		infig_position = transFigure.transform(inax_position)    
		
		x = infig_position[0]
		y = infig_position[1]
		width *= rect[2]
		height *= rect[3]
		subax = fig.add_axes([x,y,width,height])
		
		x_labelsize = subax.get_xticklabels()[0].get_size()
		y_labelsize = subax.get_yticklabels()[0].get_size()
		x_labelsize *= rect[2]**0.5
		y_labelsize *= rect[3]**0.5
		
		subax.xaxis.set_tick_params(labelsize=x_labelsize)
		subax.yaxis.set_tick_params(labelsize=y_labelsize)
		
		return subax
	
	x_data = list(range(0, all_steps))
	y_data = [0]
	
	for k in range(0,all_steps-1):
		y_data.append( (time_data[k+1] - time_data[k]) / 60 )

	def tap_plot(step):
		fig, ax = plt.subplots(figsize=(10,5))
	
		ax.grid()
		
		exp_data = {}
		sim_data = {}
		for k_names in molecules:

			exp_data[k_names] = pd.read_csv(exp_loc+'/'+k_names+'.csv',header=None)
			sim_data[k_names] = pd.read_csv(fit_loc+'/iter_'+str(step)+'_folder/flux_data/'+k_names+'.csv',header=None)
		
		for k_names in molecules:
			if yscale == 'normalized':
				ax.scatter(exp_data[k_names][0], exp_data[k_names][1]/(exp_data[k_names][1].max()),label="Exp. "+k_names,alpha=0.3)
			else:
				ax.scatter(exp_data[k_names][0], exp_data[k_names][1],label="Exp. "+k_names,alpha=0.3)
		for k_names in molecules:
			if yscale == 'normalized':
				ax.plot(sim_data[k_names][0], sim_data[k_names][1]/(exp_data[k_names][1].max()),label="Syn. "+k_names,ls='--')
			else:
				ax.plot(sim_data[k_names][0], sim_data[k_names][1],label="Syn. "+k_names,ls='--')

		ax.set_xlabel('Time (s)', fontsize=16)
		if yscale == 'normalized':
			ax.set_ylabel('Normalized Flow (1/s)', fontsize=16)
		else:
			ax.set_ylabel('Flow (nmol/s)', fontsize=16)

		props = dict(facecolor='white')

		ax.text(0.42, 0.95,'Iteration: '+str(step),transform=ax.transAxes, fontsize=14,verticalalignment='top',bbox=props)
		peak_peak = 0

		ax.legend(title='Gas Species',loc='upper right')
		if xscale == 'log':
			ax.set_xscale('log')
			ax.set_xlim(0.002,1)
		if yscale == 'normalized':
			ax.set_ylim(0,1.5)
		
		fig.canvas.draw()
		image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
		image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	
		return image
	
	kwargs_write = {'fps':4.0, 'quantizer':'nq'}
	components = list(range(0,all_steps))
	for zoo in range(0,10):
		components.append(all_steps-1) 

	imageio.mimsave(fit_loc+'/output.gif', [tap_plot(i) for i in components], fps=4)

		if reac_input['Fitting Gif'].lower() == 'true':
			print('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt')
			print(os.path.exists('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt'))

			if os.path.exists('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt') == True:
				with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt', 'r') as f:
					lines = f.readlines()
				f.close
				lines = [x.strip() for x in lines]
				times = lines[0] 
				times = times.replace('Objective Value: ','')
				times = eval(times)
				constants = lines[2]
				constants = constants.replace('Constants: ','')
				constants = eval(constants)
	
				things = len(times)
				print('things')
				for k_num in range(0,things):
					alter = pd.read_csv(input_file,header=None)
	
					variables_to_change = ['Display Graph','Fit Parameters','Sensitivity Analysis','Store Outlet Flux','Output Folder Name','Reaction_Information']
				
					for k in range(alter.shape[0]):
						if alter[0][k] in variables_to_change:
							if alter[0][k] == 'Store Outlet Flux':
								alter.iloc[k,1] = 'TRUE'
							elif alter[0][k] == 'Output Folder Name':
								alter.iloc[k,1] = reac_input['Output Folder Name']+'_folder/fitting/iter_'+str(k_num)
							elif alter[0][k] == 'Reaction_Information':
								value = 1
								current_rate = k+1
								kz = 0
								while kz < len(kVals):
									if value == 1:
										alter.iloc[current_rate,value] = constants[k_num][kz]
										
										value = 2
										kz+=1
									elif value == 2:
										if str(alter.iloc[current_rate,value]) != 'nan':
											alter.iloc[current_rate,value] = constants[k_num][kz]
											kz+=1
										else:
											pass
										current_rate += 1
										value = 1
	
							else:
								alter.iloc[k,1] = 'FALSE'
							#alter.iloc[51,1] = 'FALSE'
					
					alter.to_csv(input_file,header=None,index=False)	
				
					try:
						print()
						print('Iteration: '+str(k_num+1))
						call_sim()
					except:
						k_num = things
			
				generateGif(legend_label[:len(legend_label)], reac_input['Experimental Data Folder']+'/flux_data', './'+reac_input['Output Folder Name']+'_folder/fitting', len(constants), constants, reactor_kinetics_input['reactions_test'], times,reactor_kinetics_input['xscale'],reactor_kinetics_input['yscale'])
			
				for k_num in range(0,things):
					shutil.rmtree('./'+reac_input['Output Folder Name']+'_folder/fitting/iter_'+str(k_num)+'_folder') 
	
			user_data = pd.read_csv('./'+reac_input['Output Folder Name']+'_folder/'+input_file,header=None)
			user_data.to_csv(input_file,header=None,index=False)
	
	