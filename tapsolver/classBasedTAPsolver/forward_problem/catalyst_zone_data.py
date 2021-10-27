if reac_input['Thin-Zone Analysis'].lower() == 'true':
				print('Storing Catalyst Zone Gas and Surface Data')
				if int(reac_input['Number of Pulses']) == 1:
					cat_dataRate = {}
					#for j_species in range(0,all_molecules):
					for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
						timeBetween = time.time()
						new_values = [] 
						mol_values = []
		
						for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
							
							if z != 0:
								#value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
								new_values = np.vstack((mol_values,cat_data['rawData'][:,z*(all_molecules)+j_species]))
								mol_values = new_values
					
							else:					
								new_values = cat_data['rawData'][:,z] 
								mol_values = cat_data['rawData'][:,z]
						
						top = mp.ceil((Mesh122)*reac_input['Mesh Size'])+1
						bottom = int((Mesh122)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-1
	
						## Change so that only thin-zone data is being stored
						if thinSize == 'cat':
							print(mol_values[top-1:top,:])
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values[top-1:top,:]),axis=1)
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top-1:top,:]),axis=1))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)	
							
						elif thinSize == 'point':
							centerPoint = int((top + bottom)/2)
							cat_dataRate['convtime_'+str(j_species)] = np.transpose(mol_values[centerPoint,:])
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.transpose(mol_values[centerPoint,:]), delimiter=",")
				
						elif thinSize == 'average':
							cat_dataRate['convtime_'+str(j_species)] = np.mean(np.transpose(mol_values[top:bottom,:]),axis=1)
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.mean(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
	
						elif thinSize == 'all':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values),axis=1)
							tempSurface = pd.DataFrame(np.flip(np.transpose(mol_values[top:bottom,:])))
							tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)
							#np.savetxt('./'+reac_input['Output Folder Name']+'_folder/thin_data/'+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values),axis=1), delimiter=",")
						
					#for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
					if True == False:
						print('Storing Catalyst Zone Gas and Surface Rate Data')
						rateTest = eval(rateStrings[j_species])
						
						tempSurface = pd.DataFrame(rateTest)
						tempSurface.to_csv('./'+reac_input['Output Folder Name']+'_folder/thin_data/r_'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)						
						
				else:
					pulse_path = './'+reac_input['Output Folder Name']+'_folder/thin_data/pulse_'+str(k_pulse+1)+'/'
					generateFolder(pulse_path)
					cat_dataRate = {}
					
					for j_species in range(0,all_molecules-int(reac_input['Number of Inerts'])):
						new_values = [] 
						mol_values = []
		
						for z in range(0, int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells):
	
							if z != 0:
								#value_cat = u_n.vector().get_local()[z*(all_molecules)+jk]
								new_values = np.vstack((mol_values,cat_data['rawData'][:,z*(all_molecules)+j_species]))
								mol_values = new_values
					
							else:					
								new_values = cat_data['rawData'][:,z] 
								mol_values = cat_data['rawData'][:,z]
						
						top = mp.ceil((Mesh122)*reac_input['Mesh Size'])+1
						bottom = int((Mesh122)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-1				
	
						## Change so that only thin-zone data is being stored
						if thinSize == 'cat':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values[top:bottom,:]),axis=1)
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
				
						elif thinSize == 'point':
							centerPoint = int((top + bottom)/2)
							cat_dataRate['convtime_'+str(j_species)] = np.transpose(mol_values[centerPoint,:])
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.transpose(mol_values[centerPoint,:]), delimiter=",")
				
						elif thinSize == 'average':
							cat_dataRate['convtime_'+str(j_species)] = np.mean(np.transpose(mol_values[top:bottom,:]),axis=1)
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.mean(np.transpose(mol_values[top:bottom,:]),axis=1), delimiter=",")
	
						elif thinSize == 'all':
							cat_dataRate['convtime_'+str(j_species)] = np.flip(np.transpose(mol_values),axis=1)
							np.savetxt(pulse_path+necessary_values['reactants'][j_species]+'.csv', np.flip(np.transpose(mol_values),axis=1), delimiter=",")
	
					if True == False:	
						rateTest = eval(rateStrings[j_species])
						tempSurface = pd.DataFrame(rateTest)
						tempSurface.to_csv(pulse_path+'r_'+necessary_values['reactants'][j_species]+'.csv',header=None,index=False)



