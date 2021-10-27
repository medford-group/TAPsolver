		if reac_input['Sensitivity Analysis'].lower() == 'true':
	
			if reac_input['Uncertainty Quantification'].lower() == 'true' or (reac_input['Sensitivity Analysis'].lower() == 'true' and sens_type == 'trans'):
				if reac_input['Sensitivity Parameter'].find('Ga') > -1:
					c = r_Ga_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = r_Ga_in[reac_input['Sensitivity Parameter']]
				elif reac_input['Sensitivity Parameter'].find('dG') > -1:
					c = r_dG_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = r_dG_in[reac_input['Sensitivity Parameter']]
				elif reac_input['Sensitivity Parameter'].find('Ao') > -1:
					c = Ao_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = Ao_in[reac_input['Sensitivity Parameter']]
				elif reac_input['Sensitivity Parameter'].find('Ea') > -1:
					c = Ea_in[reac_input['Sensitivity Parameter']]
					c.tlm_value = Ea_in[reac_input['Sensitivity Parameter']]
				else:
					c = r_const[reac_input['Sensitivity Parameter']]
					c.tlm_value = r_const[reac_input['Sensitivity Parameter']]
	
				SV_du = FunctionSpace(mesh,P1)
				Sw_new = Expression('A',A=Constant(1),degree=0)
				Sw_new2 = interpolate(Sw_new,SV_du)
				Sw3 = project(Sw_new2,SV_du)
			
			elif sens_type == 'total' or sens_type == 'initial':
				pass
	
			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit()
	
		if reac_input['Sensitivity Analysis'].lower() == 'true':
			
			if sens_type == 'trans':
				sensFuncs = {}
				sensRates = {}
				for k_gasses in range(0,len(necessary_values['reactants'])):

					sensFuncs[str(k_gasses)] = []
					sensRates[str(k_gasses)] = []

			elif sens_type == 'total' or sens_type == 'initial':
				pass
			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit()
	
		sensAll = []
		sensAll2 = []
		
		sens_time_list = []
		
		if reac_input['Sensitivity Analysis'].lower() == 'true':
			if sens_type == 'trans':
				path_2 = reac_input['Output Folder Name']+'_folder/sensitivity/'
				generateFolder(path_2)
			
				path_molecules = path_2+reac_input['Sensitivity Parameter']
				generateFolder(path_molecules)
				sensFolder = path_molecules
	
			elif sens_type == 'total' or sens_type == 'initial':
				path_2 = reac_input['Output Folder Name']+'_folder/sensitivity/'
				generateFolder(path_2)			
	
			else:
				print('Sensitivity analysis is not properly defined.')
				sys.exit() 



			if sens_type == 'trans':
				sensitivity_output = {}
				RRM_der = {}
			
				for k_sens in range(monitored_gas):
					sensitivity_output[k_sens] = []
	
				for k_sens in range(len(necessary_values['reactants'])):
					RRM_der[k_sens] = []
			
			elif sens_type == 'total' or sens_type != 'initial':
				pass
	
			else:
				#print('Sensitivity analysis is not properly defined.')
				pass
				#sys.exit()


		sens_type = sensitivityType
			

				#############################################################
				################ STORE SENSITIVITY DATA #####################
				#############################################################

				if reac_input['Sensitivity Analysis'].lower() == 'true':
					
					if sens_type == 'trans':
	
						u_final = u.split(deepcopy=False)
						should_it = int(round(t*reac_input['Time Steps']/reac_input['Pulse Duration'],0))
					
						for k in range(0,monitored_gas):
							new_val = (( to_flux[k]*u.vector().get_local()[(all_molecules)+k]))
							u_graph_data['conVtime_'+str(k)].append((new_val))
						
						testrrm = rrmEqs(rate_array,rev_irr,'dx(1)',gForward,arrForward,arrBackward)

						for kGasses in range(0,len(necessary_values['reactants'])):
							#print(assemble( ( inner(u[kGasses], Sw3/Constant(0.0075000000000000015)) )* dx(1)))
							sensFuncs[str(kGasses)].append(assemble( ( inner(u[kGasses], Sw3/Constant(0.0075000000000000015)) )* dx(1)))
							
						for kGasses in range(0,len(necessary_values['reactants'])):
							#print(eval(testrrm[kGasses]))
							sensRates[str(kGasses)].append(eval(testrrm[kGasses]))


						#sys.exit()

					elif sens_type == 'total' or sens_type == 'initial':
	
						pass
	
					else:
						print('Sensitivity analysis is not properly defined.')
						sys.exit()
				
				
				sens_time_list.append(t)
				progressBar(t, reac_input['Pulse Duration'])
				#time.sleep(0.2)
				u_n.assign(u)
				
				t += dt
				constantT.assign(round(t,6))
	
			print()
			print(processTime(start_time))



			def derivCB_sens(j,dj,m):

				djv = [v.values()[0] for v in dj]

				with open('./derivSens.txt', 'w') as f:
					f.write(str(djv))
					f.close
				with open('./obj.txt', 'w') as f:
					f.write(str(j))
					f.close

			if reac_input['Sensitivity Analysis'].lower() == 'true':
				if sens_type == 'trans':
					print()
					start_time = time.time()
					print('Evaluating Tape with Tangent Linear Method. Could take some time.')
					tape2.evaluate_tlm()
	
				elif sens_type == 'total' or sens_type == 'initial':
					rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB_sens)
					rf_2.derivative()
#					dJdm = compute_gradient(jfunc_2, controls)
#					djv = [v.values()[0] for v in dJdm]
#					print(djv)
#					with open('./'+reac_input['Output Folder Name']+'_folder/sensitivity/objSens.txt', 'w') as f:
#						f.write("Parameters: "+str(legend_2))
#						f.write("Change: "+str(djv))
#						f.close
	
				else:
					print('Sensitivity analysis is not properly defined.')
					sys.exit()
	
	
			if reac_input['Sensitivity Analysis'].lower() == 'true':
	
	
				if sens_type == 'trans':
					for numEachSens,eachSens in enumerate(u_graph_data):
						np.savetxt(sensFolder+'/c_'+legend_label[numEachSens]+'.csv',u_graph_data[eachSens],delimiter=",")
	
					for numEachSens,eachSens in enumerate(sensFuncs):
						
						newList = []
						for kSensNum, kSens in enumerate(sensFuncs[eachSens]):
							newValue = kSens.block_variable.tlm_value
							newList.append(newValue)
						np.savetxt(sensFolder+'/dc_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
						
					for numEachSens,eachSens in enumerate(sensRates):
						newList = []
						for kSensNum, kSens in enumerate(sensRates[eachSens]):
							newValue = kSens.block_variable.tlm_value
							newList.append(newValue)
							np.savetxt(sensFolder+'/dr_'+necessary_values['reactants'][numEachSens]+'.csv',newList,delimiter=",")
		

				elif sens_type == 'total' or sens_type == 'initial':
					pass
	
				else:
					print('Sensitivity analysis is not properly defined.')
					sys.exit()
		
			if reac_input['Sensitivity Analysis'].lower() == 'true':
				print(processTime(start_time))
