			x_values = []
			it_times = []
			j_values = []
			dj_values = []
	
			#############################################################
			############# OPTIMIZATION ITERATION CALLBACKS ##############
			#############################################################
	
			def derivCB(j,dj,m):
				it_times.append(time.time())
				j_values.append(j)
				djv = [v.values()[0] for v in dj]
				dj_values.append(djv)
				mv = [v.values()[0] for v in m]
				x_values.append(mv)
				print("Derivative Time: "+str(time.time() - start_time))
				with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt', 'w') as f:
					f.write("Objective Value: "+str(j_values))
					f.write('\n')
					f.write("Change: "+str(dj_values))
					f.write('\n')
					f.write("Constants: "+str(x_values))

					f.close
				print(j)
				print(djv)
				print(mv)

	
			def hessCB(j,dj,m):
				it_times.append(time.time())
				j_values.append(j)
				djv = [v.values()[0] for v in dj]
				dj_values.append(djv)
				mv = [v.values()[0] for v in m]
				x_values.append(mv)
				print(time.time() - start_time)
				with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIterHess.txt', 'w') as f:
					f.write("Objective Value: "+str(j_values))
					f.write('\n')
					f.write("Change: "+str(dj_values))
					f.write('\n')
					f.write("Constants: "+str(x_values))
	
					f.close
				print(j)
				print(djv)
				print(mv)