		if reac_input['Uncertainty Quantification'].lower() == 'true':
			path_2 = reac_input['Output Folder Name']+'_folder/UQ/'
			generateFolder(path_2)
			hessFolder = path_2
			

		tape2 = Tape()
		tape2.clear_tape()
		set_working_tape(tape2)

			if reac_input['Uncertainty Quantification'].lower() == 'true':
				start_time = time.time()
				print()
				print('Calculating hessian. Could take some time.')

				rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2)# ,hessian_cb_post=hessCB

				rf_2.derivative()
				utest = []
				B = []
				for just in range(0,len(controls)):
					utest.append(Constant(0))

				for jay_z_num in range(0,len(controls)):
					utest = []
					for just in range(0,len(controls)):
						utest.append(Constant(0))
					utest[jay_z_num] = Constant(1)
					H_i = rf_2.hessian(utest)
					djv = [v.values()[0] for v in H_i]
					#print(djv)
					B.append(djv)

				hessian_array = np.array(B)

				B = hessian_array

				print('Finished generating hessian, storing now.')
				np.savetxt(hessFolder+'/hessian.csv', hessian_array, delimiter=",")
				try:
					print('The eigenvalues of the hessian are:')
					hessEigs = np.linalg.eig(hessian_array)[0]
					print(hessEigs)
					#eigenInfo = np.any((a < 0))
					#if eigenInfo == True:
					#	print('Not all eigenvalues are positive. If fitting parameters, might want to run longer.')
					np.savetxt(hessFolder+'/eigenvalues.csv', hessEigs, delimiter=",")
				except:
					print('Failed to determine eigenvalues')
				try:
					print('Generating Covariance Matrix by Inverting the Hessian')
					print(B)
					vx_new = np.linalg.inv(B)
					np.savetxt(hessFolder+'/covariance.csv', vx_new, delimiter=",")

				except:
					print('Failed to invert hessian')
				try:
					print('The first and second standard deviations of each parameter are:')
					#print('vx_new value')
					#print(vx_new)
					std_1 = np.diagonal(np.sqrt(vx_new))
					print(std_1)
					np.savetxt(hessFolder+'/std_1.csv', std_1, delimiter=",")
					std_2 = np.diagonal(2*np.sqrt(vx_new))
					print(std_2)
					np.savetxt(hessFolder+'/std_2.csv', std_2, delimiter=",")
				except:
					print('Failed to calculate confidence interval')
