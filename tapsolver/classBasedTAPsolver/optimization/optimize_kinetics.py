		# Define controls for optimization or sensitivity analysis
		if reac_input['Fit Inert'].lower() == 'true':
			controls = []
			legend_2 = []
			for j in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])]),len(legend_label)):
				print(j)
				controls.append(Control(Dout[j]))
				legend_2.append(j)
	
		if str(reac_input['Objective Species']).find(',') != -1:
			objSpecies = list(reac_input['Objective Species'].split(','))
		else:
			print(reac_input['Objective Species'])
			objSpecies = [str(int(reac_input['Objective Species']))]
	
if reac_input['Experimental Data Folder'].lower() != 'none' and reac_input['Fit Inert'].lower() != 'true' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or sampling == True or ((sens_type == 'total' or sens_type == 'initial') and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True or reac_input['Fitting Gif'].lower() == True:
			
			if reac_input['Uncertainty Quantification'].lower() == 'true':
				print("Uncertainty Quantification")
			try:
				if type(reac_input['Objective Points']) == float:
					output_fitting = pointFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
				elif reac_input['Objective Points'] == 'all':
					output_fitting = curveFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)
									
				else:
					print('Objective Points defined incorrectly')
					sys.exit()
			except TypeError:
				print('Objective Point Input Is Not Valid')
				sys.exit()
	
		if reac_input['Fit Inert'].lower() == 'true':# or reac_input['Display Objective Points'].lower() == 'true':
			try:
				if type(reac_input['Objective Points']) == int:
					output_fitting = pointFitting(legend_label[int(reac_input['Number of Inerts'])+1:],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'])
				elif reac_input['Objective Points'] == 'all':
					
					output_fitting = curveFitting(legend_label[(len(legend_label) - int(reac_input['Number of Inerts']) ) :],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies[(len(legend_label) - int(reac_input['Number of Inerts']) ) :])
					
				else:
					print('Objective Points defined incorrectly')
					sys.exit()
			except TypeError:
				print('Objective Point Input Is Not Valid')
				sys.exit()


		if reac_input['Experimental Data Folder'].lower() != 'none' and (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Display Objective Points'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true') or sampling == True or ((sens_type == 'total' or sens_type == 'initial') and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True:
			for k_fitting in range(0,len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])])):
				if objSpecies[k_fitting] == '1':
					for timeStep in range(0,len(output_fitting[legend_label[k_fitting]]['times'])):
						output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)
	
		if reac_input['Fit Inert'].lower() == 'true':
			for k_fitting in range(len(legend_label[:int(len(legend_label)-reac_input['Number of Inerts']+1)]),len(legend_label)):
				for timeStep in range(0,len(output_fitting[legend_label[k_fitting]]['times'])):
					output_fitting[legend_label[k_fitting]]['times'][timeStep] = round(output_fitting[legend_label[k_fitting]]['times'][timeStep],6)




			if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or sampling == True or fit_temperature == True:
				osub = integration_section()
				domains = MeshFunction("size_t", mesh,0)
				domains.set_all(0)
				osub.mark(domains, 1)
				dP = Measure('vertex',domain = mesh, subdomain_data=domains)



							try:
								if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
									solvertemp.solve()
		
								else:

									solvertemp.solve(annotate = False)
		
							except RuntimeError:
								print('this one')
								print('Time Step Failure')
								sys.exit()




			fitting_time = time.time()
			#print("Including temporary work around. Must remove line just below this!")
			# simplifiedTimeStep == False and 

			if (reac_input['Fit Parameters'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or fit_temperature == True):
					
				start_time = time.time()
				print()
				print('Fitting Kinetic Parameters. Will take some time!')

				######################## objective optimization (boukouvala)

				if reac_input['Optimization Method'] == 'objective':
					return float(jfunc_2)
					
				#######################

				if reac_input['Optimization Method'] == 'branch&bound':
					import PyDDSBB
					rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
					rf_2np = adReduNp.ReducedFunctionalNumPy(rf_2)
					
					def calc_loss(p):
						print('Iteration')
						print(p)
						estimate = rf_2np.__call__(p)
						return estimate
						#return rf_2np.__call__(np.array([0.5,17.892023742960912]))
					
					test_problem = PyDDSBB.DDSBBModel.Problem() ## Initialize the problem
					test_problem.add_objective(calc_loss, sense = 'minimize') ## Add objective function
					for bb in range(0,len(controls)):
						test_problem.add_variable(0,20)
					
					test = PyDDSBB.DDSBB(100,split_method = 'equal_bisection', variable_selection = 'longest_side', multifidelity = False, stop_option = {'absolute_tolerance': 100, 'relative_tolerance': 100, 'minimum_bound': 0.01, 'sampling_limit': 1000, 'time_limit': 400})
					test.optimize(test_problem)
					test.print_result()
				#######################

								######################## objective optimization (boukouvala)
				if reac_input['Optimization Method'] == 'genetic':
					from geneticalgorithm import geneticalgorithm as ga
					rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
					rf_2np = adReduNp.ReducedFunctionalNumPy(rf_2)
					
					def calc_loss(p):
						print('Iteration')
						print(p)
						estimate = rf_2np.__call__(p)
						print(estimate)
						return estimate
						#return rf_2np.__call__(np.array([0.5,17.892023742960912]))
					
					varbound = np.array([[0,10]]*len(controls))
					model=ga(function=calc_loss,dimension=len(controls),variable_type='real',function_timeout=800,variable_boundaries=varbound)
					model.run()
				#######################

				rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB)# ,hessian_cb_post=hessCB
					
				low_bounds = []
				up_bounds = []
	
				if reac_input['Fit Parameters'].lower() == 'true':
					for gt in range(0,len(controls)):
						low_bounds.append(0)
						up_bounds.append(np.inf)
				elif doe_form_pulse == True:
					for gt in range(0,len(controls)):
						low_bounds.append(0)
						up_bounds.append(5)
		
				if reac_input['Optimization Method'] == 'L-BFGS-B' or reac_input['Optimization Method'] == '':
					u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds), tol=1e-22, options={"ftol":1e-22,"gtol":1e-22})
					#u_opt_2 = minimize(rf_2, bounds = [low_bounds,up_bounds], tol=1e-9, options={"ftol":1e-9,"gtol":1e-9})
				elif reac_input['Optimization Method'] == 'Newton-CG':
					u_opt_2 = minimize(rf_2, method = 'Newton-CG',tol=1e-22, options={"xtol":1e-22})
				elif reac_input['Optimization Method'] == 'BFGS':
					u_opt_2 = minimize(rf_2, method = 'BFGS',tol=1e-22, options={"gtol":1e-22})# , "constraints":bounds
				elif reac_input['Optimization Method'] == 'SLSQP':
					u_opt_2 = minimize(rf_2, method = 'SLSQP', bounds = (low_bounds,up_bounds),tol=1e-22, options={"ftol":1e-22})
				elif reac_input['Optimization Method'] == 'CG':
					u_opt_2 = minimize(rf_2,bounds = (low_bounds,up_bounds), method = 'CG',tol=1e-22, options={"gtol":1e-22})
				elif reac_input['Optimization Method'] == 'basinhopping':
					u_opt_2 = minimize(rf_2, method = 'basinhopping', bounds = (low_bounds,up_bounds),tol=1e-22, options={"ftol":1e-22,"gtol":1e-22})
				elif reac_input['Optimization Method'] == 'nonlinear':
					print('non-linear optimization')
					print(low_bounds)
					print(up_bounds)
					problem = MinimizationProblem(rf_2) # ,bounds = (low_bounds,up_bounds)
					ipoptSolver = IPOPTSolver(problem)
					rf_2 = ipoptSolver.solve()
				else:
					print('Requested Optimization Method Does Not Exist')
					sys.exit()
				#sys.exit()
				print(processTime(start_time))
	
				optimization_success = True

