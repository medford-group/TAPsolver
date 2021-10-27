			reac_input['Taylor Test'] = 'False'

			if reac_input['Taylor Test'].lower() == 'true':
				taylorTest = False
	
				if taylorTest == True:
					h = Constant(0.0001)  # the direction of the perturbation
					for jContNum, jCont in enumerate(legend_2): #legend_2
						print('taylor test for parameter '+jCont)
						tTest = ReducedFunctional(jfunc_2, [controls[jContNum]],tape=tape2,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
						conv_rate = taylor_test(tTest, [r_const[jCont]], h)
						conv_rate = taylor_test(tTest, [r_const[jCont]], h, dJdm = 0)
