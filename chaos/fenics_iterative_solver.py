def solver_iteration(time_step,dk,solver,method):
	if method == 'simple_adaptive':
		try:
			dk.assign(time_step/1.1)
			u_temp.assign(u)
			solver.solve()
			u_1.assign(u)
			u.assign(u_temp)
			dk.assign(time_step*1.1)
			solver.solve()
			u_3.assign(u)
			u.assign(u_temp)
			dk.assign(time_step)
			solver.solve()
			u_2.assign(u)
	
			#u.assign(u_temp)
			#print('')
			#print('norms')
			##ref_norm = u_2.vector().norm("l2")
			##norm_1 = ((u_1.vector() - u_2.vector()))
			##norm1 = norm_1.norm("l2")/ref_norm
			#print(print(norm1))
			##norm_2 = (u_2.vector() - u_3.vector())
			##norm2 = norm_2.norm("l2")/ref_norm
			#print(print(norm2))
			#print('')
	
			if norm1 > 0.02:
				time_step = time_step/1.1
				dk.assign(time_step)
				u.assign(u_1)
				#print("")
				#print("decrease")
				#print("")
			elif norm2 < 0.02:
				time_step = time_step*1.1
				dk.assign(time_step)
				u.assign(u_3)
				#print("")
				#print("Increase")
			#print("")
			#time.sleep(0.4)
			return time_step
		
		except RuntimeError:
			time_step = time_step*0.5
			print(time_step)
			dk.assign(time_step)
			if time_step < 1e-5:
				print("Time step too low")
				print(time.time() - start_time)
				sys.exit()
			time_step=solver_iteration(time_step)
			return time_step
	elif method == 'None':
		try:
			solver.solve()
			return time_step
			
		except RuntimeError:
			time_step = time_step*0.5
			print(time_step)
			dk.assign(time_step)
			if time_step < 1e-5:
				print("Time step too low")
				print(time.time() - start_time)
				sys.exit()
			time_step=solver_iteration(time_step)
			return time_step