def solverIteration(time_step,method,solver,dk,dec_tim,inc_tim):

	try:
		if method == 'simple_adaptive':
	
			u_temp.assign(u)
			uout_1 = call_solver(dk.assign(time_step/dec_tim),u_temp,u,u_1,solver)
			uout_3 = call_solver(dk.assign(time_step*inc_tim),u_temp,u,u_3,solver)
			uout_2 = call_solver(dk.assign(time_step),u_temp,u,u_2,solver,keep_sol=True)
			time_step = norm_comp(uout_1,uout_2,uout_3,dec_tim,inc_tim)

			return time_step
		
		elif method == 'None':
			solver.solve(annotate = False) # ### can pass annotate = False if I don't want it to record the solution
			return time_step
	except RuntimeError:
		print('Time Step Failure')
		sys.exit()


def normComp(u1,u2,u3,d_t,i_t):
	ref_norm = u2.vector().norm("l2")
	norm1 = u1.vector().norm('l2')/ref_norm
	norm2 = u3.vector().norm('l2')/ref_norm

	if norm1 > 0.02:
		time_step = time_step/d_t
		dk.assign(time_step)
		u.assign(u1)
	elif norm2 < 0.01:
		time_step = time_step*i_t
		dk.assign(time_step)
		u.assign(u3)
	return time_step