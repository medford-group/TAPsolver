
											if thermoConstraints == True:
												w_temp = {}
												w_temp2 = {}
												w4 = {}
												tempFunc = {}

												for kip in thermoReactions.keys():
													w_temp[kip] = Expression(str(reactor_kinetics_input['thermo values'][kip]),degree=0) # deltaG = sum(-R*T*ln(kf/kb))
													w_temp2[kip] = interpolate(w_temp[kip],V_du)
													w4[kip] = project(w_temp2[kip],V_du)		
													thermoWeight = 1e-2	

													for jnum,jval in enumerate(thermoReactions[kip]):
																				
														if jnum == 0:
															tempFunc[kip] = thermoStoich[kip][jnum]*(-0.008314*reac_input['Reactor Temperature'])*ln(r_const["kf"+str(jval-1)]/r_const["kb"+str(jval-1)])
														else:
															tempFunc[kip] += thermoStoich[kip][jnum]*(-0.008314*reac_input['Reactor Temperature'])*ln(r_const["kf"+str(jval-1)]/r_const["kb"+str(jval-1)])

													#tempFunc = (-w4)*thermoWeight
													#tempFunc = (tempFunc)*thermoWeight
													#tempFunc_2 = (w4_2)*thermoWeight
													tempFunc[kip] = (tempFunc[kip]-w4[kip])*thermoWeight
													jfunc_2 += assemble(inner(tempFunc[kip],tempFunc[kip])*dx())
													#jfunc_2 += assemble(inner(tempFunc,tempFunc)*dx())
													print('Simulated Thermo Value')
													print(jfunc_2)

													#sys.exit()f
										else:
											pass