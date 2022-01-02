def errorOutput(elem_reacs):
	
	"""Return this error in the event that the user doesn't define rate constants correct"""

	print("       ")
	print("       ")
	print("NAME ERROR")
	print("       ")
	print("       ")
	print("Likely Need to enter the rate constants.")
	print("There is currently no way to generate or ")
	print("gather the rate constants automatically")
	print("       ")
	print("Must follow the following format")
	print("       ")
	print("'kf1' = forward rate constant for reaction # 1")
	print("'kb1' = reverse rate constant for reaction # 1")
	print("       ")
	print("Rate constants for the following equations")
	print("must be included")
	print("       ")
	for j_nu,k_nu in enumerate(elem_reacs):
		print("Reaction "+str(j_nu+1)+"     "+k_nu)
	print("       ")
	print("Could also have incorrect diffusion coefficient array defined")
	sys.exit()
