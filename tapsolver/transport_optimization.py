#transport_optimization


# Select inert gasses to include in the objective function
OCM_tapobject.gasses_objective = ['Ar']
# Select parameters to optimize
OCM_tapobject.parameters_of_interest = ["TAPobject_data.reactor_species.inert_gasses['Ar'].intensity","TAPobject_data.reactor_species.reference_mass"]
# Specify the optimization process
OCM_tapobject.optimize = True 
OCM_tapobject.objective_function = True

def min_function(x):
	print(x[0])
	print(x[1])
	# Redefine parameters based on current guest/iteration
	OCM_tapobject.reactor_species.reference_mass = x[0]
	OCM_tapobject.reactor_species.inert_gasses['Ar'].intensity = x[1]
	obj_new = forward_problem(0.5,1,OCM_tapobject)
	print(obj_new)
	return obj_new
# Initial guess for parameters
x0 = np.array([23,0.11])
res = minimize(min_function,x0,method='nelder-mead')
