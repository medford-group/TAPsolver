
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

class TAPobject():

	"""
	
	This class acts as a container for the three necessary components to run a TAP simulation: the reactor details, mechanism and initial conditions.
	
	Args:
		
		reactor (reactor object): The reactor specifications.

		mechanism (mechanism object): The elementary processes and their associated kinetics.

		initial_conditions (initial_conditions object): The initial conditions of the TAP reactor/experiment.

	"""


	def __init__(self):
		
		self.reactor = None
		self.mechanism = None
		self.reactor_species = None
		self.experimental_data = None

		self.mesh = 200
		self.catalyst_mesh_density = 4
		self.output_name = 'exp_new'
		
		#self.point_volume = 0
	
	#@property
	#def mesh_step_size(self):
	#	return self.total_length()/mesh_size
	#def point_volume(self):
	#	return self.mesh_step_size()*cross_sectional_area()*self.zone_voids['zone0']
	#	# Initialize the grid system, time step size, pulse size
	#def fenics_cat_location(self):	
	#	return 1 - self.reactorCenterFraction()
		




#def establishMesh(in_1,cat,in_2,mesh_size):
#
#	"""Generate the FEniCS Mesh"""
#
#	r_param = np.array((in_1,cat,in_2))
#	
#	def last_grid_point(x):
#	    return x*mesh_size
#
#	last_func = np.vectorize(last_grid_point)
#	zone_sizes = last_func(r_param)
#	grid_loc_1 = np.round(np.cumsum(zone_sizes))
#	dx_r = np.sum(r_param)/mesh_size
#	dx2_r = dx_r*dx_r
#
#	frac_length = r_param[1]/(np.sum(r_param))
#
#	cat_location = (r_param[1]*0.5+r_param[0])/(np.sum(r_param))
#	return r_param,dx_r,dx2_r,frac_length,cat_location
#
#r_param, dx_r, dx2_r, frac_length, cat_location = establishMesh(reac_input['len_inert_1'],reac_input['len_cat'],reac_input['len_inert_2'],reac_input['Mesh Size'])
#		
#dk = Constant(reac_input['Pulse Duration']/reac_input['Time Steps'])


#def defineBCs(reactor,elem_list,nec_values,V_sec,reacs_num,all_mol,reac_ratio,L_bound,R_bound,number_of_inerts):
#
#	"""Define the appropriate boundary conditions for all monitored species"""
#
#	if reactor == 'tap':
#		bcs = []
#
#		if elem_list != ['INERT_ONLY']:
#
#			for k in range(0,nec_values):
#				bcs.append(DirichletBC(V_sec.sub(k),Constant(0),R_bound))
#			
#			for k in range(0,int(number_of_inerts)):
#				bcs.append(DirichletBC(V_sec.sub(all_mol-(1+k)),Constant(0),R_bound))
#			
#		else:	
#			bcs.append(DirichletBC(V_sec,Constant(0),R_bound))
#	
#	# Vacuum BC at outlet of reactor (C = 0 @ L = L_reactor)
#	elif reactor == 't_pfr' or 't_pfr_diff':
#		bcs = []
#		newValues = reac_ratio.split(',')
#		for k in range(0,len(newValues)):
#			bcs.append(DirichletBC(V_sec.sub(k),Constant(int(newValues[k])),L_bound))
#	
#	return bcs
#
#
#		#############################################################
#		######## INITIALIZATION OF FINITE ELEMENTS AND MESH #########
#		#############################################################
#	
#		# Define the base mesh (the initial uniform mesh)
#		mesh = UnitIntervalMesh(int(reac_input['Mesh Size']))
#	
#		# Refine the mesh (depending on user specification)
#		catalystRefinement = int(reac_input['Catalyst Mesh Density'])
#		cfDict = {}
#	
#		roundedMesh2 = ((1-cat_location) + 0.5*frac_length)*reac_input['Mesh Size']	
#		Mesh2 = round(((1-cat_location) + 0.5*frac_length)*reac_input['Mesh Size'])
#		Mesh22 = round(((1-cat_location) + 0.5*frac_length)*reac_input['Mesh Size'])/reac_input['Mesh Size']
#		Mesh222 = round(((cat_location) + 0.5*frac_length)*reac_input['Mesh Size'])/reac_input['Mesh Size']
#
#		roundedMesh1 = ((1-cat_location) - 0.5*frac_length)*reac_input['Mesh Size']
#		Mesh1 = round(((1-cat_location) - 0.5*frac_length)*reac_input['Mesh Size'])
#		Mesh12 = round(((1-cat_location) - 0.5*frac_length)*reac_input['Mesh Size'])/reac_input['Mesh Size']
#		Mesh122 = round(((cat_location) - 0.5*frac_length)*reac_input['Mesh Size'])/reac_input['Mesh Size']
#
#		if Mesh2 != roundedMesh2 or Mesh1 != roundedMesh1:
#			print('Warning: Catalyst zone will be refined and rounded to the nearest whole mesh point!')
#			trueMesh = (roundedMesh2 - roundedMesh1)/reac_input['Mesh Size']
#			newMesh = (Mesh2 - Mesh1)/reac_input['Mesh Size']
#			print()
#			print('New Catalyst Fraction = '+str(newMesh))
#			print('Old Catalyst Fraction = '+str(trueMesh))
#			percentChange = abs(round(100*(trueMesh - newMesh)/trueMesh,2))
#			print('Change = '+str(percentChange)+'%')
#			print()
#			if percentChange > 4:
#				print('Consider refining the mesh to improve the accuracy of the simulation!')
#				sys.exit()
#			
#		for jayz in range(0,catalystRefinement+1):
#			class thin_zoneTest(SubDomain):
#				def inside(self, x, on_boundary):
#					return between(x[0], ((Mesh12), (Mesh22)))
#			
#			thin_zoneTest = thin_zoneTest()
#			cfDict[jayz] = MeshFunction("bool",mesh,1)
#			
#			thin_zoneTest.mark(cfDict[jayz],jayz)
#			mesh = refine(mesh,cfDict[jayz],True)
#
#		# Generate element space (for all observables)
#		P1 = FiniteElement('CG',mesh.ufl_cell(),1)
#	
#		if reac_input['reactions_test'] != ['INERT_ONLY']:
#			test_new = eval(necessary_values['element'])
#			element = MixedElement(test_new)
#			V = FunctionSpace(mesh,element)
#			V_du = FunctionSpace(mesh,P1)
#		else:
#			V = FunctionSpace(mesh,P1)
#			V_du = FunctionSpace(mesh,P1)
#
#	
#		all_molecules = necessary_values['gas_num']
#	
#		u = Function(V)
#		if runge_kutta_approach == True:
#			W = Function(V)
#			W.interpolate(Constant(0.0))
#		u_n = Function(V)
#		u_temp = Function(V)
#
#		
#		if reac_input['Advection'].lower() == 'true':
#			
#			W = VectorFunctionSpace(mesh, 'P', 1)
#			advTerm = Function(W)
#			advMulti = Constant(reac_input['Advection Value'])
#			advValue = Constant(1)
#			advTerm.vector()[:] = advValue
#			if reac_input['Fit Inert'].lower() == 'true':
#				controls = []
#				controls.append(Control(advMulti))
#			else:
#				pass
#			advTerm.vector()[totalNumCells-0] = 0
#
#meshCells = int((Mesh22)*reac_input['Mesh Size']) - mp.ceil((Mesh12)*reac_input['Mesh Size'])
#		frac_temp = (meshCells*2**(int(reac_input['Catalyst Mesh Density'])))/(int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells)
#		totalNumCells = int(reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells
#	
#		# Define subdomains (thin-zone / final cell)
#		class thin_zone(SubDomain):
#			def inside(self, x, on_boundary):
#				return between(x[0], ((Mesh12)-1/reac_input['Mesh Size'], (Mesh22)+1/reac_input['Mesh Size']))
#				
#		class singlePoint(SubDomain):
#			def inside(self,x,on_boundary):
#				return between(x[0], (((1-cat_location) - 1/totalNumCells), ((1-cat_location) + 1/totalNumCells)))
#	
#		thin_zone = thin_zone()
#		domains = MeshFunction("size_t", mesh,1)
#		thin_zone.mark(domains,1)
#	
#		newBoundaries = domains.mesh().coordinates().transpose().tolist()[0]
#		additionalCells = len(newBoundaries[int(reac_input['Mesh Size']+1):])
#	
#		centralPoint = singlePoint()
#		boundary_parts0 = MeshFunction("size_t", mesh,0)
#		boundary_parts0.set_all(0)
#		centralPoint.mark(boundary_parts0, 1)
#	
#		dx = Measure("dx",subdomain_data=domains)
#		dT = Measure("dx",subdomain_data=boundary_parts0)
#		
#		# Define inlet and outlet boundaries
#		def boundary_L(x, on_boundary):
#			return on_boundary and near(x[0],0,tol)
#		
#		def boundary_R(x, on_boundary):
#			return on_boundary and near(x[0],1,tol)
#		
#		dz = 1/reac_input['Mesh Size']
#	
#		class integration_section(SubDomain):
#			def inside(self, x, on_boundary):
#				return between(x[0], (1-dz,1.0))
#		
#		right = CompiledSubDomain("near(x[0], 1.)")
#	
#		boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
#		right.mark(boundary_parts, 1)
#	
#		monitored_gas = necessary_values['molecules_in_gas_phase']
#		
#		bcs = defineBCs(reac_input['Reactor Type'],reac_input['reactions_test'],necessary_values['molecules_in_gas_phase'],V,reac_input['Number of Reactants'],all_molecules,reac_input['Pulse Size'],boundary_L,boundary_R,reac_input['Number of Inerts'])
###############################
#		constantT = Constant(0)
#		if doe_form_pulse == True:
#			#controls = []
#			b0Test2 = Expression('x[0] < 0.002500001 ? 0.5 : 0', degree=0)
#			
#			reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))
#
#			intensityFunctions ={}
#			intensConst = {}
#			for jk in range(0,len(reactant_feed)):
#				intensConst['inten_'+str(jk)] = Constant(reactant_feed[jk]*Inert_pulse_conc)
#			if sens_type == 'initial':
#				#controls = []
#				for jk in range(0,len(reactant_feed)):
#					controls.append(Control(intensConst['inten_'+str(jk)]))
#				
#
#			reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
#			
#			timeConst = {}
#
#			for jk in range(0,len(reactant_time)):
#				timeConst['sST_'+str(jk)] = Constant(round(reactant_time[jk]+0.001,6))
#
#			fnew = eval(pulse_functions(reac_input['reactions_test'],reac_input['Number of Inerts']))
#			F += fnew
#
#		if doe_form_surf == True:
#
#			inComp = {}
#			reac_input['Initial Surface Composition'].reverse()
#
#			for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):
#				inComp[z_num] = Constant(z_sur)
#			
#			F += eval(surface_functions(reac_input['reactions_test'],reac_input['Number of Inerts']))		
#
##############################################################
#							############### DEFINE INITIAL COMPOSITION ##################
#							#############################################################
#						
#							if k_pulse == 0 and round(t,6) == 0:
#								if doe_form_surf == False:
#									print('doe_form_surf false')
#									if additionalCells == 0:
#			
#										for z in range(mp.ceil((Mesh122)*reac_input['Mesh Size'])+1,int((Mesh222)*reac_input['Mesh Size'])+2):
#											if ',' in str(reac_input['Initial Surface Composition']):
#												for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
#													if int((Mesh122)*reac_input['Mesh Size'])-1 <= 0:
#														u_n.vector()[(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
#													else:
#														u_n.vector()[z*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
#			
#											else:
#												if int((Mesh122)*reac_input['Mesh Size'])-1 <= 0:
#													u_n.vector()[(all_molecules)-(2)] = float(species_pulse_list)
#												else:
#													u_n.vector()[z*(all_molecules)-(2)] = float(species_pulse_list)
#			
#									else:
#			
#										transTest1 = mp.ceil((Mesh122)*reac_input['Mesh Size'])+1
#										if ',' in str(reac_input['Initial Surface Composition']):
#											for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
#												u_n.vector()[transTest1*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)*0.5
#			
#										transTest2 = int((Mesh122)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))
#										if ',' in str(reac_input['Initial Surface Composition']):
#											for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
#												u_n.vector()[transTest2*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)*0.5
#			
#			
#										for z in range(mp.ceil((Mesh122)*reac_input['Mesh Size'])+2,int((Mesh122)*reac_input['Mesh Size'])+meshCells*2**(int(reac_input['Catalyst Mesh Density']))):
#											if ',' in str(reac_input['Initial Surface Composition']):
#												for z_num,z_sur in enumerate(reac_input['Initial Surface Composition']):	
#													if int((Mesh122)*reac_input['Mesh Size'])-1 <= 0:
#														u_n.vector()[(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
#													else:
#														u_n.vector()[z*(all_molecules)-(int(reac_input['Number of Inerts'])+z_num+1)] = float(z_sur)
#		
#											else:
#												if int((Mesh122)*reac_input['Mesh Size'])-1 <= 0:
#													u_n.vector()[(all_molecules)-(2)] = float(species_pulse_list)
#												else:
#													u_n.vector()[z*(all_molecules)-(2)] = float(species_pulse_list)
#
#meshCells = int((Mesh222)*reac_input['Mesh Size']) - mp.ceil((Mesh122)*reac_input['Mesh Size'])