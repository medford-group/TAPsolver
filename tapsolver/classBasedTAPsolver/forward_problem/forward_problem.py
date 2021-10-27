
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
from fenics_adjoint import *
import sys
import pandas as pd
import numpy as np
import math as mp
#import dijitso
import time
import os
import scipy
import copy
import warnings

from structures import *
from file_io import *
from mechanism_construction import *
from reactor_species import *
from reference_parameters import *

warnings.simplefilter(action='ignore', category=FutureWarning)
set_log_level(30)
tol = 1E-20
runge_kutta_approach = False

def forward_problem(pulse_time, TAPobject_data: TAPobject):
	
	simplifiedTimeStep = False

	######################## DEFINING INITIAL CONSTANTS #############################

	time_steps = pulse_time*1000
	dk = Constant(pulse_time/time_steps)
	cat_location = 1 - TAPobject_data.reactor.catalyst_center_fraction
	standard_parameters = load_standard_parameters()

	dx_r = TAPobject_data.reactor.total_length/TAPobject_data.mesh
	point_volume = dx_r*TAPobject_data.reactor.cross_sectional_radius*TAPobject_data.reactor.zone_voids[0]
	reference_pulse_concentration = TAPobject_data.reactor_species.reference_pulse_size/point_volume


	mesh = UnitIntervalMesh(int(TAPobject_data.mesh))
	#kinetic_information = load_kinetics(mechanism_data)
	#print(kinetic_information)


	cfDict = {}
	
	roundedMesh2 = ((1-cat_location) + 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh
	Mesh2 = round(((1-cat_location) + 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh)
	Mesh22 = round(((1-cat_location) + 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh)/TAPobject_data.mesh
	Mesh222 = round(((cat_location) + 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh)/TAPobject_data.mesh

	roundedMesh1 = ((1-cat_location) - 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh
	Mesh1 = round(((1-cat_location) - 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh)
	Mesh12 = round(((1-cat_location) - 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh)/TAPobject_data.mesh
	Mesh122 = round(((cat_location) - 0.5*TAPobject_data.reactor.zone_lengths[0])*TAPobject_data.mesh)/TAPobject_data.mesh

	if Mesh2 != roundedMesh2 or Mesh1 != roundedMesh1:
		print('Warning: Catalyst zone will be refined and rounded to the nearest whole mesh point!')
		trueMesh = (roundedMesh2 - roundedMesh1)/TAPobject_data.mesh
		newMesh = (Mesh2 - Mesh1)/TAPobject_data.mesh
		print()
		print('New Catalyst Fraction = '+str(newMesh))
		print('Old Catalyst Fraction = '+str(trueMesh))
		percentChange = abs(round(100*(trueMesh - newMesh)/trueMesh,2))
		print('Change = '+str(percentChange)+'%')
		print()
		if percentChange > 4:
			print('Consider refining the mesh to improve the accuracy of the simulation!')
			sys.exit()
		
	for jayz in range(0,TAPobject_data.catalyst_mesh_density+1):
		class thin_zoneTest(SubDomain):
			def inside(self, x, on_boundary):
				return between(x[0], ((Mesh12), (Mesh22)))
		
		thin_zoneTest = thin_zoneTest()
		cfDict[jayz] = MeshFunction("bool",mesh,1)
		
		thin_zoneTest.mark(cfDict[jayz],jayz)
		mesh = refine(mesh,cfDict[jayz],True)
	
	P1 = FiniteElement('CG',mesh.ufl_cell(),1)

	elements = []
	for k in range(0,len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+len(TAPobject_data.reactor_species.inert_gasses)):
		elements.append(P1)

	element = MixedElement(elements)
	V = FunctionSpace(mesh,element)
	V_du = FunctionSpace(mesh,P1)

	u = Function(V)
	u_n = Function(V)
	u_temp = Function(V)

	v_d = {}
	u_d = {}
	u_nd = {}
	
	tempA = TestFunctions(V)
	tempB = split(u)
	tempC = split(u_n)
	
	for knum,k in enumerate(list(TAPobject_data.reactor_species.gasses.keys())):
		v_d['v_'+k] = tempA[knum]
		u_d['u_'+k] = tempB[knum]
		u_nd['u_n'+k] = tempC[knum]
	for knum,k in enumerate(list(TAPobject_data.reactor_species.adspecies.keys())):
		v_d['v_'+k] = tempA[len(TAPobject_data.reactor_species.gasses)+knum]
		u_d['u_'+k] = tempB[len(TAPobject_data.reactor_species.gasses)+knum]
		u_nd['u_n'+k] = tempC[len(TAPobject_data.reactor_species.gasses)+knum]
	for knum,k in enumerate(list(TAPobject_data.reactor_species.inert_gasses.keys())):
		v_d['v_'+k] = tempA[len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+knum]
		u_d['u_'+k] = tempB[len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+knum]
		u_nd['u_n'+k] = tempC[len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+knum]


	meshCells = int((Mesh22)*TAPobject_data.mesh) - mp.ceil((Mesh12)*TAPobject_data.mesh)
	totalNumCells = TAPobject_data.mesh+meshCells*2**(TAPobject_data.catalyst_mesh_density)-meshCells
	
	if TAPobject_data.reactor_species.advection > 0:
			
		W = VectorFunctionSpace(mesh, 'P', 1)
		advTerm = Function(W)
		advMulti = Constant(TAPobject_data.reactor_species.advection)
		advValue = Constant(1)
		advTerm.vector()[:] = advValue

		advTerm.vector()[totalNumCells-0] = 0

	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], ((Mesh12)-1/TAPobject_data.mesh, (Mesh22)+1/TAPobject_data.mesh))
				
	class singlePoint(SubDomain):
		def inside(self,x,on_boundary):
			return between(x[0], (((1-cat_location) - 1/totalNumCells), ((1-cat_location) + 1/totalNumCells)))
	
	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,1)
	thin_zone.mark(domains,1)
	
	newBoundaries = domains.mesh().coordinates().transpose().tolist()[0]
	additionalCells = len(newBoundaries[int(TAPobject_data.mesh+1):])
	
	centralPoint = singlePoint()
	boundary_parts0 = MeshFunction("size_t", mesh,0)
	boundary_parts0.set_all(0)
	centralPoint.mark(boundary_parts0, 1)
	
	dx = Measure("dx",subdomain_data=domains)
	dT = Measure("dx",subdomain_data=boundary_parts0)
		
	# Define inlet and outlet boundaries
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
	
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)
		
	dz = 1/TAPobject_data.mesh
	
	class integration_section(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (1-dz,1.0))
	
	right = CompiledSubDomain("near(x[0], 1.)")
	boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
	right.mark(boundary_parts, 1)
	
	bcs = []
	for k in range(0,len(TAPobject_data.reactor_species.gasses)):
		bcs.append(DirichletBC(V.sub(k),Constant(0),boundary_R))
	
	for k in range(0,len(TAPobject_data.reactor_species.inert_gasses)):
		bcs.append(DirichletBC(V.sub(len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.inert_gasses)-(1+k)),Constant(0),boundary_R))
	
	F_new = construct_f_equation(TAPobject_data)
	
	theta = 1
	Ftemp = eval(F_new)
	theta = 0.5
	F = eval(F_new)
	#except NameError:
	#	errorOutput(reac_input['reactions_test'])

	#controls = []
	b0Test2 = Expression('x[0] < 0.002500001 ? 0.5 : 0', degree=0)
			
	#pulseIntensities = [Inert_pulse_conc/2,Inert_pulse_conc,Inert_pulse_conc*0,Inert_pulse_conc/2,Inert_pulse_conc/2]
	#pulseTimes = [0.01,0.001,0.001,0.001,0.001]
	reactant_feed = list(map(float, reac_input['Pulse Size'].split(',')))

	intensityFunctions ={}
	intensConst = {}
	for jk in range(0,len(reactant_feed)):
		intensConst['inten_'+str(jk)] = Constant(reactant_feed[jk]*Inert_pulse_conc)
	if sens_type == 'initial':
		#controls = []
		for jk in range(0,len(reactant_feed)):
			controls.append(Control(intensConst['inten_'+str(jk)]))
		
	reactant_time = list(map(float, reac_input['Pulse Time'].split(',')))
			
	timeConst = {}

	for jk in range(0,len(reactant_time)):
		timeConst['sST_'+str(jk)] = Constant(round(reactant_time[jk]+0.001,6))

	fnew = eval(pulse_functions(reac_input['reactions_test'],reac_input['Number of Inerts']))
	F += fnew

	sys.exit()
	# Standard Constants
	
	

	######################## DEFINING OUTPUT FOLDER INFORMATION #############################

	def generateFolder(path_name):
		try:  
			os.mkdir(path_name)
		except OSError:  
			pass
		else: 
			pass

	path = './'+reactor_data.output_name+'_folder/'
	generateFolder(path)
			
	path_3 = './'+reactor_data.output_name+'_folder/flux_data/'
	generateFolder(path_3)

	######################## READING KINETIC PARAMETER INFORMATION ###########################







	

	######################## Define  ###########################



	# Construct the rate expression based on the microkinetic model
	necessary_values, rate_array, rev_irr = make_f_equation(reac_input['reactions_test'],reac_input['Number of Reactants'],'tap',reac_input['Number of Inerts'],reac_input['Advection'],arrForward,arrBackward,gForward,reac_input['linked names'],linkForward,linkBackard,fit_temperature)



'''
			for j_species in range(0,monitored_gas+int(reac_input['Number of Inerts'])):
				tempDict = np.transpose(dictionary_of_numpy_data[legend_label[j_species]])
				np.savetxt('./'+reac_input['Output Folder Name']+'_folder/flux_data/'+legend_label[j_species]+'.csv', tempDict, delimiter=",")


try:
			theta = 1
			Ftemp = eval(necessary_values['F'])
			theta = 0.5
			F = eval(necessary_values['F'])
		except NameError:
			errorOutput(reac_input['reactions_test'])
	
		#############################################################
		############# DEFINE VARIATIONAL PROBLEM ####################
		#############################################################
	
		J = derivative(F,u)
		Jtemp = derivative(Ftemp,u)
	
		# Define a constrained variational problem solver
		reac_input['Variational Solver'] = 'newton'
	
		if reac_input['Variational Solver'].lower() == 'constrained':
			snes_solver_parameters = {"nonlinear_solver": "snes","snes_solver": {"linear_solver": "lu","line_search":'basic',"maximum_iterations": 10,"report": False,"error_on_nonconvergence": False}}
			
			lower = Function(V)
			upper = Function(V) 
			
			ninfty = Function(V); ninfty.vector()[:] = 0
			pinfty = Function(V); pinfty.vector()[:] =  np.infty
	
			
			problem = NonlinearVariationalProblem(F,u,bcs,J)
			
			problem.set_bounds(ninfty,pinfty)
	
			solver = NonlinearVariationalSolver(problem)
			solver.parameters.update(snes_solver_parameters)
		
		# Define a newton variational problem solver
		elif reac_input['Variational Solver'].lower() == 'newton':
			problem = NonlinearVariationalProblem(F,u,bcs,J)
			solver = NonlinearVariationalSolver(problem)
	
			problemtemp = NonlinearVariationalProblem(Ftemp,u,bcs,Jtemp)
			solvertemp = NonlinearVariationalSolver(problemtemp)
	
			solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
			solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10


		for k_pulse in range(0,int(reac_input['Number of Pulses'])):
	
			start_time = time.time()
			
			if sampling == False:
				print("")
				print("Simulation Status: "+"Pulse #"+str(k_pulse+1))
	
			species_pulse_list = reactant_feed#reactant_feed[current_reac_set-1]
			
			if type(reactant_time[current_reac_set-1]) == list:
				species_time = reactant_time[current_reac_set-1]
			else:
				species_time = reactant_time
			
			if reac_input['Knudsen Test'].lower() == 'true':
				knudsenTest(legend_label[(int(len(legend_label))-int(reac_input['Number of Inerts'])):],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],reac_input['Reference Pulse Size'],species_pulse_list[(int(len(legend_label))-int(reac_input['Number of Inerts'])):])
	 
			dt = reac_input['Pulse Duration']/reac_input['Time Steps']
			
			for kTimeStep,kTime in enumerate(species_time.copy()):
				tNew = dt*round(float(kTime)/dt)
				species_time[kTimeStep] = round(tNew,6)		
			
			

	
			graph_data = {}
			u_graph_data = {}
			
			for k_gasses in range(0,necessary_values['molecules_in_gas_phase']):
				graph_data['conVtime_'+str(k_gasses)] = []
				u_graph_data['conVtime_'+str(k_gasses)] = []
			
			for kjc in range(0,int(reac_input['Number of Inerts'])):
				graph_data['conVtime_'+str((necessary_values['gas_num']+1)+necessary_values['surf_num']+kjc)] = []
			
			surf_data['timing'] = []
			
			for j_gasses in range(necessary_values['molecules_in_gas_phase'],len(necessary_values['reactants'])-1):
				surf_data['conVtime_'+str(j_gasses)] = []
			
			graph_data['timing'] = []
	
			cat_data['timing'] = []
			
			for z_gasses in range(0,all_molecules+1):
				cat_data['conVtime_'+str(z_gasses)] = []
	
			cat_data['rawData'] = []
			
			t = 0

			
			test_new = project(-u.dx(0),V)
			test_new_split = test_new.split(deepcopy=False)
	
			w_new = Expression("1", degree=0)
			w_new2 = interpolate(w_new,V_du)
	
			x_dim = list(range(0, int(reac_input['Mesh Size'])+1))
			
			cum_sum = 0
t = 0


			#############################################################
			############ SOLVE PDEs AT EACH TIME STEP ###################
			#############################################################

			# Design of experiment form
	
			while t <= reac_input['Pulse Duration']:
				if round(t/0.001,4).is_integer() == True:
					graph_data['timing'].append(t)
				
				for k in range(0,monitored_gas):
					if round(t/0.001,4).is_integer() == True:
						new_val = (to_flux[k]*( u_n.vector().get_local()[(all_molecules)+k]))
						graph_data['conVtime_'+str(k)].append((new_val))
	
				for kjc in range(0,int(reac_input['Number of Inerts'])):
					if round(t/0.001,4).is_integer() == True:
						new_val = ((to_flux[monitored_gas+kjc]*u_n.vector().get_local()[2*(all_molecules+1)-2-(int(reac_input['Number of Inerts'])-kjc)]))
						graph_data['conVtime_'+str(all_molecules-(int(reac_input['Number of Inerts'])-kjc))].append((new_val))
				
				mesh_size = reac_input['Mesh Size']
		
				#############################################################
				######## STEP FOR CONSTRUCTING OBJECTIVE FUNCTION ###########
				#############################################################



				if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or sampling == True or ((sens_type == 'total' or sens_type == 'initial') and reac_input['Sensitivity Analysis'].lower() == 'true') or fit_temperature == True:
					# Define the objective function 	
					objectiveAnalysis = True
					selectivity = False
					conversion = False
					yieldValue = False
	

				if simplifiedTimeStep == False:
					if runge_kutta_approach == False:	
						if round(t,6) not in species_time:
							#timeDiff = round(t,6)
							try:
								if reac_input['Fit Parameters'].lower() == 'true' or reac_input['Uncertainty Quantification'].lower() == 'true' or reac_input['Fit Inert'].lower() == 'true' or reac_input['Sensitivity Analysis'].lower() == 'true':
									# C.N. method
									if t > 0.0011+timeDiff:
										solver.solve()
									# B.E. method
									else:
										solvertemp.solve()
		
										if round(t,6) == round(0.001+timeDiff,6):
											dt = reac_input['Pulse Duration']/reac_input['Time Steps']
											dk.assign(dt)
											u_n.assign(u)
											solver.solve()
								else:
									if t > 0.0011+timeDiff:
										solver.solve(annotate = False)
									else:
										solvertemp.solve(annotate=False)
		
										if round(t,6) == round(0.001+timeDiff,6):
											dt = reac_input['Pulse Duration']/reac_input['Time Steps']
											dk.assign(dt)
											u_n.assign(u)
											solver.solve()
								
								#############################################################
								################### STORE THIN-ZONE DATA ####################
								#############################################################
								
								if reac_input['Thin-Zone Analysis'].lower() == 'true':
									
									arrayNew = np.vstack((cat_data['rawData'],u_n.vector().get_local()))
									cat_data['rawData'] = arrayNew
	
							except RuntimeError:
								print('Time Step Failure')
								sys.exit()
						
						else:
							if dt == 0.0001:
								pass
							else:
								dt = 0.0001
								dk.assign(dt)
		
							if round(t,6) in species_time:
								timeDiff = round(t,6)
		
							if reac_input['Reactor Type'] == 'tap':
							
								if reac_input['Fit Parameters'].lower() == 'true':
									if t == 0:
										for k in range(0,int(reac_input['Number of Reactants'])):
											u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+additionalCells))+k-(all_molecules)] = -1e-10
		
								if doe_form_pulse == False:
									for k in range(0,int(reac_input['Number of Reactants'])):
									
										if species_time[k] == round(t,6):
											u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells)))+k-(all_molecules)] = float(species_pulse_list[k])*Inert_pulse_conc#list_species_pulse[k]
									# Define: Pulse
									for k_step in range(0,int(reac_input['Number of Inerts'])):
										if species_time[-1-k_step] == round(t,6):
											u_n.vector()[int((all_molecules)*((reac_input['Mesh Size']+1)+(meshCells*2**(int(reac_input['Catalyst Mesh Density']))-meshCells))-1-k_step)] = float(species_pulse_list[-1-k_step])*Inert_pulse_conc###??? Added the porosity contribution
								

	 							#############################################################
								################### STORE THIN-ZONE DATA ####################
								#############################################################
		
								if reac_input['Thin-Zone Analysis'].lower() == 'true':
								
									if round(t,6) != 0:
										arrayNew = np.vstack((cat_data['rawData'],u_n.vector().get_local()))
										cat_data['rawData'] = arrayNew
									else:
										cat_data['rawData'] = u_n.vector().get_local()
'''