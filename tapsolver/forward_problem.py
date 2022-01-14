# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

from fenics import *
from fenics_adjoint import *
import sys
import pandas as pd
import numpy as np
import math as mp
import time
import os
import scipy
import copy
import warnings
import matplotlib.pyplot as plt

#from structures import *
from .define_adspecies import define_adspecies
from .define_gas import define_gas
from .display_gasses import display_gasses
from .display_surface import display_surface
from .experimental_data import experimental_data
from .mechanism import mechanism
from .reactor import reactor
from .reactor_species import reactor_species
#from .read_old_input import read_old_input
from .TAPobject import TAPobject

#from file_io import *
from .new_experiments import new_experiments
from .read_CSV_mechanism import read_CSV_mechanism
from .read_CSV_reactor import read_CSV_reactor
from .read_CSV_reactor_species import read_CSV_reactor_species
from .read_experimental_data_object import read_experimental_data_object
from .read_mechanism_object import read_mechanism_object
from .read_reactor_object import read_reactor_object
from .read_reactor_species_object import read_reactor_species_object 
from .read_TAPobject import read_TAPobject 
from .read_transient_sensitivity import read_transient_sensitivity 
from .save_object import save_object
#from vary_input_file import vary_input_file

#from mechanism_construction import *
#from construct_batch_equation import make_batch_equation
from .construct_f_equation import construct_f_equation
from .construct_f_equation_multiple_experiments import construct_f_equation_multiple_experiments
from .construct_rate_equations import rateEqs
from .display_elementary_processes import display_elementary_processes
from .elementary_process import elementary_process
from .elementary_process_details import elementary_process_details
from .mechanism_constructor import mechanism_constructor
from .mechanism_reactants import mechanism_reactants

#from reactor_species import *
#from reference_parameters import *
from .reference_parameters import load_standard_parameters

#from simulation_notes import *
from .timing_details import *
from .error_details import *
from .generate_folders import *

#from inverse_problem import *
from .define_fitting_species import curveFitting
#from point_objective import point_objective
from .std_objective import stdEstablishment
from .total_objective import curveFitting

import jsonpickle
import json
import ufl
import dijitso

warnings.simplefilter(action='ignore', category=FutureWarning)
set_log_level(30)
tol = 1E-20
runge_kutta_approach = False
standard_parameters = load_standard_parameters()

def forward_problem(pulse_time, pulse_number, TAPobject_data_original: TAPobject):
		
	TAPobject_data = copy.deepcopy(TAPobject_data_original)
	tape2 = Tape()
	tape2.clear_tape()
	set_working_tape(tape2)
	#parameter_scale = 120#Constant(120)
	#print(TAPobject_data.data_name)
	#print(TAPobject_data.output_name)
	#sys.exit()
	if TAPobject_data.data_name != None:
		
		TAPobject_data.experimental_data = read_experimental_data_object(TAPobject_data.data_name)#output_data
		output_data = read_experimental_data_object(TAPobject_data.data_name)

	species_time = []
	for j in TAPobject_data.reactor_species.gasses:
		species_time.append(TAPobject_data.reactor_species.gasses[j].delay)
	for j in TAPobject_data.reactor_species.inert_gasses:
		species_time.append(TAPobject_data.reactor_species.inert_gasses[j].delay)			
	
	TAPobject_data.reactor_species.reference_mass = Constant(TAPobject_data.reactor_species.reference_mass)
	if type(TAPobject_data.reactor_species.temperature) == dict:
		for j_temp_number, j_temps in enumerate(TAPobject_data.reactor_species.temperature): 
			TAPobject_data.reactor_species.temperature[j_temp_number] = Constant(TAPobject_data.reactor_species.temperature[j_temp_number])
	else:
		TAPobject_data.reactor_species.temperature = Constant(TAPobject_data.reactor_species.temperature)
	
	TAPobject_data.reactor_species.reference_temperature = Constant(TAPobject_data.reactor_species.reference_temperature)
	TAPobject_data.reactor.zone_voids[0] = Constant(TAPobject_data.reactor.zone_voids[0])
	TAPobject_data.reactor_species.catalyst_diffusion = Constant(TAPobject_data.reactor_species.catalyst_diffusion)
	TAPobject_data.reactor.zone_voids[1] = Constant(TAPobject_data.reactor.zone_voids[1])
	standard_parameters['kbt'] = Constant(standard_parameters['kbt'])
	standard_parameters['h'] = Constant(standard_parameters['h'])
	standard_parameters['Rgas'] = Constant(standard_parameters['Rgas'])
	TAPobject_data.reactor.reactor_radius = Constant(TAPobject_data.reactor.reactor_radius)
			
	species_order_dictionary = {}
	total_species = 0
	for k in TAPobject_data.reactor_species.gasses:
		species_order_dictionary[k] = total_species
		TAPobject_data.reactor_species.gasses[k].mass = Constant(TAPobject_data.reactor_species.gasses[k].mass)
		TAPobject_data.reactor_species.gasses[k].delay = Constant(TAPobject_data.reactor_species.gasses[k].delay)
		TAPobject_data.reactor_species.gasses[k].intensity = Constant(TAPobject_data.reactor_species.gasses[k].intensity)
		total_species += 1
	for k in TAPobject_data.reactor_species.adspecies:
		species_order_dictionary[k] = total_species
		TAPobject_data.reactor_species.adspecies[k].concentration = Constant(TAPobject_data.reactor_species.adspecies[k].concentration)
		total_species += 1
	for k in TAPobject_data.reactor_species.inert_gasses:
		species_order_dictionary[k] = total_species
		TAPobject_data.reactor_species.inert_gasses[k].mass = Constant(TAPobject_data.reactor_species.inert_gasses[k].mass)
		TAPobject_data.reactor_species.inert_gasses[k].delay = Constant(TAPobject_data.reactor_species.inert_gasses[k].delay)
		TAPobject_data.reactor_species.inert_gasses[k].intensity = Constant(TAPobject_data.reactor_species.inert_gasses[k].intensity)	
		total_species += 1
	
	if TAPobject_data.mechanism.reactants != []:
		for k,z in enumerate(TAPobject_data.mechanism.rate_array):
			if TAPobject_data.mechanism.elementary_processes[k].forward.use == 'G':
				TAPobject_data.mechanism.elementary_processes[k].forward.Ga = Constant(TAPobject_data.mechanism.elementary_processes[k].forward.Ga)
				TAPobject_data.mechanism.elementary_processes[k].forward.dG = Constant(TAPobject_data.mechanism.elementary_processes[k].forward.dG)
			elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'E':
				TAPobject_data.mechanism.elementary_processes[k].forward.Ao = Constant(TAPobject_data.mechanism.elementary_processes[k].forward.Ao)
				TAPobject_data.mechanism.elementary_processes[k].forward.Ea = Constant(TAPobject_data.mechanism.elementary_processes[k].forward.Ea)
			elif TAPobject_data.mechanism.elementary_processes[k].forward.use == 'k':
				TAPobject_data.mechanism.elementary_processes[k].forward.k = Constant(TAPobject_data.mechanism.elementary_processes[k].forward.k)									
			#if TAPobject_data.mechanism.elementary_processes[k].backward.use == 'G':
			#	TAPobject_data.mechanism.elementary_processes[k].backward.Ga = Constant(TAPobject_data.mechanism.elementary_processes[k].backward.Ga)
			#	TAPobject_data.mechanism.elementary_processes[k].backward.dG = Constant(TAPobject_data.mechanism.elementary_processes[k].backward.dG)
			if TAPobject_data.mechanism.elementary_processes[k].backward.use == 'E':
			#elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'E':
				TAPobject_data.mechanism.elementary_processes[k].backward.Ao = Constant(TAPobject_data.mechanism.elementary_processes[k].backward.Ao)
				TAPobject_data.mechanism.elementary_processes[k].backward.Ea = Constant(TAPobject_data.mechanism.elementary_processes[k].backward.Ea)
			elif TAPobject_data.mechanism.elementary_processes[k].backward.use == 'k':
				try:
					TAPobject_data.mechanism.elementary_processes[k].backward.k = Constant(TAPobject_data.mechanism.elementary_processes[k].backward.k)
				except:
					pass
		if TAPobject_data.mechanism.kinetic_links != {}:
			for j in TAPobject_data.mechanism.kinetic_links:
				TAPobject_data.mechanism.kinetic_links[j] = Constant(TAPobject_data.mechanism.kinetic_links[j])

	time_steps = pulse_time*1000

	dk = Constant(pulse_time/time_steps)
	cat_location = 1 - TAPobject_data.reactor.catalyst_center_fraction

	dx_r = TAPobject_data.reactor.total_length/TAPobject_data.mesh
	point_volume = dx_r*TAPobject_data.reactor.cross_sectional_radius*TAPobject_data.reactor.zone_voids[0]
	reference_pulse_concentration = TAPobject_data.reactor_species.reference_pulse_size/point_volume
	##### Making sure that I define the controls
	controls = []
	if TAPobject_data.optimize == True or TAPobject_data.tangent_linear_sensitivity == True:
		try:
			if len(TAPobject_data.parameters_of_interest) > 1:
				pass
		except:
			TAPobject_data.parameters_of_interest = TAPobject_data.parameters_of_interest[0]
		for j in TAPobject_data.parameters_of_interest:
			#print(TAPobject_data.parameters_of_interest)
			#controls.append(Control(TAPobject_data.reactor_species.temperature[j_temp_number]))
			#controls.append(Control(TAPobject_data.mechanism.elementary_processes[0].forward.k))
			controls.append(Control(eval(j)))
	#controls = [Control(TAPobject_data.reactor_species.inert_gasses['Ar'].intensity)]
	#controls = [Control(TAPobject_data.mechanism.elementary_processes[0].forward.k)]

	mesh = UnitIntervalMesh(int(TAPobject_data.mesh))

	cfDict = {}
	roundedMesh2 = ((1-cat_location) + 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh
	Mesh2 = round(((1-cat_location) + 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh)
	Mesh22 = round(((1-cat_location) + 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh)/TAPobject_data.mesh
	Mesh222 = round(((cat_location) + 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh)/TAPobject_data.mesh

	roundedMesh1 = ((1-cat_location) - 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh
	Mesh1 = round(((1-cat_location) - 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh)
	Mesh12 = round(((1-cat_location) - 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh)/TAPobject_data.mesh
	Mesh122 = round(((cat_location) - 0.5*TAPobject_data.reactor.length_fractions[1])*TAPobject_data.mesh)/TAPobject_data.mesh
	
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

	meshCells = int((Mesh22)*TAPobject_data.mesh) - mp.ceil((Mesh12)*TAPobject_data.mesh)
	totalNumCells = TAPobject_data.mesh+meshCells*2**(TAPobject_data.catalyst_mesh_density)-meshCells

	top = mp.ceil((Mesh122)*TAPobject_data.mesh)+1 
	bottom = int((Mesh122)*TAPobject_data.mesh)+meshCells*2**(TAPobject_data.catalyst_mesh_density)-1
	catalyst_center_cell = int((top+bottom)/2)
	
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
		bcs.append(DirichletBC(V.sub(len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+k),Constant(0),boundary_R))
	
	if type(TAPobject_data.reactor_species.temperature) != dict:
		F_new = construct_f_equation(TAPobject_data)
	else:
		F_new = construct_f_equation_multiple_experiments(TAPobject_data)

	constantT = Constant(0)
	remove_surface = Constant(1)
	b0Test2 = Expression('x[0] < 0.002500001 ? 0.5 : 0', degree=0)

	Fpulses = ''
	for knum,k in enumerate(TAPobject_data.reactor_species.gasses):
		##TAPobject_data.reactor_species.gasses[k].intensity = Constant(TAPobject_data.reactor_species.gasses[k].intensity)
		#Fpulses += "-TAPobject_data.reactor_species.gasses['"+k+"'].intensity*reference_pulse_concentration*b0Test2*exp(-(constantT - round(TAPobject_data.reactor_species.gasses['"+k+"'].delay+0.001,6))*(constantT - round(TAPobject_data.reactor_species.gasses['"+k+"'].delay+0.001,6))/(4*0.00000000001))*v_d['v_"+k+"']*dx"
		Fpulses += "-TAPobject_data.reactor_species.gasses['"+k+"'].intensity*reference_pulse_concentration*b0Test2*exp(-(constantT - round(TAPobject_data.reactor_species.gasses['"+k+"'].delay,6))*(constantT - round(TAPobject_data.reactor_species.gasses['"+k+"'].delay,6))/(0.00000000001))*v_d['v_"+k+"']*dx"
		if knum < len(TAPobject_data.reactor_species.gasses)-1:
			Fpulses += ' + '
	for knum,k in enumerate(TAPobject_data.reactor_species.inert_gasses):
		##TAPobject_data.reactor_species.inert_gasses[k].intensity = Constant(TAPobject_data.reactor_species.inert_gasses[k].intensity)
		#Fpulses += "-TAPobject_data.reactor_species.inert_gasses['"+k+"'].intensity*reference_pulse_concentration*b0Test2*exp(-(constantT - round(TAPobject_data.reactor_species.inert_gasses['"+k+"'].delay+0.001,6))*(constantT - round(TAPobject_data.reactor_species.inert_gasses['"+k+"'].delay+0.001,6))/(4*0.00000000001))*v_d['v_"+k+"']*dx"
		Fpulses += "-TAPobject_data.reactor_species.inert_gasses['"+k+"'].intensity*reference_pulse_concentration*b0Test2*exp(-(constantT - round(TAPobject_data.reactor_species.inert_gasses['"+k+"'].delay,6))*(constantT - round(TAPobject_data.reactor_species.inert_gasses['"+k+"'].delay,6))/(0.00000000001))*v_d['v_"+k+"']*dx"
		if knum < len(TAPobject_data.reactor_species.inert_gasses)-1:
			Fpulses += ' + '

	for knum,k in enumerate(TAPobject_data.reactor_species.adspecies):
		##TAPobject_data.reactor_species.adspecies[k].concentration = Constant(TAPobject_data.reactor_species.adspecies[k].concentration)
		Fpulses += "-remove_surface*TAPobject_data.reactor_species.adspecies['"+k+"'].concentration*exp(-(constantT)*(constantT)/(0.00000000001))*v_d['v_"+k+"']*dx"
		if knum < len(TAPobject_data.reactor_species.adspecies)-1:
			Fpulses += ' + '
	
	# dk*("+"(TAPobject_data.reactor_species.inert_diffusion*sqrt(TAPobject_data.reactor_species.reference_mass*TAPobject_data.reactor_species.temperature)/sqrt(TAPobject_data.reactor_species.reference_temperature*TAPobject_data.reactor_species.gasses['"+species_name+"'].mass))

	F_new += Fpulses

	#!#! ADDING INVERSE ANALYSIS
	if TAPobject_data.tangent_linear_sensitivity == True:
		#TAPobject_data.reactor_species.gasses['CO'].intensity = Control(TAPobject_data.reactor_species.gasses['CO'].intensity)
		#TAPobject_data.mechanism.elementary_processes[0].forward.k["value"] = Control(TAPobject_data.mechanism.elementary_processes[0].forward.k["value"])
		#c = TAPobject_data.mechanism.elementary_processes[0].forward.k["value"]
		#c.tlm_value = TAPobject_data.mechanism.elementary_processes[0].forward.k["value"]
		sens_param = eval(TAPobject_data.parameters_of_interest[0])
		
		c = sens_param
		c.tlm_value = sens_param
		#c = TAPobject_data.reactor_species.gasses['CO'].intensity
		#c.tlm_value = TAPobject_data.reactor_species.gasses['CO'].intensity
	
		SV_du = FunctionSpace(mesh,P1)
		Sw_new = Expression('A',A=Constant(1),degree=0)
		Sw_new2 = interpolate(Sw_new,SV_du)
		Sw3 = project(Sw_new2,SV_du)
	
		sensFuncs = {}
		for k_gasses in TAPobject_data.reactor_species.gasses:
			sensFuncs[k_gasses] = []

	#if TAPobject_data.adjoint_sensitivitiy == True:
	#	#print(TAPobject_data.parameters_of_interest[0])
	#	controls = [TAPobject_data.reactor_species.gasses['CO'].intensity]
	#	# output_fitting = pointFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)


	theta = Constant(1)
	Ftemp = eval(F_new)
	theta = Constant(0.5)
	F = eval(F_new)

	J = derivative(F,u)
	Jtemp = derivative(Ftemp,u)
	
	variational_solver = 'newton'
	if variational_solver == 'constrained':
		snes_solver_parameters = {"nonlinear_solver": "snes","snes_solver": {"linear_solver": "lu","line_search":'basic',"maximum_iterations": 10,"report": False,"error_on_nonconvergence": False}}
			
		lower = Function(V)
		upper = Function(V) 
		
		ninfty = Function(V); ninfty.vector()[:] = 0
		pinfty = Function(V); pinfty.vector()[:] =  np.infty
	
		problem = NonlinearVariationalProblem(F,u,bcs,J)
		
		problem.set_bounds(ninfty,pinfty)
	
		solver = NonlinearVariationalSolver(problem)
		solver.parameters.update(snes_solver_parameters)
		
	elif variational_solver == 'newton':
		problem = NonlinearVariationalProblem(F,u,bcs,J)
		solver = NonlinearVariationalSolver(problem)
	
		problemtemp = NonlinearVariationalProblem(Ftemp,u,bcs,Jtemp)
		solvertemp = NonlinearVariationalSolver(problemtemp)

		solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10

	synthetic_data = {}
	synthetic_data['time'] = {}
	if TAPobject_data.store_flux_data == True:
		for j in TAPobject_data.reactor_species.gasses:
			synthetic_data[j] =  {}
		for j in TAPobject_data.reactor_species.inert_gasses:
			synthetic_data[j] =  {}
	if TAPobject_data.store_catalyst_data == True:
		for j in TAPobject_data.reactor_species.adspecies:
			synthetic_data[j] =  {}

	if TAPobject_data.tangent_linear_sensitivity == True:
		sensitivity_output = {}
			
		for k_sens in TAPobject_data.reactor_species.gasses:
			sensitivity_output[k_sens] = []

	if TAPobject_data.tangent_linear_sensitivity == True or TAPobject_data.adjoint_sensitivitiy == True  or TAPobject_data.optimize == True:
		osub = integration_section()
		domains = MeshFunction("size_t", mesh,0)
		domains.set_all(0)
		osub.mark(domains, 1)
		dP = Measure('vertex',domain = mesh, subdomain_data=domains)
	#!#!
	#for knum,k in enumerate(TAPobject_data.reactor_species.gasses):
	#	print(TAPobject_data.reactor_species.gasses[k].inert_diffusion)
	#for knum,k in enumerate(TAPobject_data.reactor_species.inert_gasses):
	#	print(TAPobject_data.reactor_species.inert_gasses[k].inert_diffusion)
	#print(TAPobject_data.output_name)
	#if TAPobject_data.output_name != None:
	#	path_4 = './'+TAPobject_data.output_name+'/'
	#	print(path_4)
	try: 
		os.mkdir('./'+TAPobject_data.output_name)
	except OSError:  
		pass
	
	#output_data = curveFitting(pulse_time,TAPobject_data)
	#TAPobject_data.reactor_species.gasses['CO'].sigma = 0.5

	#inverse_results = pd.DataFrame(index=range(1),columns=range(len()))
	
	def derivCB(j,dj,m):

		djv = [v.values()[0] for v in dj]
		mv = [v.values()[0] for v in m]
		print('Step Time: '+str(time.time()))
		print('Objective Value: '+str(j))
		print('Derivative Values '+str(djv))
		print('Parameter Values: '+str(mv))
			
		try:
			new_addition = pd.read_csv('./'+TAPobject_data.output_name+'/optimization_results.csv')
			name_columns = ['objective']
			value_row = [j]
			for jz_num,jz in enumerate(TAPobject_data.parameters_of_interest):
				name_columns.append(jz)
				value_row.append(mv[jz_num])
			new_addition.loc[len(new_addition)] = value_row
			new_addition.to_csv('./'+TAPobject_data.output_name+'/optimization_results.csv',index=False)
		
		except:
			name_columns = ['objective']
			value_row = [j]
			for jz_num,jz in enumerate(TAPobject_data.parameters_of_interest):
				name_columns.append(jz)
				value_row.append(mv[jz_num])
			percentile_list = pd.DataFrame([value_row],columns=name_columns)
			percentile_list.to_csv('./'+TAPobject_data.output_name+'/optimization_results.csv',index=False)
		

	#def derivCB(j,dj,m):
	#	print("Derivative Time: "+str(time.time() - start_time))
	#	with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIter.txt', 'w') as f:
	#		f.write("Objective Value: "+str(j_values))
	#		f.write('\n')
	#		f.write("Change: "+str(dj_values))
	#		f.write('\n')
	#		f.write("Constants: "+str(x_values))
#
#	#		f.close
#	#	print(j)
#	#	print(djv)
	#	print(mv)

	
	def hessCB(j,dj,m):
		it_times.append(time.time())
		j_values.append(j)
		djv = [v.values()[0] for v in dj]
		dj_values.append(djv)
		mv = [v.values()[0] for v in m]
		x_values.append(mv)
		print(time.time() - start_time)
		with open('./'+reac_input['Output Folder Name']+'_folder/fitting/optIterHess.txt', 'w') as f:
			f.write("Objective Value: "+str(j_values))
			f.write('\n')
			f.write("Change: "+str(dj_values))
			f.write('\n')
			f.write("Constants: "+str(x_values))
			f.close

		print(j)
		print(djv)
		print(mv)
	test_tlm_evaluation = []

	for k_pulse in range(0,pulse_number):
		print('Pulse #: '+str(k_pulse))
		t = 0

		if k_pulse > 0:
			constantT.assign(0.0)
			remove_surface.assign(0.0)
		dt = pulse_time/time_steps
		start_time = time.time()
		
		synthetic_data['time'][k_pulse] =  []
		if TAPobject_data.store_flux_data == True:
			for j in TAPobject_data.reactor_species.gasses:
				synthetic_data[j][k_pulse] =  []
			for j in TAPobject_data.reactor_species.inert_gasses:
				synthetic_data[j][k_pulse] =  []
		if TAPobject_data.store_catalyst_data == True:
			for j in TAPobject_data.reactor_species.adspecies:
				synthetic_data[j][k_pulse] =  []

		all_species = len(TAPobject_data.reactor_species.gasses) + len(TAPobject_data.reactor_species.inert_gasses) + len(TAPobject_data.reactor_species.adspecies)
		time_step = 0

		while t <= pulse_time:

			synthetic_data['time'][k_pulse].append(round(t,6))
			if TAPobject_data.store_flux_data == True:
				for knum,k in enumerate(TAPobject_data.reactor_species.gasses):
					synthetic_data[k][k_pulse].append(2*(float(TAPobject_data.reactor_species.gasses[k].inert_diffusion) /(float(dx_r))) * (float(TAPobject_data.reactor.reactor_radius)**2)*3.14159*( float(u_n.vector().get_local()[(all_species)+knum])))
				for knum,k in enumerate(TAPobject_data.reactor_species.inert_gasses):
					synthetic_data[k][k_pulse].append(2*(float(TAPobject_data.reactor_species.inert_gasses[k].inert_diffusion) /(float(dx_r))) * (float(TAPobject_data.reactor.reactor_radius)**2)*3.14159*( float(u_n.vector().get_local()[all_species+len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+knum])))		
			if TAPobject_data.store_catalyst_data == True:
				for knum,k in enumerate(TAPobject_data.reactor_species.adspecies):
					synthetic_data[k][k_pulse].append(float(u_n.vector().get_local()[(catalyst_center_cell*(len(TAPobject_data.reactor_species.gasses) + len(TAPobject_data.reactor_species.adspecies) + len(TAPobject_data.reactor_species.inert_gasses)))+all_species+len(TAPobject_data.reactor_species.gasses)+knum]))
					
			if TAPobject_data.store_catalyst_data == True:
				pass
			#for j in TAPobject_data.reactor_species.adspecies:
			#	synthetic_data[j][k_pulse] =  {}

			if TAPobject_data.optimize == True:
				for k_fitting in (TAPobject_data.gasses_objective): #  and TAPobject_data.inert_gasses_objective
					
					if round(t,6) in output_data['time'][0]:
						c_exp = output_data[k_fitting][0][output_data['time'][0].index(round(t,6))]
						#c_exp = output_data[k_fitting]['values'][output_data[k_fitting]['times'].index(round(t,6))]
						slope = (-c_exp)/(1/TAPobject_data.mesh)
						intercept = c_exp - ((1-(1/TAPobject_data.mesh))*slope)
						w_new = Expression('A*x[0]+B',A=Constant(slope),B=Constant(intercept),degree=0)
						w_new2 = interpolate(w_new,V_du)
						w3 = project(w_new2,V_du)
						
						try:
							if k_fitting in TAPobject_data.reactor_species.gasses:
								jfunc_2 += assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.reactor_species.gasses[k_fitting].sigma)**2)
							elif k_fitting in TAPobject_data.reactor_species.inert_gasses:
								jfunc_2 += assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.inert_gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.inert_gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.reactor_species.inert_gasses[k_fitting].sigma)**2)
								
						except UnboundLocalError:
							w_temp_2 = Expression('1',degree=0) 
							w_temp2_2 = interpolate(w_temp_2,V_du)
							w4_2 = project(w_temp2_2,V_du)	
							if k_fitting in TAPobject_data.reactor_species.gasses:
								jfunc_2 = assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.reactor_species.gasses[k_fitting].sigma)**2)	
							elif k_fitting in TAPobject_data.reactor_species.inert_gasses:
								jfunc_2 = assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.inert_gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.reactor_species.inert_gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.reactor_species.inert_gasses[k_fitting].sigma)**2)	
								
				if TAPobject_data.thermodynamic_constraints == True and t == 0:
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

						tempFunc[kip] = (tempFunc[kip]-w4[kip])*thermoWeight
						jfunc_2 += assemble(inner(tempFunc[kip],tempFunc[kip])*dx())
						print('Simulated Thermo Value')
						print(jfunc_2)


			#flux_data[pulse_number]['time'].append(t)
			#for knum,k in enumerate(list(TAPobject_data.reactor_species.gasses.keys())):
			#	flux_data[pulse_number][k].append( 2*(TAPobject_data.reactor_species.gasses[k].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159*( u_n.vector().get_local()[(all_species)+knum]))
			#for knum,k in enumerate(list(TAPobject_data.reactor_species.inert_gasses.keys())):
			#	flux_data[pulse_number][k].append( 2*(TAPobject_data.reactor_species.inert_gasses[k].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159*( u_n.vector().get_local()[all_species+len(TAPobject_data.reactor_species.gasses)+len(TAPobject_data.reactor_species.adspecies)+knum]))
			
			#graph_data['conVtime_'+str(k)].append((new_val))

			if round(t,6) not in species_time:
				try:
					if TAPobject_data.tangent_linear_sensitivity == True or TAPobject_data.adjoint_sensitivitiy == True or TAPobject_data.optimize == True:
						if t > 0.0011+timeDiff:
							solver.solve()
						else:
							solvertemp.solve()
							if round(t,6) == round(0.001+timeDiff,6):
								dt = pulse_time/time_steps
								dk.assign(dt)
								u_n.assign(u)
								solver.solve()
					else:
						
						if t > 0.0011+timeDiff:
							solver.solve(annotate = False)
						else:
							solvertemp.solve(annotate=False)
							if round(t,6) == round(0.001+timeDiff,6):
								dt = pulse_time/time_steps
								dk.assign(dt)
								u_n.assign(u)
								solver.solve()
										
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
						
				try:
					if TAPobject_data.tangent_linear_sensitivity == True or TAPobject_data.adjoint_sensitivitiy == True or TAPobject_data.optimize == True:
						solvertemp.solve()
					else:
						solvertemp.solve(annotate = False)
		
				except RuntimeError:
					print('Time Step Failure')
					sys.exit()
				

			#!#! ADDING INVERSE ANALYSIS
			if TAPobject_data.tangent_linear_sensitivity == True:
				for knum,k in enumerate(TAPobject_data.reactor_species.gasses):
					new_val = 2*(TAPobject_data.reactor_species.gasses[k].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159
					#print((2*(TAPobject_data.reactor_species.gasses[k].inert_diffusion /(dx_r)) * (float(TAPobject_data.reactor.reactor_radius)**2)*3.14159*( u_n.vector().get_local()[(all_species)+knum])))
					#print(assemble(inner((new_val*(u_n[knum])), Sw3)*dP(1) )) #/Constant(0.00758*2*0.999341)
					#sensFuncs[k].append(assemble( inner( (2*(TAPobject_data.reactor_species.gasses[k].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159*( u_n[knum])), Sw3/Constant(0.00758*2*0.999341) ) * dP(1) ))
					if k == 'CO':
						test_tlm_evaluation.append(assemble(inner((new_val*(u_n[knum])), Sw3)*dP(1) ))
					sensFuncs[k].append(assemble(inner((new_val*(u_n[knum])), Sw3)*dP(1) ))
					#sensFuncs[k].append(assemble( ( inner(u[knum], Sw3/Constant(0.0075000000000000015)) )* dx(1)))

			if TAPobject_data.adjoint_sensitivitiy == True:
				pass
				#print()
				#temp11 = (assemble( ( inner(( (TAPobject_data.reactor_species.gasses[k].inert_diffusion ) * (TAPobject_data.reactor.reactor_radius**2)*3.14159*u_n[0]), Sw3/Constant(0.00758*2*0.999341)) )* dP(1)))
				#temp12 = (2*(TAPobject_data.reactor_species.gasses[k].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.reactor_radius**2)*3.14159*( u_n.vector().get_local()[(all_species)+0]))
				#try:
				#	print(temp11/temp12)
				#except:
				#	pass
				#print(sensFuncs['CO'][-1])
			#!#!

			progressBar(t, pulse_time)

			Uvector = as_backend_type(u.vector()).get_local()
			Uvector[Uvector <= DOLFIN_EPS] = DOLFIN_EPS
			u.vector().set_local(Uvector)
			
			u_n.assign(u)
				
			t += dt
			constantT.assign(round(t,6))
			time_step += 1
		print(processTime(start_time))

	if TAPobject_data.objective_return == True:
		print('test')
		return jfunc_2

	if TAPobject_data.uncertainty == True:
		start_time = time.time()
		print()
		print('Calculating hessian. Could take some time.')
		hessFolder = TAPobject_data.output_name
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
		TAPobject_data.optimize == False
		sys.exit()

	if TAPobject_data.optimize == True:
		print('test 2')
		reference_time = time.time()
		rf_2 = ReducedFunctional(jfunc_2, controls,tape=tape2,derivative_cb_post=derivCB)#,derivative_cb_post=derivCB,hessian_cb_post=hessCB)
		
		low_bounds = []
		up_bounds = []

		for j_limit in TAPobject_data.parameters_of_interest:
			if 'dG' in j_limit:
				low_bounds.append(-np.inf)
				up_bounds.append(np.inf)
			else:
				low_bounds.append(0)
				up_bounds.append(np.inf)
	
		#for gt in range(0,len(controls)):
		#	low_bounds.append(0)
		#	up_bounds.append(np.inf)

		
		u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds), tol=1e-22, options={"ftol":1e-22,"gtol":1e-22})
		sys.exit()

	#!#! ADDING INVERSE ANALYSIS
	if TAPobject_data.tangent_linear_sensitivity == True:
		print()
		start_time = time.time()
		print('Evaluating Tape with Tangent Linear Method. Could take some time.')
		tape2.evaluate_tlm()
		print(processTime(start_time))
	
	if TAPobject_data.tangent_linear_sensitivity == True:
	
		for numEachSens,eachSens in enumerate(sensFuncs):
						
			newList = []
			for kSensNum, kSens in enumerate(sensFuncs[eachSens]):
				newValue = kSens.block_variable.tlm_value
				newList.append(newValue)
				df = pd.DataFrame(newList[:-1])
				df.to_csv('./dF_'+eachSens+'.csv',header=None,index=False)
			#np.savetxt('./dF_'+eachSens+'.txt',newList[:-1])
			#sys.exit()
		print(test_tlm_evaluation)
		return sensFuncs
	#!#!

	if 	TAPobject_data.gas_noise == True:
		beta_2 = 0.00270
		w_2 = 2*3.14159*70
		for k in synthetic_data:
			if k != 'time':
				for j in synthetic_data[k]:
					for z in range(0,len(synthetic_data[k][j])):
						try:
							synthetic_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.reactor_species.gasses[k].noise# +beta_2*np.cos(w_2*(k*dt))
							
						except:
							try:
								synthetic_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.reactor_species.inert_gasses[k].noise# +beta_2*np.cos(w_2*(k*dt))
							except:
								pass

	if	TAPobject_data.surface_noise == True:
		beta_2 = 0.00270
		w_2 = 2*3.14159*70
		for k in synthetic_data:
			if k != 'time':
				for j in synthetic_data[k]:
					for z in range(0,len(synthetic_data[k][j])):
						try:
							synthetic_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.reactor_species.adspecies[k].noise# +beta_2*np.cos(w_2*(k*dt))
						except:
							pass

	if TAPobject_data.store_flux_data == True or TAPobject_data.store_catalyst_data == True:
		#print('storing Data')
		#print(synthetic_data['C3H8'][0])
		#sys.exit()
		save_object(synthetic_data,'./'+TAPobject_data.output_name+'/TAP_experimental_data.json')
		new_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')

	if TAPobject_data.show_graph == True:
		fig, ax = plt.subplots()
		ax.set_ylabel('Flow (nmol/s)')
		ax.set_xlabel('Time (s)')
		for j in TAPobject_data.reactor_species.gasses:
			for k in synthetic_data[j]:
				if k == 0:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k],label=j)
				else:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k])

		for j in TAPobject_data.reactor_species.inert_gasses:
			for k in synthetic_data[j]:
				if k == 0:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k],label=j)		
				else:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k])
		plt.legend()
		plt.show()