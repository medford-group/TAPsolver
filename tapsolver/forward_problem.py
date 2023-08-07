
# Copyright 2021, Battelle Energy Alliance, LLC All Rights Reserved

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
from .etc import *
import jsonpickle
import json
import ufl
import dijitso

### From folders

from .data import *
#from experimental_data import experimental_data
from .mechanism import *#mechanism
from .reactor import *
from .species import *
from .TAPobject import *#TAPobject
from .read_files import *
from .construct_f_equation import *
from .new_construct_f_equation import *
from .construct_f_equation_multiple_experiments import *
from .etc import *

from .flux_graph import *

from fenics import *
from fenics_adjoint import *

warnings.simplefilter(action='ignore', category=FutureWarning)
set_log_level(30)
tol = 1E-20
runge_kutta_approach = False
standard_parameters = load_standard_parameters()

def forward_problem(pulse_time, pulse_number, TAPobject_data_original: TAPobject):
	TAPobject_data = copy.deepcopy(TAPobject_data_original)
	tape2 = Tape()
	tape2.clear_tape()
	derivative_eval = False

	TAPobject_data.species.reference_diffusions = TAPobject_data.reactor.reference_diffusions

	if (TAPobject_data.tangent_linear_sensitivity == True) or (TAPobject_data.adjoint_sensitivitiy == True)  or (TAPobject_data.optimize == True):	
		derivative_eval = True

	if derivative_eval == True:	
		print('')
		print('working tape set')
		set_working_tape(tape2)

	if TAPobject_data.data_name != None:	
		TAPobject_data.experimental_data = read_experimental_data_object(TAPobject_data.data_name)
		output_data = read_experimental_data_object(TAPobject_data.data_name)

	zone_keys = list(TAPobject_data.species.reference_diffusions.keys())

	outlet_diffusions = {}
	for k in TAPobject_data.species.gasses.keys():
		outlet_diffusions[k] = TAPobject_data.species.reference_diffusions[zone_keys[-1]]*(mp.sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/mp.sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.gasses[k].mass))
	for k in TAPobject_data.species.inert_gasses.keys():
		outlet_diffusions[k] = TAPobject_data.species.reference_diffusions[zone_keys[-1]]*(mp.sqrt(TAPobject_data.species.reference_mass*TAPobject_data.species.temperature)/mp.sqrt(TAPobject_data.species.reference_temperature*TAPobject_data.species.inert_gasses[k].mass))
	
	species_time = []
	for j in TAPobject_data.species.gasses:
		species_time.append(TAPobject_data.species.gasses[j].delay)
	for j in TAPobject_data.species.inert_gasses:
		species_time.append(TAPobject_data.species.inert_gasses[j].delay)			
	
	TAPobject_data.species.reference_mass = Constant(TAPobject_data.species.reference_mass)
	if type(TAPobject_data.species.temperature) == dict:
		for j_temp_number, j_temps in enumerate(TAPobject_data.species.temperature): 
			TAPobject_data.species.temperature[j_temp_number] = Constant(TAPobject_data.species.temperature[j_temp_number])
	else:
		TAPobject_data.species.temperature = Constant(TAPobject_data.species.temperature)
	
	TAPobject_data.species.reference_temperature = Constant(TAPobject_data.species.reference_temperature)
	## EDIT
	#TAPobject_data.reactor.voids[0] = Constant(TAPobject_data.reactor.voids[0])
	TAPobject_data.species.catalyst_diffusion = Constant(TAPobject_data.species.catalyst_diffusion)
	TAPobject_data.species.inert_diffusion = Constant(TAPobject_data.species.inert_diffusion)
	#TAPobject_data.reactor.voids[1] = Constant(TAPobject_data.reactor.voids[1])
	standard_parameters['kbt'] = Constant(standard_parameters['kbt'])
	standard_parameters['h'] = Constant(standard_parameters['h'])
	standard_parameters['Rgas'] = Constant(standard_parameters['Rgas'])
	#TAPobject_data.reactor.radius = Constant(TAPobject_data.reactor.radius)
			
	species_order_dictionary = {}
	total_species = 0
	for k in TAPobject_data.species.gasses:
		species_order_dictionary[k] = total_species
		TAPobject_data.species.gasses[k].mass = Constant(TAPobject_data.species.gasses[k].mass)
		TAPobject_data.species.gasses[k].delay = Constant(TAPobject_data.species.gasses[k].delay)
		TAPobject_data.species.gasses[k].intensity = Constant(TAPobject_data.species.gasses[k].intensity)
		total_species += 1
	for k in TAPobject_data.species.adspecies:
		species_order_dictionary[k] = total_species
		TAPobject_data.species.adspecies[k].conc = Constant(TAPobject_data.species.adspecies[k].conc)
		total_species += 1
	for k in TAPobject_data.species.inert_gasses:
		species_order_dictionary[k] = total_species
		TAPobject_data.species.inert_gasses[k].mass = Constant(TAPobject_data.species.inert_gasses[k].mass)
		TAPobject_data.species.inert_gasses[k].delay = Constant(TAPobject_data.species.inert_gasses[k].delay)
		TAPobject_data.species.inert_gasses[k].intensity = Constant(TAPobject_data.species.inert_gasses[k].intensity)	
		total_species += 1
	
	if TAPobject_data.mechanism.reactants != []:
		for k,z in enumerate(TAPobject_data.mechanism.matrix):
			if TAPobject_data.mechanism.processes[k].f.use == 'G':
				TAPobject_data.mechanism.processes[k].f.Ga = Constant(TAPobject_data.mechanism.processes[k].f.Ga)
				TAPobject_data.mechanism.processes[k].f.dG = Constant(TAPobject_data.mechanism.processes[k].f.dG)
			elif TAPobject_data.mechanism.processes[k].f.use == 'E':
				TAPobject_data.mechanism.processes[k].f.Ao = Constant(TAPobject_data.mechanism.processes[k].f.Ao)
				TAPobject_data.mechanism.processes[k].f.Ea = Constant(TAPobject_data.mechanism.processes[k].f.Ea)
			elif TAPobject_data.mechanism.processes[k].f.use == 'k':
				TAPobject_data.mechanism.processes[k].f.k = Constant(TAPobject_data.mechanism.processes[k].f.k)		
			if TAPobject_data.mechanism.processes[k].b.use == 'E':
				TAPobject_data.mechanism.processes[k].b.Ao = Constant(TAPobject_data.mechanism.processes[k].b.Ao)
				TAPobject_data.mechanism.processes[k].b.Ea = Constant(TAPobject_data.mechanism.processes[k].b.Ea)
			elif TAPobject_data.mechanism.processes[k].b.use == 'k':
				try:
					TAPobject_data.mechanism.processes[k].b.k = Constant(TAPobject_data.mechanism.processes[k].b.k)
				except:
					pass
		if TAPobject_data.mechanism.links != {}:
			for j in TAPobject_data.mechanism.links:
				TAPobject_data.mechanism.links[j] = Constant(TAPobject_data.mechanism.links[j])

	time_steps = pulse_time*1000

	dk = Constant(pulse_time/time_steps)
	#cat_location = 1 - TAPobject_data.reactor.cat_center

	dx_r = TAPobject_data.reactor.length/TAPobject_data.mesh
	point_volume = dx_r*TAPobject_data.reactor.cs_area*TAPobject_data.reactor.voids[0]
	reference_pulse_concentration = TAPobject_data.species.reference_pulse_size/point_volume
	controls = []
	low_bounds = []
	up_bounds = []

	if derivative_eval == True:
		try:
			if len(TAPobject_data.poi) > 1:
				pass
		except:
			TAPobject_data.poi = TAPobject_data.poi[0]
		for j in TAPobject_data.poi:
			controls.append(Control(eval('TAPobject_data.mechanism.processes'+j)))
			
			w1 = j.split('[')[1]
			w2 = int(w1.split(']')[0])

			x1 = j.split('].')[1]
			x2 = x1.split('.')[0]

			if x2 == 'f':
				if 'dG' in j:
					low_bounds.append(-np.inf)
				elif TAPobject_data.mechanism.processes[w2].f.lower_bound != None:
					low_bounds.append(TAPobject_data.mechanism.processes[w2].f.lower_bound)
				else:
					low_bounds.append(0)
			
				if 'dG' in j:
					up_bounds.append(np.inf)
				elif TAPobject_data.mechanism.processes[w2].f.upper_bound != None:
					up_bounds.append(TAPobject_data.mechanism.processes[w2].f.upper_bound)
				else:
					up_bounds.append(np.inf)

			if x2 == 'b':
				if 'dG' in j:
					low_bounds.append(-np.inf)
				elif TAPobject_data.mechanism.processes[w2].b.lower_bound != None:
					low_bounds.append(TAPobject_data.mechanism.processes[w2].b.lower_bound)
				else:
					low_bounds.append(0)
			
				if 'dG' in j:
					up_bounds.append(np.inf)			
				elif TAPobject_data.mechanism.processes[w2].b.upper_bound != None:
					up_bounds.append(TAPobject_data.mechanism.processes[w2].b.upper_bound)
				else:
					up_bounds.append(np.inf)

	### Mesh automation

	mesh = UnitIntervalMesh(int(TAPobject_data.mesh))

	alternative = 1

	fraction_variation = 0
	
	for jk in TAPobject_data.reactor.fractions.keys():
		
		percentChange = 1000
		fraction_variation += TAPobject_data.reactor.fractions[jk]
		if fraction_variation < 0.995:
			while percentChange > 0.1:

				x_mesh = np.sort(mesh.coordinates(),axis=0)#.tolist()
				actual = x_mesh.flat[np.abs(x_mesh - fraction_variation).argmin()]
				percentChange = abs(round(100*((fraction_variation - actual)/fraction_variation),2))
				boundary_location = np.where(x_mesh == actual)[0]
				
				if percentChange > 0.1:
					cfDict = {}
					boundary_location = np.where(x_mesh == actual)[0]
					class thin_zoneTest(SubDomain):
						def inside(self, x, on_boundary):
							return between(x[0], ((x_mesh[boundary_location-2][0][0]), (x_mesh[boundary_location+2][0][0])))
			
					thin_zoneTest = thin_zoneTest()
					cfDict[0] = MeshFunction("bool",mesh,1)
			
					thin_zoneTest.mark(cfDict[0],1)
					mesh = refine(mesh,cfDict[0],True)

	cfDict = {}
	if type(TAPobject_data.refine) != {}:
		for n in TAPobject_data.refine.keys():
			class thin_zoneTest(SubDomain):
				def inside(self, x, on_boundary):
					return between(x[0], (TAPobject_data.refine[n][0][0], TAPobject_data.refine[n][0][1]))
			
			thin_zoneTest = thin_zoneTest()
			
			for jk in range(0,TAPobject_data.refine[n][1]):
				cfDict[1] = MeshFunction("bool",mesh,1)
				thin_zoneTest.mark(cfDict[1],1)
				mesh = refine(mesh,cfDict[1],True)

	if TAPobject_data.view_mesh == True:
		from matplotlib import pyplot
		plot(mesh)
		pyplot.show()

	P1 = FiniteElement('P',mesh.ufl_cell(),2)
	#P1 = FiniteElement('P',mesh.ufl_cell(),2)

	elements = []
	for k in range(0,len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+len(TAPobject_data.species.inert_gasses)):
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

	for knum,k in enumerate(list(TAPobject_data.species.gasses.keys())):
		v_d['v_'+k] = tempA[knum]
		u_d['u_'+k] = tempB[knum]
		u_nd['u_n'+k] = tempC[knum]
	for knum,k in enumerate(list(TAPobject_data.species.adspecies.keys())):
		v_d['v_'+k] = tempA[len(TAPobject_data.species.gasses)+knum]
		u_d['u_'+k] = tempB[len(TAPobject_data.species.gasses)+knum]
		u_nd['u_n'+k] = tempC[len(TAPobject_data.species.gasses)+knum]
	for knum,k in enumerate(list(TAPobject_data.species.inert_gasses.keys())):
		v_d['v_'+k] = tempA[len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+knum]
		u_d['u_'+k] = tempB[len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+knum]
		u_nd['u_n'+k] = tempC[len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+knum]
	
	if TAPobject_data.species.advection > 0:
			
		W = VectorFunctionSpace(mesh, 'P', 1)
		advTerm = Function(W)
		advMulti = Constant(TAPobject_data.species.advection)
		advValue = Constant(1)
		advTerm.vector()[:] = advValue

		advTerm.vector()[totalNumCells-0] = 0

	#class thin_zone(SubDomain):
	#	def inside(self, x, on_boundary):
	#		return between(x[0], ((Mesh12)-1/TAPobject_data.mesh, (Mesh22)+1/TAPobject_data.mesh))
	#			
	#thin_zone = thin_zone()
	#domains = MeshFunction("size_t", mesh,1)
	#thin_zone.mark(domains,1)
	
	#newBoundaries = domains.mesh().coordinates().transpose().tolist()[0]
	#additionalCells = len(newBoundaries[int(TAPobject_data.mesh+1):])
	
	##dx = Measure("dx",subdomain_data=domains)
	##dT = Measure("dx",subdomain_data=boundary_parts0)
		
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
	
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)
		
	dz = 1/TAPobject_data.mesh
	
	#class integration_section(SubDomain):
	#	def inside(self, x, on_boundary):
	#		return between(x[0], (1-dz,1.0))
	
	right = CompiledSubDomain("near(x[0], 1.)")
	boundary_parts = MeshFunction("size_t", mesh,  mesh.topology().dim()-1)
	right.mark(boundary_parts, 1)

	bcs = []
	for k in range(0,len(TAPobject_data.species.gasses)):
		bcs.append(DirichletBC(V.sub(k),Constant(0),boundary_R))
	
	for k in range(0,len(TAPobject_data.species.inert_gasses)):
		bcs.append(DirichletBC(V.sub(len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+k),Constant(0),boundary_R))
	
	if type(TAPobject_data.species.temperature) != dict:
		F_new = new_construct_f_equation(TAPobject_data)
		#F_new = construct_f_equation(TAPobject_data)
	else:
		F_new = construct_f_equation_multiple_experiments(TAPobject_data)

	constantT = Constant(0)
	remove_surface = Constant(1)
	
	if TAPobject_data.mesh < 600:
		b0Test1 = Expression('x[0] < 0.002500001 ? 0.5 : 0', degree=0)
	elif TAPobject_data.mesh < 1200:
		b0Test1 = Expression('x[0] < 0.0012500005 ? 0.5 : 0', degree=0)
	elif TAPobject_data.mesh < 1800:
		b0Test1 = Expression('x[0] < 0.00062500025 ? 0.5 : 0', degree=0)
	elif TAPobject_data.mesh < 2400:
		b0Test1 = Expression('x[0] < 0.00031250012 ? 0.5 : 0', degree=0)
		
	Fpulses = ''
	for knum,k in enumerate(TAPobject_data.species.gasses):
		Fpulses += " -TAPobject_data.species.gasses['"+k+"'].intensity*reference_pulse_concentration*b0Test1*exp(-(constantT - round(TAPobject_data.species.gasses['"+k+"'].delay,6))*(constantT - round(TAPobject_data.species.gasses['"+k+"'].delay,6))/(0.00000000001))*v_d['v_"+k+"']*dx"
		if knum < len(TAPobject_data.species.gasses)-1:
			Fpulses += ' + '
	for knum,k in enumerate(TAPobject_data.species.inert_gasses):
		Fpulses += "-TAPobject_data.species.inert_gasses['"+k+"'].intensity*reference_pulse_concentration*b0Test1*exp(-(constantT - round(TAPobject_data.species.inert_gasses['"+k+"'].delay,6))*(constantT - round(TAPobject_data.species.inert_gasses['"+k+"'].delay,6))/(0.00000000001))*v_d['v_"+k+"']*dx"
		if knum < len(TAPobject_data.species.inert_gasses)-1:
			Fpulses += ' + '


	b0Test2 = {}
	include_surface = Constant(1)
	for knum,k in enumerate(TAPobject_data.species.adspecies):
		
		if TAPobject_data.species.adspecies[k].zones == []:
			
			if TAPobject_data.species.adspecies[k].profile == 'uniform':
				b0Test2[k] = Expression('1', degree=2)
			elif TAPobject_data.species.adspecies[k].profile == 'linear':
				b0Test2[k] = Expression('A*x[0] + b',A=TAPobject_data.species.adspecies[k].p_parms['a'],b=TAPobject_data.species.adspecies[k].p_parms['b'], degree=2)
			elif TAPobject_data.species.adspecies[k].profile == 'bounded':
				eq_construct = ''
				for jk in TAPobject_data.species.adspecies[k].p_parms:
					if eq_construct != '':
						eq_construct+= ' ^ '
					eq_construct+=str(TAPobject_data.species.adspecies[k].p_parms[jk][0])
					eq_construct+=' <= x[0] && x[0] < '
					eq_construct+=str(TAPobject_data.species.adspecies[k].p_parms[jk][1])
				eq_construct+=' ? 1 : 0'
				b0Test2[k] = Expression(eq_construct,degree=2)
			Fpulses += "-include_surface*TAPobject_data.species.adspecies['"+k+"'].conc*b0Test2['"+k+"']*exp(-(constantT)*(constantT)/(0.00000000001))*v_d['v_"+k+"']*dx"
			if knum < len(TAPobject_data.species.adspecies)-1:
				Fpulses += ' + '

		else:

			eq_construct = ''
			if float(TAPobject_data.species.adspecies[k].conc) != 0:
				for jk in TAPobject_data.species.adspecies[k].zones:
					
					if eq_construct != '':
						eq_construct+= ' ^ '
					if jk == 0:
						z_start = 0
						z_end   = TAPobject_data.reactor.fractions[0]
					else:
						z_start = 0
						for z in range(0,jk):
							z_start +=  TAPobject_data.reactor.fractions[z]
						z_end = z_start + TAPobject_data.reactor.fractions[jk]
					eq_construct+=str(z_start)
					eq_construct+=' <= x[0] && x[0] < '
					eq_construct+=str(z_end)
				eq_construct+=' ? 1 : 0'
				
				b0Test2[k] = Expression(eq_construct,degree=2)
				Fpulses += "-include_surface*TAPobject_data.species.adspecies['"+k+"'].conc*b0Test2['"+k+"']*exp(-(constantT)*(constantT)/(0.00000000001))*v_d['v_"+k+"']*dx"
				
				if knum < len(TAPobject_data.species.adspecies)-1:
					Fpulses += ' + '

	eq_construct = ''
	
	for jk in TAPobject_data.reactor.voids.keys():
		if eq_construct != '':
			eq_construct+= ' ^ '
		if jk == 0:
			z_start = 0
			z_end   = TAPobject_data.reactor.fractions[0]
		else:
			z_start = 0
			for z in range(0,jk):
				z_start +=  TAPobject_data.reactor.fractions[z]
			z_end = z_start + TAPobject_data.reactor.fractions[jk]
		eq_construct+=str(z_start)
		eq_construct+=' <= x[0] && x[0] < '
		eq_construct+=str(z_end)
		eq_construct+=' ? '+str(float(TAPobject_data.reactor.voids[jk]))+' : 0'
	reference_voids = Expression(eq_construct,degree=2)


	eq_construct = ''

	for jk in TAPobject_data.reactor.voids.keys():
		if eq_construct != '':
			eq_construct+= ' ^ '
		if jk == 0:
			z_start = 0
			z_end   = TAPobject_data.reactor.fractions[0]
		else:
			z_start = 0
			for z in range(0,jk):
				z_start +=  TAPobject_data.reactor.fractions[z]
			z_end = z_start + TAPobject_data.reactor.fractions[jk]
		eq_construct+=str(z_start)
		eq_construct+=' <= x[0] && x[0] < '
		eq_construct+=str(z_end)
		eq_construct+=' ? '+str(float(TAPobject_data.species.reference_diffusions[jk]))+' : 0'
	
	zone_diffusions = Expression(eq_construct,degree=2)


	####for knum,k in enumerate(TAPobject_data.species.adspecies):
	####	##TAPobject_data.species.adspecies[k].conc = Constant(TAPobject_data.species.adspecies[k].conc)
	####	Fpulses += "-remove_surface*TAPobject_data.species.adspecies['"+k+"'].conc*exp(-(constantT)*(constantT)/(0.00000000001))*v_d['v_"+k+"']*dx"
	####	if knum < len(TAPobject_data.species.adspecies)-1:
	####		Fpulses += ' + '
	
	F_new += Fpulses

	#!#! ADDING INVERSE ANALYSIS
	if TAPobject_data.tangent_linear_sensitivity == True:
		#TAPobject_data.species.gasses['CO'].intensity = Control(TAPobject_data.species.gasses['CO'].intensity)
		#TAPobject_data.mechanism.processes[0].f.k["value"] = Control(TAPobject_data.mechanism.processes[0].f.k["value"])
		#c = TAPobject_data.mechanism.processes[0].f.k["value"]
		#c.tlm_value = TAPobject_data.mechanism.processes[0].f.k["value"]
		sens_param = eval(TAPobject_data.poi[0])
		
		c = sens_param
		c.tlm_value = sens_param
		#c = TAPobject_data.species.gasses['CO'].intensity
		#c.tlm_value = TAPobject_data.species.gasses['CO'].intensity
	
		SV_du = FunctionSpace(mesh,P1)
		Sw_new = Expression('A',A=Constant(1),degree=0)
		Sw_new2 = interpolate(Sw_new,SV_du)
		Sw3 = project(Sw_new2,SV_du)
	
		sensFuncs = {}
		for k_gasses in TAPobject_data.species.gasses:
			sensFuncs[k_gasses] = []

	#if TAPobject_data.adjoint_sensitivitiy == True:
	#	#print(TAPobject_data.parameters_of_interest[0])
	#	controls = [TAPobject_data.species.gasses['CO'].intensity]
	#	# output_fitting = pointFitting(legend_label[:int(len(legend_label)-reac_input['Number of Inerts'])],reac_input['Time Steps'],reac_input['Experimental Data Folder'],reac_input['Pulse Duration'],reac_input['Objective Points'],objSpecies)


	theta = Constant(1)
	Ftemp = eval(F_new)
	theta = Constant(0.5)
	F = eval(F_new)

	J = derivative(F,u)
	Jtemp = derivative(Ftemp,u)
	
	variational_solver = 'newton'
	if variational_solver == 'constrained':
		snes_solver_parameters = {"nonlinear_solver": "snes","snes_solver": {"linear_solver": "lu","line_search":'basic',"maximum_iterations": 100,"report": False,"error_on_nonconvergence": False}}
			
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
		#problemtemp = NonlinearVariationalProblem(Ftemp,u,bcs,Jtemp)
		solvertemp = NonlinearVariationalSolver(problemtemp)

		solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
		solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-8
		solver.parameters["newton_solver"]["maximum_iterations"] = 1000

		solvertemp.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
		solvertemp.parameters["newton_solver"]["absolute_tolerance"] = 1e-8
		solvertemp.parameters["newton_solver"]["maximum_iterations"] = 1000

	synthetic_data = {}
	synthetic_data['time'] = {}
	if TAPobject_data.store_flux_data == True:
		for j in TAPobject_data.species.gasses:
			synthetic_data[j] =  {}
		for j in TAPobject_data.species.inert_gasses:
			synthetic_data[j] =  {}

	if TAPobject_data.store_thin_data == True:
		thin_data = {}
		thin_data['time'] = {}
		if TAPobject_data.store_thin_data == True:
			for j in TAPobject_data.species.gasses:
				thin_data[j] =  {}
			for j in TAPobject_data.species.adspecies:
				thin_data[j] =  {}
			for j in TAPobject_data.species.inert_gasses:
				thin_data[j] =  {}

	if TAPobject_data.tangent_linear_sensitivity == True:
		sensitivity_output = {}
			
		for k_sens in TAPobject_data.species.gasses:
			sensitivity_output[k_sens] = []

	##if derivative_eval == True:
	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (1-dz,1.0))
	#class thin_zone(SubDomain):
	#	def inside(self, x, on_boundary):
	#		return between(x[0], ((Mesh12)-1/TAPobject_data.mesh, (Mesh22)+1/TAPobject_data.mesh))
	#			
	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,1)
	domains.set_all(0)
	thin_zone.mark(domains,1)
	dx = Measure('dx',domain = mesh, subdomain_data=domains)

	#domains = MeshFunction("size_t", mesh,0)
	#
	#osub = integration_section()
	#domains = MeshFunction("size_t", mesh,0)
	#domains.set_all(0)
	#osub.mark(domains, 1)
	#dP = Measure('vertex',domain = mesh, subdomain_data=domains)

	try: 
		os.mkdir('./'+TAPobject_data.output_name)
	except OSError:  
		pass
	
	#output_data = curveFitting(pulse_time,TAPobject_data)
	#TAPobject_data.species.gasses['CO'].sigma = 0.5

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
			for jz_num,jz in enumerate(TAPobject_data.poi):
				print(jz)
				print(jz)
				name_columns.append(jz)
				value_row.append(mv[jz_num])
			new_addition.loc[len(new_addition)] = value_row
			new_addition.to_csv('./'+TAPobject_data.output_name+'/optimization_results.csv',index=False)
		
		except:
			name_columns = ['objective']
			value_row = [j]
			for jz_num,jz in enumerate(TAPobject_data.poi):
				name_columns.append(jz)
				value_row.append(mv[jz_num])
			percentile_list = pd.DataFrame([value_row],columns=name_columns)
			percentile_list.to_csv('./'+TAPobject_data.output_name+'/optimization_results.csv',index=False)
	
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

	time_ave = []
	tz_max = []
	tz_ave = []
	tz_min = []

	for k_pulse in range(0,pulse_number):
		print('Pulse #: '+str(k_pulse))
		t = 0

		if k_pulse > 0:
			constantT.assign(0.0)
			remove_surface.assign(0.0)
		dt = pulse_time/time_steps
		start_time = time.time()
		
		
		if TAPobject_data.store_flux_data == True:
			synthetic_data['time'][k_pulse] =  []
			for j in TAPobject_data.species.gasses:
				synthetic_data[j][k_pulse] =  []
			for j in TAPobject_data.species.inert_gasses:
				synthetic_data[j][k_pulse] =  []

		if TAPobject_data.store_thin_data == True:
			thin_data['time'][k_pulse] =  []
			for j in TAPobject_data.species.gasses:
				thin_data[j][k_pulse] =  []
			for j in TAPobject_data.species.adspecies:
				thin_data[j][k_pulse] =  []
			for j in TAPobject_data.species.inert_gasses:
				thin_data[j][k_pulse] =  []

		all_species = len(TAPobject_data.species.gasses) + len(TAPobject_data.species.inert_gasses) + len(TAPobject_data.species.adspecies)
		time_step = 0
		step_number = 0
		
		#if TAPobject_data.conc_profiles != []:
		#	fig, ax = plt.subplots()ø

		#fraction_variation = 0
		#if TAPobject_data.conc_profiles != []:
		#	temp_zone_specs = {}
		#	x_mesh = np.sort(mesh.coordinates(),axis=0)#.tolist()
		#	for jl in TAPobject_data.conc_profiles:
		#		for kl,kl_num in enumerate()
		#		for kl in range(0,int(jl)):
		#			fraction_variation += TAPobject_data.reactor.fractions[kl]
		#		actual = x_mesh.flat[np.abs(x_mesh - fraction_variation).argmin()]
		#		boundary_location_1 = np.where(x_mesh == actual)#[0]
		#		fraction_variation += TAPobject_data.reactor.fractions[int(jl)]
		#		actual = x_mesh.flat[np.abs(x_mesh - fraction_variation).argmin()]
		#		boundary_location_2 = np.where(x_mesh == actual)#[0]
		#		temp_zone_specs[jl] = [boundary_location_1[0][0],boundary_location_2[0][0]]
		#	print(temp_zone_specs)
		#else:
		#	print('')
		#	pass

		x_mesh = np.sort(TAPobject_data.reactor.length*mesh.coordinates(),axis=0)

		if TAPobject_data.conc_profiles != []:
			surface_store = {}
			surface_store['mesh'] = x_mesh
			np.savetxt('./'+TAPobject_data.output_name+'/mesh.csv', surface_store['mesh'], delimiter=",")
			sys.exit()

		if TAPobject_data.conc_profiles != []:
			data_store = {}
			for jk in TAPobject_data.mechanism.reactants:
				data_store[jk] = np.empty([int(time_steps)+1,x_mesh.shape[0]])

		if TAPobject_data.objective_return == True:
			objective_value = 0
		while t <= pulse_time:
			
			if step_number%TAPobject_data.data_storage_frequency == 0:
				if TAPobject_data.store_flux_data == True:
					synthetic_data['time'][k_pulse].append(round(t,6))
					for knum,k in enumerate(TAPobject_data.species.gasses):
						synthetic_data[k][k_pulse].append(2*(float(outlet_diffusions[k]) /(float(dx_r))) * (float(TAPobject_data.reactor.radius)**2)*3.14159*( float(u_n.vector().get_local()[(all_species)+knum])))
						##synthetic_data[k][k_pulse].append(2*(float(TAPobject_data.species.gasses[k].inert_diffusion) /(float(dx_r))) * (float(TAPobject_data.reactor.radius)**2)*3.14159*( float(u_n.vector().get_local()[(all_species)+knum])))
					for knum,k in enumerate(TAPobject_data.species.inert_gasses):
						synthetic_data[k][k_pulse].append(2*(float(outlet_diffusions[k]) /(float(dx_r))) * (float(TAPobject_data.reactor.radius)**2)*3.14159*( float(u_n.vector().get_local()[(all_species)+len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+knum])))
						##synthetic_data[k][k_pulse].append(2*(float(TAPobject_data.species.inert_gasses[k].inert_diffusion) /(float(dx_r))) * (float(TAPobject_data.reactor.radius)**2)*3.14159*( float(u_n.vector().get_local()[all_species+len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+knum])))

			if step_number%TAPobject_data.data_storage_frequency == 0:
				if TAPobject_data.conc_profiles != []:
					#thin_data['time'][k_pulse].append(round(t,6))
					for jk in TAPobject_data.conc_profiles:
						try:
							data_store[jk][time_step] = np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2])
						except:
							pass
					#	thin_data[k][k_pulse].append(float(u_n.vector().get_local()[(catalyst_center_cell*(len(TAPobject_data.species.gasses) + len(TAPobject_data.species.adspecies) + len(TAPobject_data.species.inert_gasses)))+all_species+knum]))
					#for knum,k in enumerate(TAPobject_data.species.adspecies):
					#	thin_data[k][k_pulse].append(float(u_n.vector().get_local()[(catalyst_center_cell*(len(TAPobject_data.species.gasses) + len(TAPobject_data.species.adspecies) + len(TAPobject_data.species.inert_gasses)))+all_species+len(TAPobject_data.species.gasses)+knum]))
					#for knum,k in enumerate(TAPobject_data.species.inert_gasses):
					#	thin_data[k][k_pulse].append(float(u_n.vector().get_local()[(catalyst_center_cell*(len(TAPobject_data.species.gasses) + len(TAPobject_data.species.adspecies) + len(TAPobject_data.species.inert_gasses)))+all_species+len(TAPobject_data.species.gasses)+len(TAPobject_data.species.adspecies)+knum]))

			if TAPobject_data.optimize == True:
				if True == True:
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
								if k_fitting in TAPobject_data.species.gasses:
									jfunc_2 += assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3)*dx(1))/((TAPobject_data.species.gasses[k_fitting].sigma)**2)								
									##jfunc_2 += assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.species.gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(TAPobject_data.species.gasses[k_fitting].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.species.gasses[k_fitting].sigma)**2)
								elif k_fitting in TAPobject_data.species.inert_gasses:
									jfunc_2 += assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.species.inert_gasses[k_fitting].sigma)**2)
									
							except UnboundLocalError:
								w_temp_2 = Expression('1',degree=0) 
								w_temp2_2 = interpolate(w_temp_2,V_du)
								w4_2 = project(w_temp2_2,V_du)	
								if k_fitting in TAPobject_data.species.gasses:
									jfunc_2 = assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3)*dx(1))/((TAPobject_data.species.gasses[k_fitting].sigma)**2)	
								elif k_fitting in TAPobject_data.species.inert_gasses:
									jfunc_2 = assemble(inner(u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3,u_n[species_order_dictionary[k_fitting]]*(2*(outlet_diffusions[k_fitting] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159) - w3)*dP(1))/((TAPobject_data.species.inert_gasses[k_fitting].sigma)**2)	
								

					if TAPobject_data.thermodynamic_constraints == True and t == 0:
						thermoReactions = {}
						thermoStoich = {}
						for znum,z in enumerate(TAPobject_data.thermo_equations):
							thermoReactions[znum] = []
							thermoStoich[znum] = []
							tempValue = z.split('+')
							for j in tempValue:
								if j.find('*') != -1:
									thermoSplit = j.split('*')
									thermoSplit[1] = thermoSplit[1].replace('r','')
									thermoReactions[znum].append(int(thermoSplit[1]))
									thermoSplit[0] = thermoSplit[0].replace('(','')
									thermoSplit[0] = thermoSplit[0].replace(')','')
									thermoStoich[znum].append(float(thermoSplit[0]))
								else:
									j = j.replace('r','')
									thermoReactions[znum].append(int(j))
									thermoStoich[znum].append(1)

						w_temp = {}
						w_temp2 = {}
						w4 = {}
						tempFunc = {}
					
						def step_stoichiometry(kinetic_step):
							z = TAPobject_data.mechanism.matrix[kinetic_step]
							neg = []
							val_neg = []
							pos = []
							val_pos = []
							for j,v in enumerate(z):
								if v < 0:
									neg.append(TAPobject_data.mechanism.reactants[j])
									val_neg.append(v)
								elif v > 0:
									pos.append(TAPobject_data.mechanism.reactants[j])
									val_pos.append(v)
							return val_neg, val_pos

						for kip in thermoReactions.keys():

							w_temp[kip] = Expression(str(TAPobject_data.thermodynamic_free_energy[kip]),degree=0) # deltaG = sum(-R*T*ln(kf/kb))
							w_temp2[kip] = interpolate(w_temp[kip],V_du)
							w4[kip] = project(w_temp2[kip],V_du)
							thermoWeight = TAPobject_data.thermodynamic_alpha[kip]

							if TAPobject_data.mechanism.processes[0].f.use == 'k':
								for jnum,jval in enumerate(thermoReactions[kip]):
									# TAPobject_data.parameter_scale**("+str(scale_magnitude)+"
									new_neg, new_pos = step_stoichiometry(jval)
									scale_magnitude_1 = sum(new_neg)-1 
									scale_magnitude_2 = sum(new_pos)-1
									if jnum == 0:
										tempFunc[kip] = thermoStoich[kip][jnum]*(-0.008314*TAPobject_data.species.temperature)*ln( (TAPobject_data.parameter_scale**(scale_magnitude_1 - scale_magnitude_2))*(TAPobject_data.mechanism.processes[jval-1].f.k/TAPobject_data.mechanism.processes[jval-1].b.k))
									else:
										tempFunc[kip] += thermoStoich[kip][jnum]*(-0.008314*TAPobject_data.species.temperature)*ln((TAPobject_data.parameter_scale**(scale_magnitude_1 - scale_magnitude_2))*(TAPobject_data.mechanism.processes[jval-1].f.k/TAPobject_data.mechanism.processes[jval-1].b.k))

							if TAPobject_data.mechanism.processes[0].f.use == 'G':
								for jnum,jval in enumerate(thermoReactions[kip]):
									if jnum == 0:
										tempFunc[kip] = thermoStoich[kip][jnum]*mechanism.processes[jval-1].f.dG
									else:
										tempFunc[kip] += thermoStoich[kip][jnum]*r_const["dG"+str(jval-1)]

							tempFunc[kip] = (tempFunc[kip]-w4[kip])*thermoWeight
							jfunc_2 += assemble(inner(tempFunc[kip],tempFunc[kip])*dx())
							print('Simulated Thermo Value')
							print(jfunc_2)

			if round(t,6) not in species_time:
				try:
					if derivative_eval == True:
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
							solver.solve(annotate=False)
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
					if derivative_eval == True:
						solvertemp.solve()
					else:
						solvertemp.solve(annotate = False)
		
				except RuntimeError:
					print('Time Step Failure')
					sys.exit()
				

			#!#! ADDING INVERSE ANALYSIS
			if TAPobject_data.tangent_linear_sensitivity == True:
				for knum,k in enumerate(TAPobject_data.species.gasses):
					
					new_val = 2*(outlet_diffusions[k] /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159
					##new_val = 2*(TAPobject_data.species.gasses[k].inert_diffusion /(dx_r)) * (TAPobject_data.reactor.radius**2)*3.14159
					if k == 'CO':
						test_tlm_evaluation.append(assemble(inner((new_val*(u_n[knum])), Sw3)*dP(1) ))
					sensFuncs[k].append(assemble(inner((new_val*(u_n[knum])), Sw3)*dP(1) ))
			
			progressBar(t, pulse_time)

			#if TAPobject_data.plot_con_profile != []:
			#if (TAPobject_data.profile_time[0] <= t) and (t <= TAPobject_data.profile_time[1]):
			##for kj_num, kj in  enumerate(TAPobject_data.species.adspecies.keys()):
			##ax.plot(x_mesh, np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index('H2')].vector().get_local()[::2]),color='r')#,label=species_label[kj_num])
			if True == False:
				for jk in TAPobject_data.conc_profiles:
					##ax.plot(x_mesh[97:127], np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2][97:127]),color='r')#,label=species_label[kj_num])
					time_ave.append(t+k_pulse*0.8)
					tz_max.append(np.max(np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2][97:127])))
					#tz_ave.append(np.average(np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2][97:127])))
					
					temp2 = np.add(np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index('C2H5*(0)')].vector().get_local()[::2][97:127]),np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index('CH3*(0)')].vector().get_local()[::2][97:127]))
					temp3 = np.add(temp2,np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index('CH2*(0)')].vector().get_local()[::2][97:127]))
					temp4 = np.add(temp3,np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index('H*(0)')].vector().get_local()[::2][97:127]))
					
					tz_ave.append(np.average(temp4))
					#tz_ave.append(np.average(np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2][97:127])))

					tz_min.append(np.min(np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2][97:127])))
				#	data_store[jk][time_step] = np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2])
				#for jl in temp_zone_specs.keys():
				#	data_store[jk][time_step] = np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2][temp_zone_specs[jl][0]:temp_zone_specs[jl][1]])
				
			#	for jk_num, jk in enumerate(TAPobject_data.mechanism.reactants):
			#		#print(jk_num)
			#		#print(jk)
			#		data_store[jk][time_step] = np.flip(u_n.split(deepcopy=True)[TAPobject_data.mechanism.reactants.index(jk)].vector().get_local()[::2])
			#		#if jk == '*(0)':
			#			#print(specs_copy.mechanism.reactants.index(jk))  
			#			#print(np.flip(u_n.split(deepcopy=True)[specs_copy.mechanism.reactants.index(jk)].vector().get_local()[::2]))
			
			Uvector = as_backend_type(u.vector()).get_local()
			Uvector[Uvector <= DOLFIN_EPS] = DOLFIN_EPS
			u.vector().set_local(Uvector)
			
			u_n.assign(u)
				
			t += dt
			constantT.assign(round(t,6))
			time_step += 1
			step_number += 1
			include_surface.assign(0.0)# = Constant(0)
			if TAPobject_data.conc_profiles != []:
				for jk in TAPobject_data.conc_profiles:
					np.savetxt('./'+TAPobject_data.output_name+'/'+jk+'_pulse_'+str(k_pulse)+'.csv', data_store[jk], delimiter=",")

			if TAPobject_data.objective_return == True:
				objective_value += jfunc_2
		
		print(processTime(start_time))

		#if TAPobject_data.conc_profiles != []:
		#	plt.show()

	##fig,ax = plt.subplots()
	##ax.set_ylabel('Thin-zone concentration (nmol/cm3)')
	##ax.set_xlabel('Time (s)')
	###average = sum(tz_ave[:300])/len(tz_ave[:300])
	###plt.axhline(y=average,color='k')
	##ax.plot(time_ave,tz_ave,color='k')
	###ax.plot(time_ave,tz_max,linestyle='--',color='k')
	###ax.plot(time_ave,tz_min,linestyle='--',color='k')
	##plt.show()

	if TAPobject_data.objective_return == True:
		print('test')
		return objective_value

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

		#self.lower_bound = None
		#self.upper_bound = None

		#for gt in range(0,len(controls)):
		#	low_bounds.append(0)
		#	up_bounds.append(np.inf)
		print('Lower parameter bounds are:')
		print(low_bounds)
		print('Upper parameter bounds are:')
		print(up_bounds)
		if TAPobject_data.optimization_method == 'L-BFGS-B':
			##u_opt_2 = minimize(rf_2)
			u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds), tol=TAPobject_data.ftol, options={"ftol":TAPobject_data.ftol,"gtol":TAPobject_data.gtol})
			#u_opt_2 = minimize(rf_2, bounds = (low_bounds,up_bounds), tol=1e-22, options={"ftol":1e-22,"gtol":1e-22})
		else:
			u_opt_2 = minimize(rf_2, method = TAPobject_data.optimization_method, tol=TAPobject_data.ftol, options={"ftol":TAPobject_data.ftol,"gtol":TAPobject_data.gtol})
			#u_opt_2 = minimize(rf_2, method = TAPobject_data.optimization_method, tol=1e-22, options={"ftol":1e-22,"gtol":1e-22})

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
							synthetic_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.species.gasses[k].noise# +beta_2*np.cos(w_2*(k*dt))
							
						except:
							try:
								synthetic_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.species.inert_gasses[k].noise# +beta_2*np.cos(w_2*(k*dt))
							except:
								pass
	# Edit more
	if TAPobject_data.store_thin_data == True and TAPobject_data.surface_noise == True:
		beta_2 = 0.00270
		w_2 = 2*3.14159*70
		for k in thin_data:
			if k != 'time':
				for j in thin_data[k]:
					for z in range(0,len(thin_data[k][j])):
						if k in TAPobject_data.species.gasses:
							try:
								thin_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.species.gasses[k].noise# +beta_2*np.cos(w_2*(k*dt))
							except:
								pass
						if k in TAPobject_data.species.adspecies:
							try:
								thin_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.species.adspecies[k].noise# +beta_2*np.cos(w_2*(k*dt))
							except:
								pass
						if k in TAPobject_data.species.inert_gasses:						
							try:
								thin_data[k][j][z] += np.random.normal(0,1)*TAPobject_data.species.inert_gasses[k].noise# +beta_2*np.cos(w_2*(k*dt))
							except:
								pass
						

	if TAPobject_data.store_flux_data == True:
		save_object(synthetic_data,'./'+TAPobject_data.output_name+'/TAP_experimental_data.json')
		new_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_experimental_data.json')

	if TAPobject_data.store_thin_data == True:
		save_object(thin_data,'./'+TAPobject_data.output_name+'/TAP_thin_data.json')
		new_data = read_experimental_data_object('./'+TAPobject_data.output_name+'/TAP_thin_data.json')		
	
	if TAPobject_data.show_graph == True:
		fig, ax = plt.subplots()
		ax.set_ylabel('Flow (nmol/s)')
		ax.set_xlabel('Time (s)')
		for j in TAPobject_data.species.gasses:
			for k in synthetic_data[j]:
				if k == 0:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k],label=j)
				else:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k])

		for j in TAPobject_data.species.inert_gasses:
			for k in synthetic_data[j]:
				if k == 0:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k],label=j)		
				else:
					plt.plot(synthetic_data['time'][0], synthetic_data[j][k])
		plt.legend()
		plt.show()