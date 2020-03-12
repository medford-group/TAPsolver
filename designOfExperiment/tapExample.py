from fenics import *
from fenics_adjoint import *
import mpmath
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import csv
import ast
import shutil
import imageio
import sys
import os
import scipy
import pip
import pkg_resources
import ufl
from ufl import sqrt,exp,ln
import warnings

kwargs_write = {'fps':4.0, 'quantizer':'nq'}
figList = []


set_log_level(30)
############# Unit Interval mesh #################
def callSim(currentTimes,pulseIntensities,initialCompositions):
	tape2 = Tape()
	tape2.clear_tape()
	set_working_tape(tape2)
	
	mesh_size = 400
	mesh = UnitIntervalMesh(mesh_size)
	
	
	catalystRefinement = 0
	cfDict = {}
	
	for jayz in range(0,catalystRefinement+1):
		class thin_zoneTest(SubDomain):
			def inside(self, x, on_boundary):
				return between(x[0], (((1-0.5) - 0.5*0.01), ((1-0.5) + 0.5*0.01))) # between(x[0], (((1-cat_location) - 0.5*frac_length), ((1-cat_location) + 0.5*frac_length)))
	
		thin_zoneTest = thin_zoneTest()
		cfDict[jayz] = MeshFunction("bool",mesh,1)
		
		thin_zoneTest.mark(cfDict[jayz],jayz)
		mesh = refine(mesh,cfDict[jayz],True)
	
	
	####### For multiple species ###########
	
	CG1 = FiniteElement('CG',mesh.ufl_cell(),1) # Continuous galerkin
	P1 = FiniteElement('P',mesh.ufl_cell(),1)
	combElements = [P1,P1,P1,P1,P1,P1]
	elements = MixedElement(combElements)
	V = FunctionSpace(mesh,elements)
	V_du = FunctionSpace(mesh,P1)


	u = Function(V)
	u_n = Function(V)
	
	tempA = TestFunctions(V)
	tempB = split(u)
	tempC = split(u_n)
	
	v_1 = tempA[0]
	v_2 = tempA[1]
	v_3 = tempA[2]
	v_4 = tempA[3]
	v_5 = tempA[4]
	v_6 = tempA[5]
	
	u_1 = tempB[0]
	u_2 = tempB[1]
	u_3 = tempB[2]
	u_4 = tempB[3]
	u_5 = tempB[4]
	u_6 = tempB[5]
	
	un_1 = tempC[0]
	un_2 = tempC[1]
	un_3 = tempC[2]
	un_4 = tempC[3]
	un_5 = tempC[4]
	un_6 = tempC[5]
	
	# f = Expression('x[0]>=0.4 && x[0]<=.6 ? 1 : 0', degree=1)
	
	cat_location = 0.5
	frac_length = 0.01
	
	## First, define the boundaries at the entrance of the reactor
	tol = 1e-14
	
	# Boundary at the entrance
	def boundary_L(x, on_boundary):
		return on_boundary and near(x[0],0,tol)
	
	# Boundary at exit
	def boundary_R(x, on_boundary):
		return on_boundary and near(x[0],1,tol)
	
	
	class thin_zone(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (((1-0.5) - 0.5*0.01) - 0.0025, ((1-0.5) + 0.5*0.01) + 0.0025))
	
	
	BC1 = DirichletBC(V.sub(0),Constant(0),boundary_R)
	BC2 = DirichletBC(V.sub(1),Constant(0),boundary_R)
	BC3 = DirichletBC(V.sub(2),Constant(0),boundary_R)

	bcs = [BC1,BC2,BC3]
	
	dz = 1/mesh_size
	
	class integration_section(SubDomain):
		#def inside(self, x, on_boundary):
		#	return between(x[0], (1-dz,1.0))
		def inside(self, x, on_boundary):
			return between(x[0], (((1-0.5) - 0.5*0.01) - 0.0025, ((1-0.5) + 0.5*0.01) + 0.0025))
	
	osub = integration_section()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	osub.mark(domains, 1)
	dP = Measure('vertex',domain = mesh, subdomain_data=domains)
	
	
	class entrance_section(SubDomain):
		def inside(self, x, on_boundary):
			return between(x[0], (0,dz))
		
	esub = entrance_section()
	domains = MeshFunction("size_t", mesh,0)
	domains.set_all(0)
	esub.mark(domains, 1)
	dE = Measure('vertex',domain = mesh, subdomain_data=domains)
	
	
	thin_zone = thin_zone()
	domains = MeshFunction("size_t", mesh,1)
	thin_zone.mark(domains,1)
	#dT = Measure('dx',domain = mesh, subdomain_data=domains)
	
	k1f = Constant(2.24)
	k2f = Constant(0.4)
	k3f = Constant(10.0385)
	k4f = Constant(20.157)
	k1b = Constant(40.49)
	k2b = Constant(0.0)
	k3b = Constant(0.006297)
	k4b = Constant(0.005289207)
	
	
	#k1f = Constant(0)
	#k2f = Constant(0)
	#k3f = Constant(0)
	#k4f = Constant(0)
	#
	#k1b = Constant(0)
	#k2b = Constant(0)
	#k3b = Constant(0)
	#k4b = Constant(0)
	
	
	cntrl1 = Control(k1f)
	cntrl2 = Control(k2f)
	
	t = 0
	T = 0.2
	dt = 0.001
	dk = Constant(dt)
	constantT = Constant(t)
	
	class Delta(Expression):
		def __init__(self, eps):
			self.eps = eps
		def eval(self, values, x):
			eps = self.eps
			values[0] = eps/pi/(x[0]**2 + eps**2)

	eps = 1E-4
	
	toFlux = [1155.8851453677419, 1081.2315481068708, 922.077362410682, 967.0828963937186]
	
	# dk*(Dout["+str(k)+"]/(eb[0]*np.sum(r_param)**2) )
	
	diff1 = Constant(66.37200365)
	diff2 = Constant(62.08532435)
	diff3 = Constant(52.94654251)
	diff4 = Constant(55.5308022)
	
	eb = Constant(0.4)
	
	length2 = Constant(33.50994)
	
	# First define the diffusion for all terms
	F = (u_1 - un_1)*v_1*dx + dk*(diff1/(eb*length2))*dot(grad(u_1),grad(v_1))*dx
	#F = (u_1 - un_1)*v_1*dx(0) + dk*(diff1/(eb*length2))*dot(grad(u_1),grad(v_1))*dx(0)
	#F += (u_1 - un_1)*v_1*dx(1) + dk*(diff1/(eb*length2))*dot(grad(u_1),grad(v_1))*dx(1)
	
	F += (u_2 - un_2)*v_2*dx + dk*(diff2/(eb*length2))*dot(grad(u_2),grad(v_2))*dx
	#F += (u_2 - un_2)*v_2*dx(0) + dk*(diff2/(eb*length2))*dot(grad(u_2),grad(v_2))*dx(0)
	#F += (u_2 - un_2)*v_2*dx(1) + dk*(diff2/(eb*length2))*dot(grad(u_2),grad(v_2))*dx(1)
	
	F += (u_3 - un_3)*v_3*dx + dk*(diff3/(eb*length2))*dot(grad(u_3),grad(v_3))*dx
	#F += (u_3 - un_3)*v_3*dx(0) + dk*(diff3/(eb*length2))*dot(grad(u_3),grad(v_3))*dx(0)
	#F += (u_3 - un_3)*v_3*dx(1) + dk*(diff3/(eb*length2))*dot(grad(u_3),grad(v_3))*dx(1)
	
	
	F += (u_4 - un_4)*v_4*dx
	F += (u_5 - un_5)*v_5*dx
	F += (u_6 - un_6)*v_6*dx
	
	
	#F += (u_4 - un_4)*v_4*dx(1)
	#F += (u_5 - un_5)*v_5*dx(1)
	#F += (u_6 - un_6)*v_6*dx(1)
	
	#### Old F equations
	
	#F += + dk*k1f*u_1*u_6*v_1*dx + dk*k1f*u_1*u_6*v_6*dx - dk*k1f*u_1*u_6*v_4*dx
	
	#F += + dk*k1b*u_4*v_4*dx - dk*k1b*u_4*v_1*dx - dk*k1b*u_4*v_6*dx
	
	#F += + dk*k2f*u_2*(u_6**2)*v_2*dx + dk*2*k2f*u_2*(u_6**2)*v_6*dx - dk*2*k2f*u_1*u_6*v_5*dx
	
	#F += + dk*2*k2b*(u_5**2)*v_5*dx - dk*k2b*(u_5**2)*v_2*dx - dk*2*k2b*(u_5**2)*v_6*dx 
	
	#### dx(1) section
	
	#F += - dk*(1.0* k1b*(u_4**1.0)*v_1*dx)+ dk*(1.0* k1f*(u_1**1.0)*(u_6**1.0)*v_1*dx(1))
	#F += - dk*(1.0* k1b*(u_4**1.0)*v_6*dx(1))+ dk*(1.0* k1f*(u_1**1.0)*(u_6**1.0)*v_6*dx(1))
	#F += + dk*(1.0* k1b*(u_4**1.0)*v_4*dx(1))- dk*(1.0* k1f*(u_1**1.0)*(u_6**1.0)*v_4*dx(1))
	#F += - dk*(1.0* k2b*(u_5**2.0)*v_2*dx(1))+ dk*(1.0* k2f*(u_2**1.0)*(u_6**2.0)*v_2*dx(1))
	#F += - dk*(2.0* k2b*(u_5**2.0)*v_6*dx(1))+ dk*(2.0* k2f*(u_2**1.0)*(u_6**2.0)*v_6*dx(1))
	#F += + dk*(2.0* k2b*(u_5**2.0)*v_5*dx(1))- dk*(2.0* k2f*(u_2**1.0)*(u_6**2.0)*v_5*dx(1))
	#F += - dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_1*dx(1))+ dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_1*dx(1))
	#F += - dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_5*dx(1))+ dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_5*dx(1))
	#F += + dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_3*dx(1))- dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_3*dx(1))
	#F += + dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_6*dx(1))- dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_6*dx(1))
	#F += - dk*(1.0* k4b*(u_3**1.0)*(u_6**2.0)*v_4*dx(1))+ dk*(1.0* k4f*(u_4**1.0)*(u_5**1.0)*v_4*dx(1))
	#F += - dk*(1.0* k4b*(u_3**1.0)*(u_6**2.0)*v_5*dx(1))+ dk*(1.0* k4f*(u_4**1.0)*(u_5**1.0)*v_5*dx(1))
	#F += + dk*(1.0* k4b*(u_3**1.0)*(u_6**2.0)*v_3*dx(1))- dk*(1.0* k4f*(u_4**1.0)*(u_5**1.0)*v_3*dx(1))
	#F += + dk*(2.0* k4b*(u_3**1.0)*(u_6**2.0)*v_6*dx(1))- dk*(2.0* k4f*(u_4**1.0)*(u_5**1.0)*v_6*dx(1))
	#
	#F += - dk*(1.0* k1b*(un_4**1.0)*v_1*dx(1))+ dk*(1.0* k1f*(un_1**1.0)*(un_6**1.0)*v_1*dx(1))
	#F += - dk*(1.0* k1b*(un_4**1.0)*v_6*dx(1))+ dk*(1.0* k1f*(un_1**1.0)*(un_6**1.0)*v_6*dx(1))
	#F += + dk*(1.0* k1b*(un_4**1.0)*v_4*dx(1))- dk*(1.0* k1f*(un_1**1.0)*(un_6**1.0)*v_4*dx(1))
	#F += - dk*(1.0* k2b*(un_5**2.0)*v_2*dx(1))+ dk*(1.0* k2f*(un_2**1.0)*(un_6**2.0)*v_2*dx(1))
	#F += - dk*(2.0* k2b*(un_5**2.0)*v_6*dx(1))+ dk*(2.0* k2f*(un_2**1.0)*(un_6**2.0)*v_6*dx(1))
	#F += + dk*(2.0* k2b*(un_5**2.0)*v_5*dx(1))- dk*(2.0* k2f*(un_2**1.0)*(un_6**2.0)*v_5*dx(1))
	#F += - dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_1*dx(1))+ dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_1*dx(1))
	#F += - dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_5*dx(1))+ dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_5*dx(1))
	#F += + dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_3*dx(1))- dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_3*dx(1))
	#F += + dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_6*dx(1))- dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_6*dx(1))
	#F += - dk*(1.0* k4b*(un_3**1.0)*(un_6**2.0)*v_4*dx(1))+ dk*(1.0* k4f*(un_4**1.0)*(un_5**1.0)*v_4*dx(1))
	#F += - dk*(1.0* k4b*(un_3**1.0)*(un_6**2.0)*v_5*dx(1))+ dk*(1.0* k4f*(un_4**1.0)*(un_5**1.0)*v_5*dx(1))
	#F += + dk*(1.0* k4b*(un_3**1.0)*(un_6**2.0)*v_3*dx(1))- dk*(1.0* k4f*(un_4**1.0)*(un_5**1.0)*v_3*dx(1))
	#F += + dk*(2.0* k4b*(un_3**1.0)*(un_6**2.0)*v_6*dx(1))- dk*(2.0* k4f*(un_4**1.0)*(un_5**1.0)*v_6*dx(1))
	
	
	#### CO Adsorption Process
	
	F += - dk*(1.0* k1b*(u_4**1.0)*v_1*dx)+ dk*(1.0* k1f*(u_1**1.0)*(u_6**1.0)*v_1*dx)
	F += - dk*(1.0* k1b*(u_4**1.0)*v_6*dx)+ dk*(1.0* k1f*(u_1**1.0)*(u_6**1.0)*v_6*dx)
	F += + dk*(1.0* k1b*(u_4**1.0)*v_4*dx)- dk*(1.0* k1f*(u_1**1.0)*(u_6**1.0)*v_4*dx)
	
	F += - dk*(1.0* k1b*(un_4**1.0)*v_1*dx)+ dk*(1.0* k1f*(un_1**1.0)*(un_6**1.0)*v_1*dx)
	F += - dk*(1.0* k1b*(un_4**1.0)*v_6*dx)+ dk*(1.0* k1f*(un_1**1.0)*(un_6**1.0)*v_6*dx)
	F += + dk*(1.0* k1b*(un_4**1.0)*v_4*dx)- dk*(1.0* k1f*(un_1**1.0)*(un_6**1.0)*v_4*dx)
	
	#### O2 Adsorption Process
	
	F += - dk*(1.0* k2b*(u_5**2.0)*v_2*dx)+ dk*(1.0* k2f*(u_2**1.0)*(u_6**2.0)*v_2*dx)
	F += - dk*(2.0* k2b*(u_5**2.0)*v_6*dx)+ dk*(2.0* k2f*(u_2**1.0)*(u_6**2.0)*v_6*dx)
	F += + dk*(2.0* k2b*(u_5**2.0)*v_5*dx)- dk*(2.0* k2f*(u_2**1.0)*(u_6**2.0)*v_5*dx)
	
	F += - dk*(1.0* k2b*(un_5**2.0)*v_2*dx)+ dk*(1.0* k2f*(un_2**1.0)*(un_6**2.0)*v_2*dx)
	F += - dk*(2.0* k2b*(un_5**2.0)*v_6*dx)+ dk*(2.0* k2f*(un_2**1.0)*(un_6**2.0)*v_6*dx)
	F += + dk*(2.0* k2b*(un_5**2.0)*v_5*dx)- dk*(2.0* k2f*(un_2**1.0)*(un_6**2.0)*v_5*dx)
	
	#### Eley-Rideal Reaction
	
	F += - dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_1*dx)+ dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_1*dx)
	F += - dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_5*dx)+ dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_5*dx)
	F += + dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_3*dx)- dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_3*dx)
	F += + dk*(1.0* k3b*(u_3**1.0)*(u_6**1.0)*v_6*dx)- dk*(1.0* k3f*(u_1**1.0)*(u_5**1.0)*v_6*dx)
	
	F += - dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_1*dx)+ dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_1*dx)
	F += - dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_5*dx)+ dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_5*dx)
	F += + dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_3*dx)- dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_3*dx)
	F += + dk*(1.0* k3b*(un_3**1.0)*(un_6**1.0)*v_6*dx)- dk*(1.0* k3f*(un_1**1.0)*(un_5**1.0)*v_6*dx) 
	
	#### Lang-Rideal Reaction

	F += - dk*(1.0* k4b*(u_3**1.0)*(u_6**2.0)*v_4*dx)+ dk*(1.0* k4f*(u_4**1.0)*(u_5**1.0)*v_4*dx)
	F += - dk*(1.0* k4b*(u_3**1.0)*(u_6**2.0)*v_5*dx)+ dk*(1.0* k4f*(u_4**1.0)*(u_5**1.0)*v_5*dx)
	F += + dk*(1.0* k4b*(u_3**1.0)*(u_6**2.0)*v_3*dx)- dk*(1.0* k4f*(u_4**1.0)*(u_5**1.0)*v_3*dx)
	F += + dk*(2.0* k4b*(u_3**1.0)*(u_6**2.0)*v_6*dx)- dk*(2.0* k4f*(u_4**1.0)*(u_5**1.0)*v_6*dx)
	
	F += - dk*(1.0* k4b*(un_3**1.0)*(un_6**2.0)*v_4*dx)+ dk*(1.0* k4f*(un_4**1.0)*(un_5**1.0)*v_4*dx)
	F += - dk*(1.0* k4b*(un_3**1.0)*(un_6**2.0)*v_5*dx)+ dk*(1.0* k4f*(un_4**1.0)*(un_5**1.0)*v_5*dx)
	F += + dk*(1.0* k4b*(un_3**1.0)*(un_6**2.0)*v_3*dx)- dk*(1.0* k4f*(un_4**1.0)*(un_5**1.0)*v_3*dx)
	F += + dk*(2.0* k4b*(un_3**1.0)*(un_6**2.0)*v_6*dx)- dk*(2.0* k4f*(un_4**1.0)*(un_5**1.0)*v_6*dx)

	f = Expression('x[0]>=((1-0.5) - 0.5*0.01) - 2*0.0025 && x[0]<=((1-0.5) + 0.5*0.01) + 2*0.0025 ? 1 : 0', degree=0)
	#f = Expression('x[0]>=((1-0.5) - 0.5*0.01) && x[0]<=((1-0.5) + 0.5*0.01) ? 1 : 0', degree=0)

	cntrls = []

	initialComp1 = Constant(initialCompositions[0])
	initialComp2 = Constant(initialCompositions[1])
	#initialComp3 = Constant(277.77778)
	initialComp3 = Constant(initialCompositions[2])

	#cntrls.append(Control(initialComp1))
	#cntrls.append(Control(initialComp2))
	#cntrls.append(Control(initialComp3))

	F += -initialComp1*f*exp(-(constantT)*(constantT)/(4*0.000000001))*v_4*dx#*b0Test2#
	F += -initialComp2*f*exp(-(constantT)*(constantT)/(4*0.000000001))*v_5*dx#*b0Test2#
	F += -initialComp3*f*exp(-(constantT)*(constantT)/(4*0.000000001))*v_6*dx#*b0Test2#
	
	sST1 = Constant(currentTimes[0]) 
	sSI1 = Constant(pulseIntensities[0]) 

	#cntrls.append(Control(sST1))
	cntrls.append(Control(sSI1))

	b0Test2 = Expression('x[0]<0.02 ? 1 : 0', degree=0)# 

	F += -sSI1*b0Test2*exp(-(constantT - sST1)*(constantT - sST1)/(4*0.00000001))*v_1*dx#*b0Test2

	sST2 = Constant(currentTimes[1]) 
	sSI2 = Constant(pulseIntensities[1]) 
	#cntrls.append(Control(sST2))
	cntrls.append(Control(sSI2))

	F += -sSI2*b0Test2*exp(-(constantT - sST2)*(constantT - sST2)/(4*0.00000001))*v_2*dx#*b0Test2

	sST3 = Constant(0.0)
	sSI3 = Constant(0)
	#cntrls.append(Control(sST3))
	#cntrls.append(Control(sSI3))

	F += -sSI3*b0Test2*exp(-(constantT - sST3)*(constantT - sST3)/(4*0.00000001))*v_3*dx#*b0Test2
	
	J = derivative(F,u)
	
	problem = NonlinearVariationalProblem(F,u,bcs,J)
	solver = NonlinearVariationalSolver(problem)
		
	solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
	solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10
	
	fig,ax = plt.subplots()
	
	newList = np.arange(0, 1, 1/401).tolist()
	
	Y1 = V.sub(0).collapse()
	Y2 = V.sub(1).collapse()
	Y3 = V.sub(2).collapse()
	
	newTest1 = project(u_1,Y1)
	y1 = interpolate(newTest1, Y1)
	
	newTest2 = project(u_2,Y2)
	y2 = interpolate(newTest2, Y2)
	
	newTest3 = project(u_3,Y3)
	y3 = interpolate(newTest3, Y3)
	
	data1 = []
	data2 = []
	data3 = []
	timeTotal = []
	
	startTime = time.time()

	eleyRate = []
	langRate = []
	
	while t <= T:
		timeTotal.append(t)
		#print(type(b0Test))
		newTest1 = project(u_6,Y1)
		newData = interpolate(newTest1, Y1).compute_vertex_values()
		#data1.append(newData[0])
		#ax.plot(newData,color='r')

		eleyTest = project((1.0* k3f*(u_1**1.0)*(u_5**1.0)),Y1)
		newData = interpolate(eleyTest, Y1).compute_vertex_values()
		eleyRate.append(newData[200])

		langTest = project((1.0* k4f*(u_4**1.0)*(u_5**1.0)),Y1)
		newData = interpolate(langTest, Y1).compute_vertex_values()
		langRate.append(newData[200])

		#eleyRate.append()
		#langRate.append()

		newTest1 = project(u_1,Y1)
		newData = interpolate(newTest1, Y1).compute_vertex_values()
		data1.append(newData[399]*toFlux[0])
		#ax.plot(newData,color='r')

		newTest2 = project(u_2,Y2)
		newData = interpolate(newTest2, Y2).compute_vertex_values()
		data2.append(newData[399]*toFlux[1])
		#ax.plot(newList,newData,color='b')

		newTest3 = project(un_3,Y3)
		newData = interpolate(newTest3, Y3).compute_vertex_values()
		data3.append(newData[399]*toFlux[2])
		#ax.plot(newData,color='g')

		solver.solve()

		if t > 0 and t < 0.07:
			##jfunc_2 += assemble(inner((1.0* k3f*(u_1**1.0)*(u_5**1.0)),Constant(1))*dx)
			#jfunc_2 += assemble(inner((1.0* k3f*(u_1**1.0)*(u_5**1.0) / (1.0* k4f*(u_4**1.0)*(u_5**1.0))),Constant(1))*dP)
			jfunc_2 += assemble(inner((1.0* k4f*(u_4**1.0)*(u_5**1.0)),Constant(1))*dx)

		else:
			##jfunc_2 = assemble(inner((1.0* k3f*(u_1**1.0)*(u_5**1.0)),Constant(1))*dx)
			#jfunc_2 = assemble(inner((1.0* k3f*(u_1**1.0)*(u_5**1.0) / (1.0* k4f*(u_4**1.0)*(u_5**1.0))) ,Constant(1))*dP) # 
			jfunc_2 = assemble(inner((1.0* k4f*(u_4**1.0)*(u_5**1.0)) ,Constant(1))*dx)
	
		#print(jfunc_2)
		#plot(u_2)
		#plot(u_3)
		
		t += dt
		constantT.assign(t)
		u_n.assign(u)

	print(jfunc_2)


	ax.plot(timeTotal,data1,color='r',label='CO')
	ax.plot(timeTotal,data2,color='g',label='O2')
	ax.plot(timeTotal,data3,color='b',label='CO2')

	def add_subplot_axes(ax,rect,axisbg='w'):
		
		"""
		Generates the subplot to help visualize some other quantitiy
		"""

		fig = plt.gcf()
		box = ax.get_position()
		width = box.width
		height = box.height
		inax_position  = ax.transAxes.transform(rect[0:2])
		transFigure = fig.transFigure.inverted()
		infig_position = transFigure.transform(inax_position)    
		
		x = infig_position[0]
		y = infig_position[1]
		width *= rect[2]
		height *= rect[3]
		subax = fig.add_axes([x,y,width,height])
		
		x_labelsize = subax.get_xticklabels()[0].get_size()
		y_labelsize = subax.get_yticklabels()[0].get_size()
		x_labelsize *= rect[2]**0.5
		y_labelsize *= rect[3]**0.5
		
		subax.xaxis.set_tick_params(labelsize=x_labelsize)
		subax.yaxis.set_tick_params(labelsize=y_labelsize)
		
		return subax


	#ax = fig.add_subplot(111)
	rect = [0.4,0.68,0.4,0.3]
	ax1 = add_subplot_axes(ax,rect)

	ax1.plot(timeTotal,eleyRate,color='r',label='eley')
	ax1.plot(timeTotal,langRate,color='g',label='lang')

	ax.legend()
	ax1.legend()

	ax.set_ylim(0,0.12)
	ax1.set_ylim(0,5)

	fig.canvas.draw()
	image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
	image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

	figList.append(image)

	#plt.show()
	#sys.exit()


	#print("FINAL TIME")
	#print(time.time() - startTime)

	# Define the callback
	

	#rf_2 = ReducedFunctional(-jfunc_2, cntrls,derivative_cb_post=derivCB)

	return jfunc_2,cntrls

pulseTimes = [0.0,0.0]
pulseIntensities = [1,1]
initialCompositions = [0,0,150]

def derivCB(j,dj,m):
	djv = [v.values()[0] for v in dj]
	mv = [v.values()[0] for v in m]
	print(j)
	print(djv)
	print(mv)

#jfunc_2,cntrls = callSim(pulseTimes,pulseIntensities,pulseIntensities)
#
#low_bounds = []
#up_bounds = []
#
#for gt in range(0,len(cntrls)):
#	low_bounds.append(0)
#	if gt == 0:
#		up_bounds.append(0.1)
#	else:
#		up_bounds.append(5)
#
#rf_2 = ReducedFunctional(-jfunc_2, cntrls,derivative_cb_post=derivCB)
#
#u_opt_2 = minimize(rf_2, method = 'L-BFGS-B',bounds = (low_bounds,up_bounds),tol=1e-12, options={"ftol":1e-12,"gtol":1e-12})
#
#sys.exit()



for k in range(0,8):
	print(pulseTimes)
	print(pulseIntensities)
	jfunc_2,cntrls = callSim(pulseTimes,pulseIntensities,initialCompositions)
	pulseIntensities[1] += 0.5

	#imageio.mimsave('./output.gif', figList, fps=4)

#for jnum,j in enumerate(djv):
#		if j < 0:
#			pulseTimes[jnum] += 0.002
#		else:
#			pulseTimes[jnum] -= 0.002
#
#		if pulseTimes[1] < 0.0:
#			imageio.mimsave('./output.gif', figList, fps=4)
#			sys.exit()
#	pulseTimes[0] = 0.0

#while 1 != 0:
#
#	print('Pulse Time')
#	print(pulseTimes[0])
#	print(pulseTimes[1])
#	jfunc_2,cntrls = callSim(pulseTimes,pulseIntensities,initialCompositions)
#	dJdm = compute_gradient(-jfunc_2, cntrls)
#	djv = [v.values()[0] for v in dJdm]
#	print('derivative value')
#	print(djv)
#	for jnum,j in enumerate(djv):
#		if j < 0:
#			pulseTimes[jnum] += 0.002
#		else:
#			pulseTimes[jnum] -= 0.002
#
#		if pulseTimes[1] < 0.0:
#			imageio.mimsave('./output.gif', figList, fps=4)
#			sys.exit()
#	pulseTimes[0] = 0.0


imageio.mimsave('./output_intensity.gif', figList, fps=4)
sys.exit()
u_opt_2 = minimize(rf_2, method = 'L-BFGS-B',bounds = (low_bounds,up_bounds),tol=1e-12, options={"ftol":1e-12,"gtol":1e-12})
