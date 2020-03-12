# FEniCS & Dolfin-adjoint Tutorial

#This is a condensed and adjusted tutorial for FEniCS and Dolfin-adjoint use.
#The primary difference is that all discussions are based around chemical 
#engineering applications.
#
#Further examples can be found at:
#
#[FEniCS Basics](https://fenicsproject.org/tutorial/)
#
#[Dolfin-adjoint Basics](http://www.dolfin-adjoint.org/en/latest/documentation/examples.html)
#
#*Written by Adam C. Yonge*
#
#*January 31, 2020*

## Importing packages

# Both FEniCS and fenics_adjoint (dolfin_adjoint) must be imported. fenics_adjoint overwrites some functions in fenics, 
# so swapping the order would result in errors. Though you can delete any of the imports below 'fenics' 'fenics_adjoint', 
# I've generally found that importing all of these packages is helpful, so I usually just copy and paste them to the top 
# of each fenics script I generate.

from fenics import *
from fenics_adjoint import *
import mpmath
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
#import dijitso
import time
import csv
import ast
import shutil
import sys
import os
import scipy
import pip
import pkg_resources
import ufl
from ufl import sqrt,exp,ln
import warnings

## Defining the mesh

#[Several options for mesh generation are available](https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/built-in_meshes/python/documentation.html).
#
#Feel free to play with and visualize the available meshes below. A unit interval mesh will be used throughout this tutorial, but references in how to manipulate 2D and 3D meshes will be outlined.
#
#The visualization will not work with 3D structures, though.

############# Unit Interval mesh #################

numberOfCells = 30
mesh = UnitIntervalMesh(numberOfCells)

############# Unit Square mesh ###################

## cells in horizontal direction
#nx = 10
## cells in vertical direction
#ny = 15
## Other possible diagonal directions =  “left”, “right”, “left/right”, “crossed”
#diagDirection = 'crossed'
## Generate the mesh
#mesh = UnitSquareMesh(nx, ny, diagonal=diagDirection)

############### Rectangle mesh ###################

#x1 = Point(0,0)
#y1 = Point(3,5)
#nx = 7
##ny = 10
##diagDirection = "crossed"
## Other possible diagonal directions =  “left”, “right”, “left/right”, “crossed”
#mesh = RectangleMesh(x1, y1, nx, ny, diagonal=diagDirection)

############### Unit Cube Mesh ###################

#nx = 3
#ny = 3
#nz = 3
#mesh = UnitCubeMesh(nx, ny, nz)

################## Box Mesh ######################

#p1 = Point(0,0,0)
#p2 = Point(2,2,2)
#nx = 3
#ny = 3
#nz = 3
#mesh = BoxMesh(p1, p2, nx, ny, nz)

################## View Mesh ######################

#print("Plotting a UnitIntervalMesh")
#plot(mesh, title="Unit interval")
#plt.show()

## Defining the elements for trial and test functions

# Next the element space being used must be defined. There are many options that 
# can be found in the [Periodic Table of Finite Elements](http://femtable.org/). 
# The basic one to use is the 'lagrange' family. The optional elements do not always work immediatly.

######### Basic, single monitored u (species or solution). #########

#V = FunctionSpace(mesh,'P',1)
## 'P' is for the lagrange element
## 1 is for the degree of the element

####### For multiple species ###########

CG1 = FiniteElement('CG',mesh.ufl_cell(),1) # Continuous galerkin
P1 = FiniteElement('P',mesh.ufl_cell(),1)
combElements = [P1,P1,P1] # an element for each one of the species
elements = MixedElement(combElements)
V = FunctionSpace(mesh,elements)
V_du = FunctionSpace(mesh,P1)

# Now the test and trial functions can be defined in terms of the function space, V. 
# There are currently only three species that will be manually defined. It is possible 
# to increase the flexibility of the code by using dictionaries and input files.

## Trial Function
## Since we will be working with a time dependent problem, we must define a new and previous function space.
u = Function(V)
u_n = Function(V)

tempA = TestFunctions(V) # generate the test and trial functions
tempB = split(u)
tempC = split(u_n)

v_1 = tempA[0]
v_2 = tempA[1]
v_3 = tempA[2]

u_1 = tempB[0]
u_2 = tempB[1]
u_3 = tempB[2]

un_1 = tempC[0]
un_2 = tempC[1]
un_3 = tempC[2]

## Defining the boundary conditions

# In this example, there are three species of interest. Let's call them A, B and C. We will work with a simple 
# diffusive process with non-catalytic reactions incorporated. Neumann (flux) and direclet (value) boundary conditions. 
# These can be fixed values or defined as expressions that vary with time. It is also possible to define spatially and 
# temporally dependent initial distributions (defined through expressions).
#
#If no boundary is defined, FEniCS will assume a neumann BC of zero.

## First, define the boundaries at the entrance of the reactor
tol = 1e-14

# Boundary at the entrance
def boundary_L(x, on_boundary):
	return on_boundary and near(x[0],0,tol)

# Boundary at exit
def boundary_R(x, on_boundary):
	return on_boundary and near(x[0],1,tol)


#In the above examples, x[0] is included to define the dimension of interest. 0 and 1 are included to define at what point in the reactor (remember, it is a unit interval mesh). 0 for the entrance and 1 for the exit.
#
#If it were a 2D mesh, defining the boundary in the second dimension would look something like this:
#
#def boundary_top(x, on_boundary):
#    
#    return on_boundary and near(x[1],1,tol)

## Set the outlet concentration to 0
BC1f = DirichletBC(V.sub(0),Constant(1),boundary_L)
BC1 = DirichletBC(V.sub(0),Constant(0),boundary_R)
BC2 = DirichletBC(V.sub(1),Constant(0),boundary_R)
BC3 = DirichletBC(V.sub(2),Constant(0),boundary_R)

bcs = [BC1f,BC1,BC2,BC3]

# When working with neumann boundary conditions, the component must be added directly the equation formulation 
# (which will be discussed shortly). A brief, unused example of this is as follows:

## Constructing the diffusion-reaction equation

# Deriving the variational form can be tricky from scratch, but observing the patterns in the expressions 
# and directly implementing them is also possible and often more efficient.

dk = Constant(1)
diff = Constant(0.01)

k1 = Constant(0.1)
k2 = Constant(0.2)

cntrl1 = Control(k1)
cntrl2 = Control(k2)

# First define the diffusion for all terms
F = (u_1 - un_1)*v_1*dx + diff*dot(grad(u_1),grad(v_1))*dx

F += (u_2 - un_2)*v_2*dx + diff*dot(grad(u_2),grad(v_2))*dx

F += (u_3 - un_3)*v_3*dx + diff*dot(grad(u_3),grad(v_3))*dx

F += k1*u_1*v_1*dx - k1*u_1*v_2*dx

F += k2*u_2*v_2*dx - k2*u_2*v_3*dx


## Define the problem to be solved by fenics


J = derivative(F,u)

problem = NonlinearVariationalProblem(F,u,bcs,J)
solver = NonlinearVariationalSolver(problem)

solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10


## Define the objective function zone

dz = 1/numberOfCells

class integration_section(SubDomain):
	def inside(self, x, on_boundary):
		return between(x[0], (1-dz,1.0))

osub = integration_section()
domains = MeshFunction("size_t", mesh,0)
domains.set_all(0)
osub.mark(domains, 1)
dP = Measure('vertex',domain = mesh, subdomain_data=domains)

## Run the simulation & assemble objective function

t = 0
T = 2
dt = 0.1

Y1 = V.sub(0).collapse()
Y2 = V.sub(1).collapse()
Y3 = V.sub(2).collapse()

newTest1 = project(u_1,Y1)
y1 = interpolate(newTest1, Y1)

newTest2 = project(u_2,Y2)
y2 = interpolate(newTest2, Y2)

newTest3 = project(u_3,Y3)
y3 = interpolate(newTest3, Y3)

fig,ax = plt.subplots()

while t <= T:
	print(t)
	solver.solve()
	
	newTest1 = project(u_1,Y1)
	newData = interpolate(newTest1, Y1).compute_vertex_values()
	ax.plot(newData,color='r')
	
	newTest1 = project(u_2,Y1)
	newData = interpolate(newTest1, Y1).compute_vertex_values()
	ax.plot(newData,color='g')
	
	newTest2 = project(u_3,Y2)
	newData = interpolate(newTest2, Y2).compute_vertex_values()    
	ax.plot(newData,color='b')
	
	if t > 0:
		jfunc_2 += assemble(inner(u_1,u_1)*dP(1))
	else:
		jfunc_2 = assemble(inner(u_1,u_1)*dP(1))
	
	t += dt
	u_n.assign(u)
	
plt.show()


## Defining the reduced functional & minimization process

# Define the callback
def derivCB(j,dj,m):
	djv = [v.values()[0] for v in dj]
	mv = [v.values()[0] for v in m]
	print(j)
	print(djv)
	print(mv)

# Define the gradient calculation function
rf_2 = ReducedFunctional(jfunc_2, [cntrl1,cntrl2],derivative_cb_post=derivCB)

## Run minimization

u_opt_2 = minimize(rf_2, method = 'L-BFGS-B')
print('Done')