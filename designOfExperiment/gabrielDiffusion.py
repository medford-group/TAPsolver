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

numberOfCells = 100
mesh = UnitIntervalMesh(numberOfCells)

print("Plotting a UnitIntervalMesh")
plot(mesh, title="Unit interval")
plt.show()

V = FunctionSpace(mesh,'P',1)
## 'P' is for the lagrange element
## 1 is for the degree of the element

v = TestFunction(V)
u = Function(V)
u_n = Function(V)

## First, define the boundaries at the entrance of the reactor
tol = 1e-14

# Boundary at the entrance
def boundary_L(x, on_boundary):
	return on_boundary and near(x[0],0,tol)

# Boundary at exit
def boundary_R(x, on_boundary):
	return on_boundary and near(x[0],1,tol)


## Set the outlet concentration to 0
BCL = DirichletBC(V,Constant(1),boundary_L)
BCR = DirichletBC(V,Constant(1),boundary_R)

bcs = [BCR]


dz = 1/numberOfCells
dk = Constant(1)
diff = Constant(0.001)

# First define the diffusion for all terms
F = (u - u_n)*v*dx + diff*dot(grad(u),grad(v))*dx


J = derivative(F,u)

problem = NonlinearVariationalProblem(F,u,bcs,J)
solver = NonlinearVariationalSolver(problem)

solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8
solver.parameters["newton_solver"]["absolute_tolerance"] = 1e-10


Y1 = V
newTest1 = project(u,Y1)
y1 = interpolate(newTest1, Y1)

t = 0
T = 2
dt = 0.1

while t <= T:
	#print(u)
	solver.solve()
	newTest1 = project(u,Y1)
	newData = interpolate(newTest1, Y1).compute_vertex_values()
	print(newData)
	print(type(newData))
	plot(u)
	t += dt
	u_n.assign(u)
	
plt.show()