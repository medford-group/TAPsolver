from fenics import *
from fenics_adjoint import *
from mpmath import nsum, exp, inf
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import time
import math as mp
import dijitso
import csv
import sys
import os
import ufl
import warnings
from pyadjoint.enlisting import Enlist



#Define the 
parameters["std_out_all_processes"] = False                                                                                         
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}                                            
#set_log_active(False)
warnings.filterwarnings("ignore", category=DeprecationWarning)


tape = Tape()
tape.clear_tape()
set_working_tape(tape)

# Define the Mesh
mesh_size = 10000
mesh = UnitIntervalMesh(mesh_size) 
### Similarly, you could define a square mesh with 'mesh = UnitSquareMesh(8,8)'
### or mesh = UnitCubeMesh(10, 10, 10)

# Define the function space
#V = FunctionSpace(mesh,'P',1)
P1 = FiniteElement('P',mesh.ufl_cell(),1)
V = FunctionSpace(mesh,MixedElement([P1,P1]))

#K = Constant(10)


#Define the boundary conditions 

tol = 1E-14
def boundary_L(x, on_boundary):
    return on_boundary and near(x[0],0,tol)
### The x[0] is basically just the first dimension (if necessary, you could define x[1] if you had a y dimension, too)
### The 0 just means at the start of the mesh (whereas 1 below means at the end of the mesh)
### tol just indicates how close to the BC (little more nuanced, but it has little impact)

def boundary_R(x, on_boundary):
    return on_boundary and near(x[0],1,tol)

### Since batch reactors have a zero flux boundary condition, nothing is specified

bcs=[]

sim_time = 25
time_steps = 50
#establish the time step for each simulation
dt = sim_time/time_steps
dk = Constant(dt)


# These are the trial functions at the current and previous time steps.
# If you are simulating more than one species, you defined a unique function and 
# associated trial function for each unknown (concentration)

# u at next step (gas species)
u = Function(V)
u_1,u_2 = split(u)

print(type(u_1))
print(type(u_2))
print(type(u[0]))


# u from previous step
u_n = Function(V)
#u_0 = Expression('1',degree=1)
u_n1, u_n2 = split(u_n)
for z in range(0,mesh_size+1):
    u_n.vector()[z*2] = float(1)
    #print(z)
    u_n.vector()[z*2+1] = float(0)
    #print(z+(mesh_size+1))
v = TestFunction(V)
v1, v2 = split(v)

    # Variational form of the PDE


A = Constant(1)
dG = Constant(0.5)
heatR = Constant(0.1)

F = ((u_1 - u_n1)/dk)*v1*dx  + ((u_2 - u_n2)/dk)*v2*dx + A*mp.exp(dG-heatR)*u_1*v1*dx - A*mp.exp(dG - heatR)*u_1*v2*dx - A*mp.exp(dG+heatR)*u_2*v1*dx + A*mp.exp(dG+heatR)*u_2*v2*dx
### The '((u - u_n))*v*dx' term is defining how u (concentration 1 is changing with time)
### The 'dk*K*u_1*u_2*v*dx' term is defining the rate of consumption of C1
### The '((u_2 - u_n2))*v1*dx' term is defining how u2 (concentration 2 is changing with time)
### The 'K*u*u_2*v1*dx' term is defining the rate of formation of C2

#Define the problem to be solved 
J = derivative(F,u)
problem = NonlinearVariationalProblem(F,u,bcs,J)
solver = NonlinearVariationalSolver(problem)


#Define a mesh for the graph (not crucial)
x1=0
x2=mesh_size+1
space_step = 1/mesh_size
list2 = [x/mesh_size for x in range(x1, x2)]
list2 = list(reversed(list2))

#Establish the figure
fig2, ax2 = plt.subplots()

#Step through and solve for each time step
t = 0


#J2 = assemble(inner(u, u) * dx)
#direction = interpolate(Constant(1),V)
#m = [Control(K),Control(K2)]

timeAll = []
dj1 = []
dj2 = []
d2j1 = []
d2j2 = []
J2 = assemble(inner(u, u) * dx)
n = FacetNormal(mesh)

c = dG
Js = []
Jk = []
#print(type(J2))
#sys.exit()
V_du = FunctionSpace(mesh,P1)
w_new = Expression('A',A=Constant(1),degree=0)#Expression("1", degree=0)
w_new2 = interpolate(w_new,V_du)
w3 = project(w_new2,V_du)

w_new = Expression("1", degree=0)
w_new2 = interpolate(w_new,V_du)

amountTime = []
for n in range(time_steps):
    xTime = time.time()
    timeAll.append(n*dt)
    _u_1 = []
    _u_2 = []
    #for z in range(0,mesh_size+1):
    #    _u_1.append(u_n.vector()[z*2])
    #    _u_2.append(u_n.vector()[z*2+1])
    #ax2.plot(list2,_u_1,color='b')
    #ax2.plot(list2,_u_2,color='r')
                                                        
    #Adjoint Stuff
    
    #if n < 5:
    #Js.append(assemble(inner(u, u) * dx))
    Js.append(assemble(2*(inner((A*mp.exp(dG)*u[0] - A*mp.exp(dG)*u[1]), w3) )* dx))
    #Js.append(assemble((inner((A*mp.exp(dG-heatR)*u[0] - A*mp.exp(dG+heatR)*u[1]), w3) )* dx))
    print(assemble((inner((A*mp.exp(dG)*u[0] - A*mp.exp(dG)*u[1]), w3) )* dx))
    #Jk.append(assemble(( K2*inner(u[1],u[1]) - K*inner(u[0],u[0])*dx ) ) )# - K*u[0]), (K2*u[1] - K*u[0])) )
    
    #print(assemble((inner((K2*u[1] - K*u[0]), w3) )* dx))
    #print(assemble((inner((-K2*u[1] + K*u[0]), w3) )* dx))
    #print(assemble((inner((K2*u[1] - K*u[0]), (K2*u[1] - K*u[0])) )* dx))
    #print(assemble(inner(u[0],w3)*dx))
    
    #print(0.5*assemble(inner(u[1],u[1])*dx))
    print(u.vector()[0])
    print(u.vector()[1])
    print(u.vector()[2])
    print(u.vector()[3])
    #time.sleep(1)
    #print(type(assemble(inner(u, u) * dx)))
    #J2 = assemble(inner(u, u) * dx)
    #J3 = assemble(inner(u, u) * dx)
    #dJdm = compute_gradient(J2,m)
    #dj1.append(float(dJdm[0]))
    #dJdm2 = compute_gradient(J3,m)
    #dj2.append(float(dJdm[0]))
    
    #print(float(dJdm[0]))
    ##d2j1.append(float(dJdm[1]))
    ##print(float(dJdm[1]))
                                                                                
    #Hessian Stuff
    #test = Enlist(m)
    ##dJdm2 = compute_hessian(J2, m, dJdm)
    ##dj1.append(float(dJdm2[0]))
    ##dj2.append(float(dJdm2[1]))  

    #TLM Stuff
    #c.tlm_value = K2
    #tape.evaluate_tlm()
    #c.tlm_value = K2
    #tape.evaluate_tlm()
    #J2.block_variable.tlm_value
    #dj2.append(J3.block_variable.tlm_value)
    
    #print(J2.block_variable.tlm_value)
    #print(J3.block_variable.tlm_value)
    #print('next')
    t += dt
    solver.solve()
    u_n.assign(u)
    
    amountTime.append(time.time() - xTime)
    #print(amountTime)
    #print(dj1)
    #print(dj2)
    #print(time.time() - xTime)

c.tlm_value = dG
tape.evaluate_tlm()

for k in Js:
    print(k.block_variable)
sys.exit()
for k in Jk:
    print(k.block_variable.tlm_value)

sys.exit()
c.tlm_value = K
tape.evaluate_tlm()
for k in Jk:
    print(k.block_variable.tlm_value)


#ax2.set_xlabel('$Dimensionless\ Reactor\ Length$')
#ax2.set_ylabel('$Concentration$')
#plt.show()
#sys.exit()
#Establish the figure
fig2, ax2 = plt.subplots()
ax2.plot(timeAll,dj1,color='b')
ax2.plot(timeAll,dj2,color='r')
plt.show()
sys.exit()
#Establish the figure
fig2, ax2 = plt.subplots()
ax2.plot(timeAll,d2j1,color='b')
ax2.plot(timeAll,d2j2,color='r')
plt.show()