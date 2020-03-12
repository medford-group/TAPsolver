from fenics import *
from fenics_adjoint import *
import time

set_log_level(30)
tol = 1E-14

#### Defining the initial parameters for each of the three species

controls = []

K_list =[1,1,1]
for j in range(0,len(K_list)):
	K_list[j] = Constant(K_list[j])

bc_list = [1,1,1]
M_list = [3,3,3]
for j in range(0,len(M_list)):
	M_list[j] = Constant(M_list[j])

ka_list = [0.05,0.03,0.01]
for j in range(0,len(ka_list)):
	ka_list[j] = Constant(ka_list[j])
	controls.append(Control(ka_list[j]))

V_list = [0.49,0.5,0.51]
for j in range(0,len(V_list)):
	V_list[j] = Constant(V_list[j])
	controls.append(Control(V_list[j]))

a_list = [3,2,1] 
for j in range(0,len(a_list)):
	a_list[j] = Constant(a_list[j])

b_list =  [3,2,1]
for j in range(0,len(b_list)):
	b_list[j] = Constant(b_list[j])

monitored_gases = 3
length = 12.8
mesh_size = 400

#### Define the elements

graph_data = {}
v_d = {}
u_d = {}
u_nd = {}

mesh = UnitIntervalMesh(mesh_size)
P1 = FiniteElement('P',mesh.ufl_cell(),1)

element = '['

for k in range(0,monitored_gases*2):
	if (k+1) < monitored_gases*2:
		element = element + 'P1,'
	else:
		element = element + 'P1]'

test_new = eval(element)

element = MixedElement(test_new)
V = FunctionSpace(mesh,element)
v_du = FunctionSpace(mesh,P1)
	
#### Define boundary conditions

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0],0,tol)
		
def boundary_R(x, on_boundary):
	return on_boundary and near(x[0],1,tol)

bcs=[]
for k in range(0,monitored_gases):
	bcs.append(DirichletBC(V.sub(k),Constant(float(bc_list[k])),boundary_L))

all_molecules = monitored_gases*2
u = Function(V)
u_n = Function(V)
tempA = TestFunctions(V)
tempB = split(u)
tempC = split(u_n)

for kit in range(0,monitored_gases*2):
	v_d['v_'+str(kit)] = tempA[kit]
	u_d['u_'+str(kit)] = tempB[kit]
	u_nd['u_n'+str(kit)] = tempC[kit]

dz = 1/mesh_size

class integration_section(SubDomain):
	def inside(self, x, on_boundary):
		return between(x[0], (1-dz,1.0))

osub = integration_section()
domains = MeshFunction("size_t", mesh,0)
domains.set_all(0)
osub.mark(domains, 1)
dP = Measure('vertex',domain = mesh, subdomain_data=domains)

#### Define the PDEs in variational form

w_list = []

for k in range(0,monitored_gases):
	W = VectorFunctionSpace(mesh, 'P', 1)
	w = Function(W)
	w.vector()[:] = float(M_list[k])
	w_list.append(w)

sim_time = 25
time_steps = 400

dt = sim_time/time_steps
dk = Constant(dt)

F_str = ''
for k in range(0,monitored_gases):
	F_str += ' (u_d["u_'+str(k)+'"] - u_nd["u_n'+str(k)+'"])*v_d["v_'+str(k)+'"]*dx +'
	F_str += '(Constant(length))*dot(w_list['+str(k)+'], grad(u_d["u_'+str(k)+'"]))*v_d["v_'+str(k)+'"]*dx' 
	F_str += ' + ((u_d["u_'+str(k+monitored_gases)+'"] - u_nd["u_n'+str(k+monitored_gases)+'"]))*v_d["v_'+str(k+monitored_gases)+'"]*dx' 
	F_str += '+ dk*ka_list['+str(k)+']*((V_list['+str(k)+']*b_list['+str(k)+']*u_d["u_'+str(k)+'"])/(1 + a_list['+str(k)+']*u_d["u_'+str(k)+'"]) - u_d["u_'+str(k+monitored_gases)+'"])*v_d["v_'+str(k)+'"]*dx'
	F_str += '- dk*ka_list['+str(k)+']*((V_list['+str(k)+']*b_list['+str(k)+']*u_d["u_'+str(k)+'"])/(1 + a_list['+str(k)+']*u_d["u_'+str(k)+'"]) - u_d["u_'+str(k+monitored_gases)+'"])*v_d["v_'+str(k+monitored_gases)+'"]*dx'
	if k != monitored_gases-1:
		F_str += '+'

#### Defining the controls

tx = Constant(0)
tstart = Constant(0.5)
intensity = Constant(1)

x =  SpatialCoordinate(mesh)
coeff = tstart/intensity/(x[0]**2 + (tx-tstart)**2)

m = Control(tstart)
m_intensity = Control(intensity)

F = eval(F_str)- coeff*v_d["v_0"]*dx
J = derivative(F,u)
problem = NonlinearVariationalProblem(F ,u,bcs,J)
solver = NonlinearVariationalSolver(problem)
solver.parameters["newton_solver"]["relative_tolerance"] = 1.0e-8

#### Solving the forward problem

time_list = []
t = 0
print('Solving Forward Problem')
startTime = time.time()

for n in range(int(time_steps)):
	solver.solve()
	for k_fitting in range(0,monitored_gases):
		if t > 0:
			jfunc_2 += assemble(inner(u_n[k_fitting],u_n[k_fitting])*dP(1))
		else:
			jfunc_2 = assemble(inner(u_n[k_fitting],u_n[k_fitting])*dP(1))

	time_list.append(t)
	t += dt
	tx.assign(t)
	u_n.assign(u)

print('Forward Time: '+str(round(time.time() - startTime,4)))

#### Calculate the gradient through the adjoint approach

currentTime = time.time()
dJdm = compute_gradient(jfunc_2,[m,m_intensity])
print('Gradient Time: '+str(round(time.time() - currentTime,4)))
djv = [v.values()[0] for v in dJdm]
print(djv)