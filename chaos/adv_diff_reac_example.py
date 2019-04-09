from fenics import *
from dolfin_adjoint import *
import time

T = 5.0            # final time
num_steps = 30    # number of time steps
dt = T / num_steps # time step size
eps = 0.01         # diffusion coefficient
K = 10.0           # reaction rate

parameters["std_out_all_processes"] = False																							
cache_params = {"cache_dir":"*","lib_dir":"*","log_dir":"*","inc_dir":"*","src_dir":"*"}											
set_log_active(False)
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Read mesh from file
mesh = Mesh('./cylinder.xml.gz')

# Define function space for velocity
W = VectorFunctionSpace(mesh, 'P', 2)

# Define function space for system of concentrations
P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

# Define test functions
v_1, v_2 = TestFunctions(V)

# Define functions for velocity and concentrations
w = Function(W)
u = Function(V)
u_n = Function(V)

# Split system functions to access components
u_1, u_2 = split(u)
u_n1, u_n2 = split(u_n)

# Define source terms
f_1 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.1,2)<0.05*0.05 ? 0.1 : 0',
                 degree=1)
f_2 = Expression('pow(x[0]-0.1,2)+pow(x[1]-0.3,2)<0.05*0.05 ? 0.1 : 0',
                 degree=1)
f_3 = Constant(0)

# Define expressions used in variational forms
k = Constant(dt)
K = Constant(K)
eps = Constant(eps)

# Define variational problem
F = ((u_1 - u_n1) / k)*v_1*dx + eps*dot(grad(u_1), grad(v_1))*dx + ((u_2 - u_n2) / k)*v_2*dx + eps*dot(grad(u_2), grad(v_2))*dx# + K*u_1*u_2*v_2*dx - K*u_1*u_2*v_3*dx + K*u_3*v_3*dx - f_1*v_1*dx - f_2*v_2*dx - f_3*v_3*dx

# Create time series for reading velocity data
#timeseries_w = TimeSeries('navier_stokes_cylinder/velocity_series')

# Create VTK files for visualization output
#vtkfile_u_1 = File('reaction_system/u_1.pvd')
#vtkfile_u_2 = File('reaction_system/u_2.pvd')
#vtkfile_u_3 = File('reaction_system/u_3.pvd')

# Create progress bar
progress = Progress('Time-stepping')
set_log_level(PROGRESS)

# Time-stepping
t = 0

start_time = time.time()

for n in range(num_steps):

    # Update current time
    t += dt

    # Read velocity from file
    #timeseries_w.retrieve(w.vector(), t)

    # Solve variational problem for time step
    solve(F == 0, u)
    
    # Save solution to file (VTK)
    #_u_1, _u_2, _u_3 = u.split()
    #vtkfile_u_1 << (_u_1, t)
    #vtkfile_u_2 << (_u_2, t)
    #vtkfile_u_3 << (_u_3, t)

    print("TIME STEP")
    print(t)

    # Update previous solution
    u_n.assign(u)

    # Update progress bar
    progress.update(t / T)

print(time.time() - start_time)


control1 = Control(eps)
control2 = Control(K)

#sens_func = assemble(inner(u,u)*dP(1))
sens_func = assemble(inner(u,u)*dx)
#ens_func = assemble(test_new[k][0]*dP(1))
start_time = time.time()
print("")
print('sensitivity time')
dJdm1,dJdm2 = compute_gradient(sens_func,[control1,control2])#compute_gradient(sens_func,[control1,control2])
print(time.time() - start_time)
print("dJdm for each control")
print(dJdm1.values(),dJdm2.values())
sys.exit()

