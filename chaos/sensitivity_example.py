from fenics import *
from fenics_adjoint import *
import time

n = 80
mesh = UnitSquareMesh(n, n)
V = VectorFunctionSpace(mesh, "CG", 2)

u = project(Expression(("sin(2*pi*x[0])", "cos(2*pi*x[1])"), degree=2),  V)
control = Control(u)

u_next = Function(V)
v = TestFunction(V)

nu = Constant(0.0001)

timestep = Constant(0.01)

F = (inner((u_next - u)/timestep, v)
     + inner(grad(u_next)*u_next, v)
     + nu*inner(grad(u_next), grad(v)))*dx

bc = DirichletBC(V, (0.0, 0.0), "on_boundary")
start_time = time.time()
t = 0.0
end = 0.1
while (t <= end):
    solve(F == 0, u_next, bc)
    u.assign(u_next)
    t += float(timestep)
print(time.time() - start_time)

start_time = time.time()
J = assemble(inner(u, u)*dx)
dJdu, dJdnu = compute_gradient(J, [control, Control(nu)])
print(time.time() - start_time)

print(dJdu,dJdnu)