import matplotlib.pyplot as plt
from dolfin import *

class CustomProblem(NonlinearProblem):
    def __init__(self, a, L, bcs):
        NonlinearProblem.__init__(self)
        self.a = a
        self.L = L
        self.bcs = bcs

    def F(self, b, x):
        assembler = SystemAssembler(self.a, self.L, self.bcs)
        assembler.assemble(b, x)

    def J(self, A, x):
        assembler = SystemAssembler(self.a, self.L, self.bcs)
        assembler.assemble(A)


class CustomSolver(NewtonSolver):
    def converged(self, residual, problem, iteration):
        # Some custom convergence criterion here
        rnorm = residual.norm("l2")
        if rnorm < 1e-6:
            print("%d Residual norm ;) %.3e" % (iteration, rnorm))
            return True
        print("%d Residual norm :( %.3e" % (iteration, rnorm))
        return False

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return abs(x[0] - 1.0) < DOLFIN_EPS and on_boundary

mesh = UnitSquareMesh(32, 32)

V = FunctionSpace(mesh, "CG", 1)

g = Constant(1.0)
bc = DirichletBC(V, g, DirichletBoundary())
bcs = []
bcs.append(bc)
bcs.append(bc)
u = Function(V)
v = TestFunction(V)
f = Expression("x[0]*sin(x[1])", degree=2)
F = inner((1 + u**2)*grad(u), grad(v))*dx - f*v*dx

solver = CustomSolver()
problem = CustomProblem(derivative(F, u), F, bcs)
print(type(u.vector()))
solver.solve(problem, u.vector())

plot(u, title="Solution")
#plt.show()​