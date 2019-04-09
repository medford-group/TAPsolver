# -*- coding: utf-8 -*-
"""This module provides a small Python layer on top of the C++
VariationalProblem/Solver classes as well as the solve function."""

# Copyright (C) 2011 Anders Logg
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Marie E. Rognes, 2011.
# Modified by Johan Hake, 2011.

# Import SWIG-generated extension module (DOLFIN C++)
import dolfin.cpp as cpp

# Import UFL
import ufl

# Local imports
from dolfin.fem.form import Form
from dolfin.fem.formmanipulations import derivative
import six

__all__ = ["LinearVariationalProblem",
           "LinearVariationalSolver",
           "LocalSolver",
           "NonlinearVariationalProblem",
           "NonlinearVariationalSolver",
           "solve"]

# Problem classes need special handling since they involve JIT
# compilation


class LinearVariationalProblem(cpp.LinearVariationalProblem):

    # Reuse C++ doc-string
    __doc__ = cpp.LinearVariationalProblem.__doc__

    def __init__(self, a, L, u, bcs=None,
                 form_compiler_parameters=None):
        """
        Create linear variational problem a(u, v) = L(v).

        An optional argument bcs may be passed to specify boundary
        conditions.

        Another optional argument form_compiler_parameters may be
        specified to pass parameters to the form compiler.
        """

        # Extract and check arguments
        u = _extract_u(u)
        bcs = _extract_bcs(bcs)

        # Store input UFL forms and solution Function
        self.a_ufl = a
        self.L_ufl = L
        self.u_ufl = u

        # Store form compiler parameters
        form_compiler_parameters = form_compiler_parameters or {}
        self.form_compiler_parameters = form_compiler_parameters

        # Wrap forms (and check if linear form L is empty)
        if L.empty():
            L = cpp.Form(1, 0)
        else:
            L = Form(L, form_compiler_parameters=form_compiler_parameters)
        a = Form(a, form_compiler_parameters=form_compiler_parameters)

        # Initialize C++ base class
        cpp.LinearVariationalProblem.__init__(self, a, L, u, bcs)


class LocalSolver(cpp.LocalSolver):

    # Reuse C++ doc-string
    __doc__ = cpp.LocalSolver.__doc__

    def __init__(self, a, L=None, solver_type=cpp.LocalSolver.SolverType_LU):
        """Create a local (cell-wise) solver for a linear variational problem
        a(u, v) = L(v).

        """

        # Store input UFL forms and solution Function
        self.a_ufl = a
        self.L_ufl = L

        # Wrap as DOLFIN forms
        a = Form(a)
        if L is None:
            # Initialize C++ base class
            cpp.LocalSolver.__init__(self, a, solver_type)
        else:
            if L.empty():
                L = cpp.Form(1, 0)
            else:
                L = Form(L)

            # Initialize C++ base class
            cpp.LocalSolver.__init__(self, a, L, solver_type)


class NonlinearVariationalProblem(cpp.NonlinearVariationalProblem):

    # Reuse C++ doc-string
    __doc__ = cpp.NonlinearVariationalProblem.__doc__

    def __init__(self, F, u, bcs=None, J=None,
                 form_compiler_parameters=None):
        """
        Create nonlinear variational problem F(u; v) = 0.

        Optional arguments bcs and J may be passed to specify boundary
        conditions and the Jacobian J = dF/du.

        Another optional argument form_compiler_parameters may be
        specified to pass parameters to the form compiler.
        """

        # Extract and check arguments
        u = _extract_u(u)
        bcs = _extract_bcs(bcs)

        # Store input UFL forms and solution Function
        self.F_ufl = F
        self.J_ufl = J
        self.u_ufl = u

        # Store form compiler parameters
        form_compiler_parameters = form_compiler_parameters or {}
        self.form_compiler_parameters = form_compiler_parameters

        # Wrap forms
        F = Form(F, form_compiler_parameters=form_compiler_parameters)
        if J is not None:
            J = Form(J, form_compiler_parameters=form_compiler_parameters)

        # Initialize C++ base class
        cpp.NonlinearVariationalProblem.__init__(self, F, u, bcs, J)


# Solver classes are imported directly
from dolfin.cpp import LinearVariationalSolver, NonlinearVariationalSolver
from dolfin.fem.adaptivesolving import AdaptiveLinearVariationalSolver
from dolfin.fem.adaptivesolving import AdaptiveNonlinearVariationalSolver


# Solve function handles both linear systems and variational problems
def solve(*args, **kwargs):
    """Solve linear system Ax = b or variational problem a == L or F == 0.

    The DOLFIN solve() function can be used to solve either linear
    systems or variational problems. The following list explains the
    various ways in which the solve() function can be used.

    *1. Solving linear systems*

    A linear system Ax = b may be solved by calling solve(A, x, b),
    where A is a matrix and x and b are vectors. Optional arguments
    may be passed to specify the solver method and preconditioner.
    Some examples are given below:

    .. code-block:: python

        solve(A, x, b)
        solve(A, x, b, "lu")
        solve(A, x, b, "gmres", "ilu")
        solve(A, x, b, "cg", "hypre_amg")

    Possible values for the solver method and preconditioner depend
    on which linear algebra backend is used and how that has been
    configured.

    To list all available LU methods, run the following command:

    .. code-block:: python

        list_lu_solver_methods()

    To list all available Krylov methods, run the following command:

    .. code-block:: python

        list_krylov_solver_methods()

    To list all available preconditioners, run the following command:

    .. code-block:: python

        list_krylov_solver_preconditioners()

    To list all available solver methods, including LU methods, Krylov
    methods and, possibly, other methods, run the following command:

    .. code-block:: python

        list_linear_solver_methods()

    *2. Solving linear variational problems*

    A linear variational problem a(u, v) = L(v) for all v may be
    solved by calling solve(a == L, u, ...), where a is a bilinear
    form, L is a linear form, u is a Function (the solution). Optional
    arguments may be supplied to specify boundary conditions or solver
    parameters. Some examples are given below:

    .. code-block:: python

        solve(a == L, u)
        solve(a == L, u, bcs=bc)
        solve(a == L, u, bcs=[bc1, bc2])

        solve(a == L, u, bcs=bcs,
              solver_parameters={"linear_solver": "lu"},
              form_compiler_parameters={"optimize": True})

    For available choices for the 'solver_parameters' kwarg, look at:

    .. code-block:: python

        info(LinearVariationalSolver.default_parameters(), True)

    *3. Solving nonlinear variational problems*

    A nonlinear variational problem F(u; v) = 0 for all v may be
    solved by calling solve(F == 0, u, ...), where the residual F is a
    linear form (linear in the test function v but possibly nonlinear
    in the unknown u) and u is a Function (the solution). Optional
    arguments may be supplied to specify boundary conditions, the
    Jacobian form or solver parameters. If the Jacobian is not
    supplied, it will be computed by automatic differentiation of the
    residual form. Some examples are given below:

    .. code-block:: python

        solve(F == 0, u)
        solve(F == 0, u, bcs=bc)
        solve(F == 0, u, bcs=[bc1, bc2])

        solve(F == 0, u, bcs, J=J,
              solver_parameters={"linear_solver": "lu"},
              form_compiler_parameters={"optimize": True})


    For available choices for the 'solver_parameters' kwarg, look at:

    .. code-block:: python

        info(NonlinearVariationalSolver.default_parameters(), True)

    *4. Solving linear/nonlinear variational problems adaptively

    Linear and nonlinear variational problems maybe solved adaptively,
    with automated goal-oriented error control. The automated error
    control algorithm is based on adaptive mesh refinement in
    combination with automated generation of dual-weighted
    residual-based error estimates and error indicators.

    An adaptive solve may be invoked by giving two additional
    arguments to the solve call, a numerical error tolerance and a
    goal functional (a Form).

    .. code-block:: python

        M = u*dx()
        tol = 1.e-6

        # Linear variational problem
        solve(a == L, u, bcs=bc, tol=tol, M=M)

        # Nonlinear problem:
        solve(F == 0, u, bcs=bc, tol=tol, M=M)

    """

    assert(len(args) > 0)

    # Call adaptive solve if we get a tolerance
    if "tol" in kwargs:
        _solve_varproblem_adaptive(*args, **kwargs)

    # Call variational problem solver if we get an equation (but not a
    # tolerance)
    elif isinstance(args[0], ufl.classes.Equation):
        _solve_varproblem(*args, **kwargs)

    # Default case, just call the wrapped C++ solve function
    else:
        if kwargs:
            cpp.dolfin_error("solving.py",
                             "solve linear algebra problem",
                             "Not expecting keyword arguments when solving "
                             "linear algebra problem")

        return cpp.la_solve(*args)


def _solve_varproblem(*args, **kwargs):
    "Solve variational problem a == L or F == 0"

    # Extract arguments
    eq, u, bcs, J, tol, M, form_compiler_parameters, solver_parameters \
        = _extract_args(*args, **kwargs)

    # Solve linear variational problem
    if isinstance(eq.lhs, ufl.Form) and isinstance(eq.rhs, ufl.Form):

        # Create problem
        problem = LinearVariationalProblem(eq.lhs, eq.rhs, u, bcs,
                                           form_compiler_parameters=form_compiler_parameters)

        # Create solver and call solve
        solver = LinearVariationalSolver(problem)
        solver.parameters.update(solver_parameters)
        solver.solve()

    # Solve nonlinear variational problem
    else:

        # Create Jacobian if missing
        if J is None:
            cpp.info("No Jacobian form specified for nonlinear variational problem.")
            cpp.info("Differentiating residual form F to obtain Jacobian J = F'.")
            F = eq.lhs
            J = derivative(F, u)

        # Create problem
        problem = NonlinearVariationalProblem(eq.lhs, u, bcs, J,
                                              form_compiler_parameters=form_compiler_parameters)

        # Create solver and call solve
        solver = NonlinearVariationalSolver(problem)
        solver.parameters.update(solver_parameters)
        solver.solve()


def _solve_varproblem_adaptive(*args, **kwargs):
    "Solve variational problem a == L or F == 0 adaptively"

    # Extract arguments
    eq, u, bcs, J, tol, M, form_compiler_parameters, \
        solver_parameters = _extract_args(*args, **kwargs)

    # Check that we received the goal functional
    if M is None:
        cpp.dolfin_error("solving.py",
                         "solve variational problem adaptively",
                         "Missing goal functional")

    # Solve linear variational problem
    if isinstance(eq.lhs, ufl.Form) and isinstance(eq.rhs, ufl.Form):

        # Create problem
        problem = LinearVariationalProblem(eq.lhs, eq.rhs, u, bcs,
                                           form_compiler_parameters=form_compiler_parameters)

        # Create solver and call solve
        solver = AdaptiveLinearVariationalSolver(problem, M)
        solver.parameters.update(solver_parameters)
        solver.solve(tol)

    # Solve nonlinear variational problem
    else:

        # Create Jacobian if missing
        if J is None:
            cpp.info("No Jacobian form specified for nonlinear variational problem.")
            cpp.info("Differentiating residual form F to obtain Jacobian J = F'.")
            F = eq.lhs
            J = derivative(F, u)

        # Create problem
        problem = NonlinearVariationalProblem(eq.lhs, u, bcs, J,
                                              form_compiler_parameters=form_compiler_parameters)

        # Create solver and call solve
        solver = AdaptiveNonlinearVariationalSolver(problem, M)
        solver.parameters.update(solver_parameters)
        solver.solve(tol)


def _extract_args(*args, **kwargs):
    "Common extraction of arguments for _solve_varproblem[_adaptive]"

    # Check for use of valid kwargs
    valid_kwargs = ["bcs", "J", "tol", "M",
                    "form_compiler_parameters", "solver_parameters"]
    for kwarg in six.iterkeys(kwargs):
        if kwarg not in valid_kwargs:
            cpp.dolfin_error("solving.py",
                             "solve variational problem",
                             "Illegal keyword argument \"%s\"; valid keywords are %s" %
                             (kwarg,
                              ", ".join("\"%s\"" % kwarg for kwarg in valid_kwargs)))

    # Extract equation
    if not len(args) >= 2:
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Missing arguments, expecting solve(lhs == rhs, "
                         "u, bcs=bcs), where bcs is optional")
    if len(args) > 3:
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Too many arguments, expecting solve(lhs == rhs, "
                         "u, bcs=bcs), where bcs is optional")

    # Extract equation
    eq = _extract_eq(args[0])

    # Extract solution function
    u = _extract_u(args[1])

    # Extract boundary conditions
    if len(args) > 2:
        bcs = _extract_bcs(args[2])
    elif "bcs" in kwargs:
        bcs = _extract_bcs(kwargs["bcs"])
    else:
        bcs = []

    # Extract Jacobian
    J = kwargs.get("J", None)
    if J is not None and not isinstance(J, ufl.Form):
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Expecting Jacobian J to be a UFL Form")

    # Extract tolerance
    tol = kwargs.get("tol", None)
    if tol is not None and not (isinstance(tol, (float, int)) and tol >= 0.0):
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Expecting tolerance tol to be a non-negative number")

    # Extract functional
    M = kwargs.get("M", None)
    if M is not None and not isinstance(M, ufl.Form):
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Expecting goal functional M to be a UFL Form")

    # Extract parameters
    form_compiler_parameters = kwargs.get("form_compiler_parameters", {})
    solver_parameters = kwargs.get("solver_parameters", {})

    return eq, u, bcs, J, tol, M, form_compiler_parameters, solver_parameters


def _extract_eq(eq):
    "Extract and check argument eq"
    if not isinstance(eq, ufl.classes.Equation):
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Expecting first argument to be an Equation")
    return eq


def _extract_u(u):
    "Extract and check argument u"
    if not isinstance(u, cpp.Function):
        cpp.dolfin_error("solving.py",
                         "solve variational problem",
                         "Expecting second argument to be a Function")
    return u


def _extract_bcs(bcs):
    "Extract and check argument bcs"
    if bcs is None:
        bcs = []
    elif not isinstance(bcs, (list, tuple)):
        bcs = [bcs]
    for bc in bcs:
        if not isinstance(bc, cpp.DirichletBC):
            cpp.dolfin_error("solving.py",
                             "solve variational problem",
                             "Unable to extract boundary condition arguments")
    return bcs
