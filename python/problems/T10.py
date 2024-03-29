"""
T10 is a bi-objective mixed-integer convex test instance with a
non-quadratic objective function and quadratic constraint functions.
"""

import numpy as np
from pyomo.environ import *

class structure():
    pass

def build_model(_):

    model = ConcreteModel()
    model_parameters = structure()

    # Dimension of decision and criterion space
    n = 4         # Continuous variables
    m = 4         # Integer variables
    p = 2         # Dimension criterion space
    q = 4         # Number of constraints

    model_parameters.n = n
    model_parameters.m = m
    model_parameters.p = p
    model_parameters.q = q

    # Problem type
    model_parameters.is_convex = True
    model_parameters.is_quadratic = True

    # Variables
    model.x = Var(range(n+m))
    for i in range(n):
        model.x[i].domain = Reals
        model.x[i].bounds = (-20,20)
    for j in range(n,n+m):
        model.x[j].domain = Integers
        model.x[j].bounds = (-20,20)

    # Objective functions
    model.f = ObjectiveList()
    model.f.add(expr = model.x[0] + model.x[2] + model.x[4] + exp(model.x[6]) - 1).deactivate()
    model.f.add(expr = model.x[1] + model.x[3] + model.x[5] + model.x[7]).deactivate()

    # Derivatives of objective functions
    model.Df = Expression(range(p), range(n+m))
    model.Df[0, 0] = 1
    model.Df[0, 1] = 0
    model.Df[0, 2] = 1
    model.Df[0, 3] = 0
    model.Df[0, 4] = 1
    model.Df[0, 5] = 0
    model.Df[0, 6] = exp(model.x[6])
    model.Df[0, 7] = 0

    model.Df[1, 0] = 0
    model.Df[1, 1] = 1
    model.Df[1, 2] = 0
    model.Df[1, 3] = 1
    model.Df[1, 4] = 0
    model.Df[1, 5] = 1
    model.Df[1, 6] = 0
    model.Df[1, 7] = 1

    # Constraint functions
    model.g = ConstraintList()
    model.g.add(expr = model.x[0] ** 2 + model.x[1] ** 2 - 1 <= 0)
    model.g.add(expr = model.x[2] ** 2 + model.x[3] ** 2 - 1 <= 0)
    model.g.add(expr = (model.x[4] - 2) ** 2 + (model.x[5] - 5) ** 2 - 10 <= 0)
    model.g.add(expr = (model.x[6] - 3) ** 2 + (model.x[7] - 8) ** 2 - 10 <= 0)

    # Derivatives of constraint functions
    model.Dg = Expression(range(q), range(n+m))
    model.Dg[0, 0] = 2 * model.x[0]
    model.Dg[0, 1] = 2 * model.x[1]
    model.Dg[0, 2] = 0
    model.Dg[0, 3] = 0
    model.Dg[0, 4] = 0
    model.Dg[0, 5] = 0
    model.Dg[0, 6] = 0
    model.Dg[0, 7] = 0

    model.Dg[1, 0] = 0
    model.Dg[1, 1] = 0
    model.Dg[1, 2] = 2 * model.x[2]
    model.Dg[1, 3] = 2 * model.x[3]
    model.Dg[1, 4] = 0
    model.Dg[1, 5] = 0
    model.Dg[1, 6] = 0
    model.Dg[1, 7] = 0

    model.Dg[2, 0] = 0
    model.Dg[2, 1] = 0
    model.Dg[2, 2] = 0
    model.Dg[2, 3] = 0
    model.Dg[2, 4] = 2 * (model.x[4] - 2)
    model.Dg[2, 5] = 2 * (model.x[5] - 5)
    model.Dg[2, 6] = 0
    model.Dg[2, 7] = 0

    model.Dg[3, 0] = 0
    model.Dg[3, 1] = 0
    model.Dg[3, 2] = 0
    model.Dg[3, 3] = 0
    model.Dg[3, 4] = 0
    model.Dg[3, 5] = 0
    model.Dg[3, 6] = 2 * (model.x[6] - 3)
    model.Dg[3, 7] = 2 * (model.x[7] - 8)

    # Solver
    model_parameters.solver = 'gurobi'
    model_parameters.patch_solver = 'ipopt'

    # Initial enclosure
    model_parameters.L = np.array([[-3,5]])
    model_parameters.U = np.array([[12,22]])

    return model, model_parameters