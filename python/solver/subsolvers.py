import numpy as np
from pyomo.environ import *
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr

def init_enclosure(model, model_parameters, OFFSET):
    """
    Computes an initial enclosure consisting of a single upper and lower bound based on interval arithmetic.

    Parameters
    ----------
    model : pyomo model
        a model respresenting the current optimization problem
    model_parameters : structure
        a structure respresenting the model parameters
    OFFSET : float
        offset for the lower and upper bounds to ensure that the area of interest is contained in the interior of the enclosure

    Returns
    -------
    L : numpy array
        set of lower bounds for the initial enclosure
    U : numpy array
        set of upper bounds for the initial enclosure

    """

    p = model_parameters.p
    L = model_parameters.L
    U = model_parameters.U
    if not (np.any(L) and np.any(U)):
        for i in range(p):
            L[0,i], U[0,i] = compute_bounds_on_expr(model.f[i+1].expr)
    return L - OFFSET, U + OFFSET

def WSNLP(build_model, params, x_int, weights):
    # Load the original model
    model,model_param = build_model(params)
    n = model_param.n
    m = model_param.m

    # Fix the integer variables
    for i in range(m):
        model.x[n+i].fix(x_int[i])

    # Initialization of return values and number of iterations for weighted sum approach to solve the subproblem
    iterations = weights.shape[0]
    sol_x = np.empty((iterations,n+m))

    # Solve the Pyomo model
    solver = SolverFactory(model_param.patch_solver)
    for i in range(iterations):
        model.obj = Objective(expr=sum(weights[i,j] * model.f[j+1].expr for j in range(len(model.f))))
        results = solver.solve(model, tee=False)
        sol_x[i,:] = np.array([value(model.x[i]) for i in range(len(model.x))])
        model.del_component(model.obj)

    return sol_x

def FeasibilityCheckPatch(build_model, params, x_int):
    # Load the original model
    model,model_param = build_model(params)
    n = model_param.n
    m = model_param.m
    q = model_param.q

    # Fix the integer variables
    for i in range(m):
        model.x[n+i].fix(x_int[i])

    # Add the additional variable alpha
    model.add_component('alpha', Var(within=Reals))
    alpha = model.component('alpha')

    # Modify the constraints
    for j in range(q):
        con = model.g[j+1]
        con.set_value(con.body - alpha <= con.upper)
        
    # Objective function
    model.del_component(model.f)
    model.obj = Objective(expr=alpha)

    # Solve the Pyomo model
    solver = SolverFactory(model_param.patch_solver)
    results = solver.solve(model, tee=False)

    # Extract solution
    err = value(alpha)
    sol_x = np.array([value(model.x[i]) for i in range(len(model.x))])
    if not results.solver.termination_condition == TerminationCondition.optimal:
        exitflag = -1
    else:
        exitflag = 0

    return sol_x, err

def SUP(build_model, params, x_int, a, r):
    """
    Performs the search for an update point on the patch level

    Parameters
    ----------
    build_model : function
        function that returns the pyomo model of the optimization problem and the model paramaters
    params : array
        parameters for the model of the optimization problem
    x_int : numpy array
        integer assignment (that describes the current patch problem)
    a : numpy array
        point from where to start the search for an update point
    r : numpy array
        direction in which to search for an update point

    Returns
    -------
    sol_x : numpy array
        part of an optimal solution of SUP corresponding to the x variable of the original optimization problem
    sol_t : float
        optimal value of the optimization problem SUP
    exitflag : int
        value that indicates whether the problem SUP was solved to optimality

    """

    # Load the original model
    model,model_param = build_model(params)
    n = model_param.n
    m = model_param.m
    p = model_param.p

    # Fix the integer variables
    for i in range(m):
        model.x[n+i].fix(x_int[i])

    # Add the additional variable t
    model.add_component('t', Var(within=Reals))
    t = model.component('t')

    # Add the SUP constraints
    for i in range(p):
        obj = model.f[i+1]
        model.add_component('sup_con_'+str(i+1),Constraint(expr=obj.expr - t * r[i] <= a[i]))

    # Objective function
    model.del_component(model.f)
    model.obj = Objective(expr=t)

    # TODO: Add warm start initialization for x variables

    # Solve the Pyomo model
    solver = SolverFactory(model_param.patch_solver)
    results = solver.solve(model, tee=False)

    # Extract solution
    sol_x = np.array([value(model.x[i]) for i in range(len(model.x))])  
    sol_t = value(t)

    if not results.solver.termination_condition == TerminationCondition.optimal:
        exitflag = -1
    else:
        exitflag = 0

    return sol_x, sol_t, exitflag

def RSUP(build_model, params, X, a, r):
    """
    Performs the search for an update point on the global level based on the relaxed problem R(X)

    Parameters
    ----------
    build_model : function
        function that returns the pyomo model of the optimization problem and the model paramaters
    params : array
        parameters for the model of the optimization problem
    X : numpy array
        array of linearization points (with each linearization point being a numpy array itself)
    a : numpy array
        point from where to start the search for an update point
    r : numpy array
        direction in which to search for an update point

    Returns
    -------
    sol_x : numpy array
        part of an optimal solution of RSUP corresponding to the x variable of the relaxed problem R(X)
    sol_eta : float
        part of an optimal solution of RSUP corresponding to the eta variable of the relaxed problem R(X)
    exitflag : int
        value that indicates whether the problem RSUP was solved to optimality

    """
    # Load the original model
    model,model_param = build_model(params)
    n = model_param.n
    m = model_param.m
    p = model_param.p
    q = model_param.q

    # Add the additional variable t
    model.add_component('t', Var(within=Reals))
    t = model.component('t')

    # Add eta to shift objectives to constraints
    model.add_component('eta', Var(range(p),within=Reals))
    eta = model.component('eta')

    # Add the linearized objective and constraint functions
    for k in range(X.shape[0]):
        x_lin = X[k,:]
        for i in range(n+m):
            model.x[i].value = x_lin[i]
        for i in range(p):
            obj = model.f[i+1]
            model.add_component('f_lin_'+str(k)+'_'+str(i+1),Constraint(expr=obj.expr() + sum(model.Df[i,j].expr() * (model.x[j] - x_lin[j]) for j in range(len(model.x))) - eta[i] <= 0))
        for i in range(q):
            con = model.g[i+1]
            model.add_component('g_lin_'+str(k)+'_'+str(i+1),Constraint(expr=con.body() + sum(model.Dg[i,j].expr() * (model.x[j] - x_lin[j]) for j in range(len(model.x))) <= con.upper))

    # Add the SUP constraints
    for i in range(p):
        model.add_component('sup_con_'+str(i+1),Constraint(expr=eta[i] - t * r[i] <= a[i]))

    # Remove old objective and constraint functions
    model.del_component(model.f)
    model.del_component(model.g)
    
    # Objective function
    model.obj = Objective(expr=t)

    # Solve the Pyomo model
    solver = SolverFactory(model_param.solver)
    results = solver.solve(model, tee=False)

    # Extract solution
    sol_x = np.array([value(model.x[i]) for i in range(len(model.x))])
    sol_eta = np.array([value(eta[i]) for i in range(len(eta))])  
    sol_t = value(t)

    if not results.solver.termination_condition == TerminationCondition.optimal:
        exitflag = -1
    else:
        exitflag = 0

    return sol_x, sol_eta, exitflag

def GSUP(build_model, params, a, r):
    """
    Performs the search for an update point on the global level based on the continuous relaxation

    Parameters
    ----------
    build_model : function
        function that returns the pyomo model of the optimization problem and the model paramaters
    params : array
        parameters for the model of the optimization problem
    a : numpy array
        point from where to start the search for an update point
    r : numpy array
        direction in which to search for an update point

    Returns
    -------
    sol_x : numpy array
        part of an optimal solution of GSUP corresponding to the x variable of the original optimization problem
    sol_t : float
        optimal value of the optimization problem GSUP
    exitflag : int
        value that indicates whether the problem GSUP was solved to optimality

    """

    # Load the original model
    model,model_param = build_model(params)
    n = model_param.n
    m = model_param.m
    p = model_param.p

    # Change integer variables to reals
    for j in range(n,n+m):
        model.x[j].domain = Reals

    # Add the additional variable t
    model.add_component('t', Var(within=Reals))
    t = model.component('t')

    # Add the SUP constraints
    for i in range(p):
        obj = model.f[i+1]
        model.add_component('sup_con_'+str(i+1),Constraint(expr=obj.expr - t * r[i] <= a[i]))

    # Objective function
    model.del_component(model.f)
    model.obj = Objective(expr=t)

    # TODO: Add warm start initialization for x variables

    # Solve the Pyomo model
    solver = SolverFactory(model_param.patch_solver)
    results = solver.solve(model, tee=False)

    # Extract solution
    sol_x = np.array([value(model.x[i]) for i in range(len(model.x))])  
    sol_t = value(t)

    if not results.solver.termination_condition == TerminationCondition.optimal:
        exitflag = -1
    else:
        exitflag = 0

    return sol_x, sol_t, exitflag
