import numpy as np
from solver.subsolvers import *

def snia_boxes_init(m, int_lb, int_ub, num_splits):
    """
    Initializes the boxes to search for a new integer assignment for the fixed boxes approach

    Parameters
    ----------
    m : int
        number of integer variables
    int_lb : numpy array
        lower bounds for the integer variables
    int_ub : numpy array
        upper bounds for the integer variables
    num_splits : int
        number of splits

    Returns
    -------
    boxesLB : numpy array
        lower bounds of all boxes (with each lower bound being a numpy array itself)
    boxesUB : numpy array
        upper bounds of all boxes (with each upper bound being a numpy array itself)
    int_num : int
        number of possible integer assignments
    current_int : numpy array
        array that indicates the position of the last integer assignment that has been computed for each box (which is 0 after the initialization)

    """

    # Extracting integer bounds
    width = int_ub - int_lb
    int_num = np.prod(width + 1)
    boxesLB = np.array([int_lb])
    boxesUB = np.array([int_ub])

    for i in range(num_splits):
        max_width, k = np.max(width), np.argmax(width)
        if max_width < 0.75:
            break
        
        offset = np.zeros(m)
        offset[k] = np.ceil(max_width / 2 - 0.25)
        
        boxesLB = np.row_stack((boxesLB, boxesLB + offset))
        boxesUB = np.row_stack((boxesUB - offset, boxesUB))
        
        width[k] = np.floor(max_width / 2 + 0.25)

    current_int = np.zeros(boxesLB.shape[0], dtype=int)

    return boxesLB, boxesUB, int_num, current_int


def snia_boxes_fixed(build_model, params, ids, X, boxesLB, boxesUB, int_num, current_int, a, r, OFFSET ,GUESS=False):
    """
    Performs the search for a new integer assignment based on the fixed boxes approach

    Parameters
    ----------
    build_model : function
        function that returns the pyomo model of the optimization problem and the model paramaters
    params : array
        parameters for the model of the optimization problem
    ids : numpy array
        point from where to start the search for an update point
    X : numpy array
        direction in which to search for an update point
    boxesLB : numpy array
        lower bounds of all boxes (with each lower bound being a numpy array itself)
    boxesUB : numpy array
        upper bounds of all boxes (with each upper bound being a numpy array itself)
    int_num : int
        number of possible integer assignments
    current_int : numpy array
        array that indicates the position of the last integer assignment that has been computed for each box
    a : numpy array
        parameter for GSUP when using the buest guess based on continuous relaxation
    r : numpy array
        parameter for GSUP when using the buest guess based on continuous relaxation
    OFFSET : float
        tolerance for comparisons and the evaluation of inequalities
    GUESS : boolean
        indicates whether to first try to compute a new integer assignment based on continuous relaxation

    Returns
    -------
    sol_x : numpy array
        part of an optimal solution of SUP corresponding to the x variable of the original optimization problem
    flag : int
        value that indicates why the procedure terminated
    allSolved : boolean
        indicates whether all integer assignments have been visited
    current_int : numpy array
        updated array that indicates the position of the last integer assignment that has been computed for each box

    """

    # Searches a new integer assignment within predefined boxes

    # Load the model parameters
    model, model_param = build_model(params)
    n = model_param.n
    m = model_param.m
    q = model_param.q

    # Initialization
    intArray = np.row_stack((np.row_stack([entry.x_int for entry in ids]), X[:, n:n+m]))
    allSolved = False
    flag = -1

    # Check if all integers are already solved
    if not intArray.shape[0] < int_num:
        allSolved = True
        sol_x = np.array([model.x[i].lb for i in model.x])
        flag = 0
        return sol_x, flag, allSolved, current_int

    # Buest guess based on relaxation
    if GUESS:
        gsup_x, _, _ = GSUP(build_model, params, a, r)
        x_int = np.round(gsup_x[n:n+m])
        sol_x,err = FeasibilityCheckPatch(build_model, params, x_int)
        if err < OFFSET:
            flag = 1
        else:
            flag = 0
        if not np.any(np.all(intArray == x_int, axis = 1)):
            return sol_x, flag, allSolved, current_int

    # Find next integer assignment
    boxes_num = boxesLB.shape[0]
    count = np.zeros(boxes_num)

    for i in range(boxes_num):
        count[i] = np.sum(np.all((boxesLB[i,:] <= intArray) & (intArray <= boxesUB[i,:]), axis=1))

    index = np.argmin(count)
    int_lb = boxesLB[index, :]
    int_ub = boxesUB[index, :]
    width = int_ub - int_lb
    width_prod = np.prod(np.triu(np.ones((m, m)),1) * (width + 1) + np.tril(np.ones((m, m))), axis=1)

    for i in range(current_int[index], int(np.prod(width + 1))):
        x_int = int_lb + np.mod(np.floor(i / width_prod), width + 1)
        current_int[index] = i
        if not np.any(np.all(intArray == x_int, axis = 1)):
            sol_x,err = FeasibilityCheckPatch(build_model, params, x_int)
            if err < OFFSET:
                flag = 1
            else:
                flag = 0
            return sol_x, flag, allSolved, current_int