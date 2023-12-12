import numpy as np
import time
import warnings
from solver.update_bounds import *
from solver.subsolvers import *
from solver.snia_functions import *

class IdsEntry:
    def __init__(self, x_int, L, active_state, nd_candidates, eff_candidates):
        self.x_int = x_int
        self.L = L
        self.S = active_state
        self.N = nd_candidates
        self.E = eff_candidates

def f_eval(model, x):
    """
    Evaluates the objective functions of the model at a certain point

    Parameters
    ----------
    model : pyomo model
        a model respresenting the current optimization problem
    x : numpy array
        point where to evaluate the objective functions

    Returns
    -------
    value_obj : numpy array
        an array consisting of the objective function values

    """

    dim = len(x)
    #vars = list(model.component_objects(Var))
    for i in range(dim):
        #vars[i].value = x[i]
        model.x[i].value = x[i]
    #value_obj = np.array([obj.expr() for obj in model.component_objects(Objective)])
    value_obj = np.array([model.f[i+1].expr() for i in range(len(model.f))])

    return value_obj

def compute_nds(list):
    """
    Computes the nondominated points of a given list

    Parameters
    ----------
    list : numpy array
        array of all points

    Returns
    -------
    indexlist : numpy array
        array of boolean values indicating which of the points in the list are nondominated (True) or dominated (False)

    """

    size_list = list.shape[0]
    indexlist = np.ones(size_list, dtype=bool)

    for i in range(size_list):
        current_point = list[i,:]
        indexlist[i] = False
        if not np.any(np.all(current_point >= list[indexlist,:], axis = 1)):
            indexlist[i] = True
    
    return indexlist

def HyPaD(build_model, params, L, U, EPSILON, OFFSET):
    """
    Main procedure for the Hybrid Patch Decomposition algorithm (HyPaD)

    Parameters
    ----------
    build_model : function
        function that returns the pyomo model of the optimization problem and the model paramaters
    params : array
        parameters for the model of the optimization problem
    L : numpy array
        initial array of lower bounds for the enclosure (with each lower bound being a numpy array itself)
    U : numpy array
        initial array of upper bounds for the enclosure (with each upper bound being a numpy array itself)
    EPSILON : float
        tolerance for the width of the enclosure
    OFFSET : float
        tolerance for comparisons and the evaluation of inequalities

    Returns
    -------
    L : numpy array
        array of lower bounds for the enclosure (with each lower bound being a numpy array itself)
    U : numpy array
        array of upper bounds for the enclosure (with each upper bound being a numpy array itself)
    N : numpy array
        array of all potentially nondominated points computed by the algorithm
    ids : array
        integer data structure containing the data computed on the patch level
    it : int
        number of iterations
    exitflag : int
        value that indicates why the algorithm terminated

    """
    
    # Load model parameters
    model,model_parameters = build_model(params)
    n = model_parameters.n
    m = model_parameters.m
    p = model_parameters.p

    # Initialization phase    
    N = np.empty((0,p))
    E = np.empty((0,n+m))
    X = np.empty((0,n+m))
    ids = []
    it = 0
    done = False
    exitflag = -3
    TIME_LIMIT = 3600
    ALL_INTEGERS_VISITED = False
    boxesLB,boxesUB,int_num,current_int = snia_boxes_init(m, np.array([model.x[i].lb for i in range(n,n+m)]), np.array([model.x[i].ub for i in range(n,n+m)]), 4)

    # Test a first solution
    l = np.min(L, axis = 0)
    r = np.max(U, axis = 0) - l
    gsup_x, _, _ = GSUP(build_model, params, l, r)
    gsup_x_int = np.round(gsup_x[n:n+m])
    x, err = FeasibilityCheckPatch(build_model, params, gsup_x_int)
    if err < OFFSET:
        eff_candidates = WSNLP(build_model, params,gsup_x_int, np.eye(p))
        nd_candidates = np.array([f_eval(model,eff_candidates[i,:]) for i in range(eff_candidates.shape[0])])
        z = np.min(nd_candidates, axis = 0)
        N = np.row_stack([N, nd_candidates])
        E = np.row_stack([E, eff_candidates])
        for k in range(nd_candidates.shape[0]):
            U = updateLUB3(U,nd_candidates[k,:])
        ids.append(IdsEntry(gsup_x_int,np.array([z-OFFSET]),True,nd_candidates,eff_candidates))
    else:
        X = np.row_stack([X, x])
    
    # Main Loop
    start_time = time.time()
    while (time.time() - start_time < TIME_LIMIT) and not done:
        it += 1
        size_L = L.shape[0]
        finished_bounds = 0
        rsup_error = True
        for l in L:
            if (time.time() - start_time > TIME_LIMIT):
                break  # Recognize time limit early
            U_temp = np.all(l < (U - EPSILON + OFFSET), axis=1)    
            if np.any(U_temp):
                directions = U[U_temp,:] - l
                u_index = np.argmax(np.min(directions, axis=1))
                r = directions[u_index,:]
                u = l + r
                rsup_x, rsup_eta, rsup_flag = RSUP(build_model, params, np.vstack((X, E)), l, r) 
                if rsup_flag < 0:
                    warnings.warn("RSUP terminated with non-optimal solution", category=UserWarning)
                else:
                    rsup_x_int = np.round(rsup_x[n:n+m])
                    rsup_error = False
                    if np.all(rsup_eta > u + OFFSET):
                        # Try to improve the lower bound set using continuous relaxation
                        warnings.warn("RSUP returned eta > u", category=UserWarning)
                        gsup_x, gsup_t, gsup_flag = GSUP(build_model, params, l, r)
                        if gsup_flag < 0:
                            warnings.warn("GSUP terminated with non-optimal solution", category=UserWarning)
                        else:
                            L = updateLLB3(L, l + gsup_t * r)
                    else:
                        L = updateLLB3(L, rsup_eta)
                    if ~ALL_INTEGERS_VISITED:
                        x, err = FeasibilityCheckPatch(build_model, params, rsup_x_int)
                        if err < OFFSET:
                            if any(ids):
                                ids_index = next((i for i, entry in enumerate(ids) if np.all(entry.x_int == rsup_x_int)), None)
                            else:
                                ids_index = None
                            if ids_index == None:
                                # InitIDS: Initialize entry of ids
                                eff_candidates = WSNLP(build_model, params,rsup_x_int, np.eye(p))
                                nd_candidates = np.array([f_eval(model,eff_candidates[i,:]) for i in range(eff_candidates.shape[0])])
                                z = np.min(nd_candidates, axis = 0)
                                N = np.row_stack([N, nd_candidates])
                                E = np.row_stack([E, eff_candidates])
                                for k in range(nd_candidates.shape[0]):
                                    U = updateLUB3(U,nd_candidates[k,:])
                                ids.append(IdsEntry(rsup_x_int,np.array([z-OFFSET]),True,nd_candidates,eff_candidates))
                            else:
                                if not ids[ids_index].S:
                                    # Select another active patch since the current one is inactive
                                    ids_index = next((i for i, entry in enumerate(ids) if entry.S), None)
                                if not ids_index == None:
                                    # UpdateIDS: Update entry of IDS
                                    patch_L = ids[ids_index].L
                                    patch_done = True
                                    for patch_l in patch_L:
                                        if (time.time() - start_time > TIME_LIMIT):
                                            break  # Recognize time limit early
                                        patch_U_temp = np.all(patch_l < (U - EPSILON + OFFSET), axis=1)
                                        
                                        if np.any(patch_U_temp):
                                            patch_done = False
                                            patch_directions = U[patch_U_temp,:] - patch_l
                                            patch_u_index = np.argmax(np.min(patch_directions, axis=1))
                                            patch_r = patch_directions[patch_u_index,:]
                                            
                                            sup_x, sup_t, sup_flag = SUP(build_model, params, ids[ids_index].x_int, patch_l, patch_r)
                                            
                                            if sup_flag < 0:
                                                warnings.warn("SUP terminated with non-optimal solution", category=UserWarning)
                                            else:
                                                y_lub = f_eval(model,sup_x)
                                                y_llb = patch_l + sup_t * patch_r
                                                U = updateLUB3(U, y_lub)
                                                ids[ids_index].L = updateLLB3(ids[ids_index].L, y_llb)
                                                ids[ids_index].N = np.row_stack([ids[ids_index].N, y_lub])
                                                ids[ids_index].E = np.row_stack([ids[ids_index].E, sup_x])
                                                N = np.row_stack([N, y_lub])
                                                E = np.row_stack([E, sup_x])
                                    if patch_done:
                                        ids[ids_index].S = False
                                else:
                                    # SNIA: Compute new integer assignment
                                    #TODO: Implement more SNIA modes
                                    snia_x,snia_flag,ALL_INTEGERS_VISITED,current_int = snia_boxes_fixed(build_model, params, ids, X, boxesLB, boxesUB, int_num, current_int, l, r, OFFSET)
                                    if snia_flag > 0:
                                        # InitIDS: Initialize entry of ids
                                        snia_x_int = snia_x[n:n+m]
                                        eff_candidates = WSNLP(build_model, params,snia_x_int, np.eye(p))
                                        nd_candidates = np.array([f_eval(model,eff_candidates[i,:]) for i in range(eff_candidates.shape[0])])
                                        z = np.min(nd_candidates, axis = 0)
                                        N = np.row_stack([N, nd_candidates])
                                        E = np.row_stack([E, eff_candidates])
                                        for k in range(nd_candidates.shape[0]):
                                            U = updateLUB3(U,nd_candidates[k,:])
                                        ids.append(IdsEntry(snia_x_int,np.array([z-OFFSET]),True,nd_candidates,eff_candidates))
                                        #print('INFO: SNIA computed a new feasible integer assignment.')
                                    else:
                                        if ALL_INTEGERS_VISITED:
                                            # Terminate HyPaD and compute global lower bound set based on patch level lower bounds
                                            print('END:  All integer assigments visited! Switching solver mode.')
                                            L = np.vstack([entry.L for entry in ids])
                                            L_indexlist = compute_nds(L)
                                            L = L[L_indexlist,:]
                                            exitflag = 0
                                            done = True
                                            return L, U, N, ids, it, exitflag
                                        else:
                                            if snia_flag < 0:
                                                print('INFO: SNIA computed an invalid integer assignment.')
                                                pass
                                            else:
                                                X = np.row_stack([X, snia_x])
                                                print('INFO: SNIA computed a new infeasible integer assignment.')
                        else:
                            X = np.row_stack([X, x])
                indexlist = [i for i, entry in enumerate(ids) if entry.S]
            else:
                finished_bounds += 1
        if finished_bounds > size_L-1:
            print("END:  Finished! Width is <= EPSILON.")
            exitflag = 1
            done = True
            break
        if rsup_error:
            print("ERR:  Algorithm terminated due to RSUP issues.")
            exitflag = -2
            done = True
            break
        if not np.any(indexlist):
            #print('INFO: Iteration ended with all patches inactive');
            pass
    if time.time() - start_time > TIME_LIMIT:
        exitflag = -1
    elif done:
        exitflag = 1
    
    return L, U, N, ids, it, exitflag