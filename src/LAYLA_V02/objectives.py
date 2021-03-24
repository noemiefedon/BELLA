# -*- coding: utf-8 -*-
"""
objective functions
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def objectives(
        lampam,
        targets,
        lampam_weightings,
        constraints,
        parameters,
        mat_prop=None):
    '''
    returns the objective function evaluated for a lamination
    parameter

    INPUTS

    - lampam: lamiation parameters
    - targets.lampam: target lamiation parameters
    - lampam_weightings: set of coefficients to be applied to the square terms
    of the norm 2, specific for each step of the algorithm
    - constraints: design and manufacturing constraints
    - parameters: optimiser parameters
    - mat_prop: material properties of the laminae
    '''
    if lampam.ndim == 1:

        if parameters.type_obj_func == 2:
            return sum(lampam_weightings*(lampam - targets.lampam)**2)
        elif parameters.type_obj_func == 1:
            return sum(lampam_weightings*abs(lampam - targets.lampam))
        else:
            raise Exception('TYpe of objective function not recognised')

    result = np.zeros((lampam.shape[0],), float)

    for ind in range(lampam.shape[0]):
        result[ind] = objectives(
            lampam=lampam[ind],
            targets=targets,
            lampam_weightings=lampam_weightings,
            constraints=constraints,
            parameters=parameters,
            mat_prop=mat_prop)
    return result


def calc_obj_multi_ss(
        objective,
        penalty_10=0,
        penalty_bal_ipo=0,
        penalty_oopo=0,
        coeff_10=0,
        coeff_bal_ipo=0,
        coeff_oopo=0):
    """
    calculates objective function for single-panel structures considering the
    penalties for the design and manufacturing guidelines

    INPUTS

    - objective: objective function value of each panel with no
    consideration of the penalties for the design rules
    - penalty_10: penalty of each panel for the 10% rule constraint
    - penalty_bal_ipo: penalty for in-plane orthotropy/balance
    - penalty_oopo: penalty of each panel for out-of-plane orthotropy
    - coeff_10: weight of the penalty for the 10% rule
    - coeff_bal_ipo: weight of the penalty for in-plane orthotropy/balance
    - coeff_oopo: weight of the penalty for out-of-plane orthotropy
    """
    return objective * (1 + coeff_10 * penalty_10) \
           * (1 + coeff_bal_ipo * penalty_bal_ipo) \
           * (1 + coeff_oopo * penalty_oopo)
