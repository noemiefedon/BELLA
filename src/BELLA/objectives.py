# -*- coding: utf-8 -*-
"""
objective functions

- calc_obj_one_panel
    calculates objective function value for a single panel during beam_search
    with no consideration of the penalties for the design rules

- calc_obj_each_panel
    calculates objective function of multi-panel structures during beam search
    with no consideration of the penalties for the design rules

- calc_obj_multi_panel
    calculates objective function of a multi_panel structure considering the
    penalties for the design and manufacturing guidelines

- calc_unconst_obj_multi_panel_A
    calculates the unconstrained in-plane objective function values for
    multi-panel structures

- calc_unconst_obj_multi_panel_D
    calculates the unconstrained out-of-plane objective function values for
    multi-panel structures
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
from src.buckling.buckling import buckling_factor

def calc_obj_one_panel(
        lampam,
        lampam_target,
        lampam_weightings=np.array([])):
    """
    calculates objective function value for a single panel during beam_search
    with no consideration of the penalties for the design rules

    INPUTS

    - lampam: panel lamination parameters
    - lampam_target: target lamination parameters
    - lampam_weightings: weights of the lamination parameters in the objective
    function
    - obj_func_param: objective function parameters
    """
    return lampam_weightings@((lampam - lampam_target)**2)


def calc_obj_each_panel(
        multipanel, lampam, obj_func_param, mat=0, lampam_weightings=[]):
    """
    calculates objective function of multi-panel structures during beam search
    with no consideration of the penalties for the design rules

    INPUTS

    - multipanel: multi-panel class instance
    - lampam: group partial lamination parameters
    - lampam_weightings: weightings for each lamination parameter
    - mat: material properties of the laminae
    - obj_func_param: objective function parameters
    """
    objectives = np.zeros((multipanel.reduced.n_panels), dtype=float)

    for ind_panel, panel in enumerate(multipanel.reduced.panels):

#        print('lampam[ind_panel]', lampam[ind_panel])
#        print('panel.lampam_target', panel.lampam_target)
#        print('lampam_weightings[ind_panel]', lampam_weightings[ind_panel])

        objectives[ind_panel] = calc_obj_one_panel(
            lampam=lampam[ind_panel],
            lampam_target=panel.lampam_target,
            lampam_weightings=lampam_weightings[ind_panel])
        
    return objectives


def calc_obj_multi_panel(
        objective,
        actual_panel_weightings,
        penalty_diso=0,
        penalty_contig=0,
        penalty_10=0,
        penalty_bal_ipo=0,
        penalty_oopo=0,
        penalty_weight=0,
        coeff_diso=0,
        coeff_contig=0,
        coeff_10=0,
        coeff_bal_ipo=0,
        coeff_oopo=0,
        coeff_weight=0,
        with_Nones=False):
    """
    calculates objective function of a multi_panel structure considering the
    penalties for the design and manufacturing guidelines

    INPUTS

    - objective: objective function value of each panel with no
    consideration of the penalties for the design rules
    - actual_panel_weightings: weightings of the different panels in the
    objective function
    - penalty_diso: penalty of each panel for the disorientation constraint
    - penalty_contig: penalty of each panel for the contiguity constraint
    - penalty_10: penalty of each panel for the 10% rule constraint
    - penalty_bal_ipo: penalty of each panel for in-plane orthotropy/balance
    - penalty_oopo: penalty of each panel for out-of-plane orthotropy
    - penalty_weight: penalty of each panel for weight
    - coeff_diso: weight of the penalty for the disorientation rule
    - coeff_contig: weight of the penalty for the contiguity rule
    - coeff_10: weight of the penalty for the 10% rule
    - coeff_bal_ipo: weight of the penalty for in-plane orthotropy/balance
    - coeff_oopo: weight of the penalty for out-of-plane orthotropy
    - coeff_weight: weight of the penalty for weight
    """
#    print(penalty_weight)
#    print(penalty_diso)
#    print(penalty_contig)
#    print(penalty_bal_ipo)
#    print(penalty_oopo)
#    print(penalty_10)
#    print(objective)
    if not with_Nones:
        return sum(actual_panel_weightings * objective \
                   * (1 + coeff_diso * penalty_diso) \
                   * (1 + coeff_contig * penalty_contig) \
                   * (1 + coeff_10 * penalty_10) \
                   * (1 + coeff_bal_ipo * penalty_bal_ipo) \
                   * (1 + coeff_oopo * penalty_oopo) \
                   * (1 + coeff_weight * penalty_weight))
    my_sum = 0
    for ind in range(actual_panel_weightings.size):
        if objective[ind] is not None:
            to_add = actual_panel_weightings[ind] * objective[ind]
            
            if penalty_diso is not None:
                if (isinstance(penalty_diso, list) \
                    and len(penalty_diso) > 1) or penalty_diso.size > 1:
                    to_add *= (1 + coeff_diso * penalty_diso[ind])
                else:
                    to_add *= (1 + coeff_diso * penalty_diso)
            
            if penalty_contig is not None:
                if (isinstance(penalty_contig, list) \
                    and len(penalty_contig) > 1) or  penalty_contig.size > 1:
                    to_add *= (1 + coeff_contig * penalty_contig[ind])
                else:
                    to_add *= (1 + coeff_contig * penalty_contig)

            if penalty_10 is not None:
                if (isinstance(penalty_10, list) and len(penalty_10) > 1) \
                or  penalty_10.size > 1:
                    to_add *= (1 + coeff_10 * penalty_10[ind])
                else:
                    to_add *= (1 + coeff_10 * penalty_10)

            if penalty_oopo is not None:
                if (isinstance(penalty_oopo, list) \
                    and len(penalty_oopo) > 1) or  penalty_oopo.size > 1:
                    to_add *= (1 + coeff_oopo * penalty_oopo[ind])
                else:
                    to_add *= (1 + coeff_oopo * penalty_oopo)

            if penalty_bal_ipo is not None:
                to_add *= (1 + coeff_bal_ipo * penalty_bal_ipo[ind])

            if penalty_weight is not None:
                to_add *= (1 + coeff_weight * penalty_weight[ind])

            my_sum += to_add
    return my_sum

def calc_unconst_obj_multi_panel_A(
        multipanel, lampamA, obj_func_param, inner_step=-1, mat=0):
    """
    calculates unconstrained in-plane objective function for multi-panel
    structures

    INPUTS

    - multipanel: multi-panel class instance
    - lampamA: in-plane partial lamination parameters
    - inner_step: inner loop step number
    - mat: material properties of the laminae
    - obj_func_param: objective function parameters
    """
    objectives = np.zeros((multipanel.reduced.n_panels), dtype=float)
    if inner_step == -1:
        for ind_reduced_panel in range(multipanel.reduced.n_panels):
            objectives[ind_reduced_panel] = calc_obj_one_panel(
                lampam=lampamA[ind_reduced_panel],
                lampam_target=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_target[:4],
                lampam_weightings=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_weightings[:4])
    else:
        for ind_reduced_panel in range(multipanel.reduced.n_panels):
            objectives[ind_reduced_panel] = calc_obj_one_panel(
                lampam=lampamA[ind_reduced_panel],
                lampam_target=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_target[:4],
                lampam_weightings=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_weightings[:4])
    return sum(obj_func_param.reduced_actual_panel_weightingsA * objectives)

def calc_unconst_obj_multi_panel_D(
        multipanel, lampamD, obj_func_param, inner_step=-1, mat=0):
    """
    calculates the unconstrained out-of-plane objective function values for
    multi-panel structures

    INPUTS

    - multipanel: multi-panel class instance
    - lampamD: out-of-plane partial lamination parameters
    - inner_step: inner loop step number
    - mat: material properties of the laminae
    - obj_func_param: objective function parameters
    """
    objectives = np.zeros((multipanel.reduced.n_panels), dtype=float)
    if inner_step == -1:
        for ind_reduced_panel in range(multipanel.reduced.n_panels):
            objectives[ind_reduced_panel] = calc_obj_one_panel(
                lampam=lampamD[ind_reduced_panel],
                lampam_target=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_target[8:12],
                lampam_weightings=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_weightings[8:12])
    else:
        for ind_reduced_panel in range(multipanel.reduced.n_panels):
            objectives[ind_reduced_panel] = calc_obj_one_panel(
                lampam=lampamD[ind_reduced_panel],
                lampam_target=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_target[8:12],
                lampam_weightings=multipanel.panels[
                    multipanel.ind_for_reduc[
                        ind_reduced_panel]].lampam_weightings[8:12])
    return sum(obj_func_param.reduced_actual_panel_weightingsD * objectives)

