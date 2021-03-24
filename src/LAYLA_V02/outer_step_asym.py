# -*- coding: utf-8 -*-
"""
inner loop optimiser for asymmetric laminates
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import random
import numpy as np

sys.path.append(r'C:\BELLA_and_LAYLA')
from src.guidelines.ten_percent_rule import calc_penalty_10_pc
from src.guidelines.balance import calc_penalty_bal
from src.guidelines.ipo_oopo import calc_penalty_ipo_oopo_ss
from src.LAYLA_V02.objectives import calc_obj_multi_ss
from src.LAYLA_V02.objectives import objectives
from src.LAYLA_V02.beam_search import beam_search
from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

def outer_step_asym(
        parameters,
        constraints,
        targets,
        lampam_assumed,
        lampam_weightings,
        n_plies_in_groups,
        levels_in_groups,
        n_groups,
        cummul_mom_areas,
        delta_lampams,
        mat_prop=None,
        not_constraints=None):
    '''
    performs an inner loop for asymmetric and unbalanced laminates

    OUTPUTS

    - outer_step_result: instance of the class OuterStepResults

    INPUTS

    - parameters: input parameters for the tuning of the algorithm
    - constraints: set of constraints
    - targets: target lamination parameters and ply counts
    - n_plies: total number of plies
    - cummul_mom_areas: cummulated ply moments of areas
    - delta_lampams: ply partial lamination parameters
    - lampam_assumed: the assumption of the remaining partial
    lamination parameters
    - lampam_weightings: lamination parameter weightings at each search level
    - levels_in_groups: indices of the plies for each ply group optimisation
    - n_groups: number of groups
    - mat_prop: material properties of the laminae
    - not_constraints: design guidelines that should not be satisfied
    '''
    outer_step_result = OuterStepResults()

    lampam_current = sum(lampam_assumed) # initial lamination parameters
    n_plies_per_angle = np.zeros(constraints.n_set_of_angles, float)
    ss_top = np.array([], dtype='int16')
    ss_bot = np.array([], dtype='int16')

    # details # do not consider
    if not_constraints is not None and not_constraints.rule_10_percent:
        random_for_10 = random.randint(0, 4)
    else:
        random_for_10 = None

    for inner_step in range(n_groups):

#        print('inner_step', inner_step)

        if inner_step < n_groups - 1: # not last ply group
            lampam_current -= lampam_assumed[inner_step]

            result = beam_search(
                levels=levels_in_groups[inner_step],
                lampam_current=lampam_current,
                lampam_weightings=lampam_weightings,
                group_size=n_plies_in_groups[inner_step],
                targets=targets,
                parameters=parameters,
                constraints=constraints,
                n_plies_per_angle=n_plies_per_angle,
                cummul_mom_areas=cummul_mom_areas,
                delta_lampams=delta_lampams,
                last_group=False,
                mat_prop=mat_prop,
                not_constraints=not_constraints,
                random_for_10=random_for_10,
                ss_top=ss_top,
                ss_bot=ss_bot)

            lampam_current = result.lampam_best
            n_plies_per_angle = result.ply_counts
#            outer_step_result.n_obj_func_calls += result.n_obj_func_calls

            ss_bot = np.hstack((ss_bot, result.ss_bot_best))
            ss_top = np.hstack((result.ss_top_best, ss_top))
#            print('result.ss_top_best', result.ss_top_best)
#            print('result.ss_bot_best', result.ss_bot_best)
#            print('ss_top', ss_top)
#            print('ss_bot', ss_bot)

        elif inner_step == n_groups - 1: # last ply group

            lampam_current -= lampam_assumed[inner_step]

            result = beam_search(
                levels=levels_in_groups[inner_step],
                lampam_current=lampam_current,
                lampam_weightings=lampam_weightings,
                group_size=n_plies_in_groups[inner_step],
                targets=targets,
                parameters=parameters,
                constraints=constraints,
                n_plies_per_angle=n_plies_per_angle,
                cummul_mom_areas=cummul_mom_areas,
                delta_lampams=delta_lampams,
                last_group=True,
                mat_prop=mat_prop,
                not_constraints=not_constraints,
                random_for_10=random_for_10,
                ss_top=ss_top,
                ss_bot=ss_bot)

            outer_step_result.lampam_best = result.lampam_best
#            outer_step_result.n_obj_func_calls += result.n_obj_func_calls
            outer_step_result.n_designs_last_level = result.n_designs_last_level
            outer_step_result.n_designs_repaired = result.n_designs_repaired
            outer_step_result.n_designs_repaired_unique \
            = result.n_designs_repaired_unique
            outer_step_result.ss_best = result.ss_best

    obj_no_const = objectives(
        outer_step_result.lampam_best,
        targets=targets,
        lampam_weightings=lampam_weightings[-1],
        constraints=constraints,
        parameters=parameters,
        mat_prop=mat_prop)

    n_plies_per_angle = np.zeros((constraints.n_set_of_angles,), float)
    for ind in range(outer_step_result.ss_best.size):
        index = constraints.ind_angles_dict[outer_step_result.ss_best[ind]]
        n_plies_per_angle[index] += 1

    # calculation the penalties for the in-plane and out-of-plane
    # orthotropy requirements based on lamination parameters
    penalty_ipo_lampam, penalty_oopo = calc_penalty_ipo_oopo_ss(
        outer_step_result.lampam_best,
        constraints=constraints,
        parameters=parameters)
#    print('penalty_ipo_lampam', penalty_ipo_lampam)
#    print('penalty_oopo', penalty_oopo)

    # calculation the penalties for the in-plane orthotropy
    # requirements based on ply counts
    penalty_ipo_pc = 0
    if constraints.ipo and parameters.penalty_bal_switch:
        penalty_ipo_pc = calc_penalty_bal(
            n_plies_per_angle,
            constraints)
#    print('penalty_ipo_pc', penalty_ipo_pc)

    penalty_10 = 0
    if constraints.rule_10_percent:
        penalty_10 = calc_penalty_10_pc(n_plies_per_angle, constraints)

    penalty_bal_ipo = max(penalty_ipo_pc, penalty_ipo_lampam)

#    print('obj_no_const', obj_no_const)
#    print('penalty_10', penalty_10)
#    print('penalty_ipo_lampam', penalty_ipo_lampam)
#    print('penalty_ipo_pc', penalty_ipo_pc)
#    print('penalty_oopo', penalty_oopo)

    # calculation of the bounds
    outer_step_result.obj_const = calc_obj_multi_ss(
        objective=obj_no_const,
        penalty_10=penalty_10,
        penalty_bal_ipo=penalty_bal_ipo,
        penalty_oopo=penalty_oopo,
        coeff_10=parameters.coeff_10,
        coeff_bal_ipo=parameters.coeff_bal_ipo,
        coeff_oopo=parameters.coeff_oopo)

    # if repair failed
    if ((constraints.rule_10_percent and penalty_10)\
    or (constraints.ipo and penalty_ipo_lampam \
        and constraints.penalty_ipo_switch) \
    or (constraints.ipo and penalty_ipo_pc \
        and constraints.penalty_bal_switch)):
         outer_step_result.obj_const = 1e10

#    print('outer_step_result.ss_best', result.ss_best)
#    print('outer_step_result.lampam_best', outer_step_result.lampam_best)
#    print('outer_step_result.obj_const', outer_step_result.obj_const)

    return outer_step_result

class OuterStepResults():
    " An object for storing the results of a outer step in LAYLA"
    def __init__(self):
        "Initialise the results of a outer step in LAYLA"
        # solution stacking sequence
        self.ss_best = None
        # solution lamination parameters
        self.lampam_best = None
        # solution ply counts in each fibre direction
        self.ply_counts = None
        # solution constrained objective function
        self.obj_const = None
#        # number of objective function calls during outer step
#        self.n_obj_func_calls = 0
        # number of nodes at the last level of the search tree
        self.n_designs_last_level = 0
        # number of repaired nodes
        self.n_designs_repaired = 0
        # number of unique repaired nodes
        self.n_designs_repaired_unique = 0
