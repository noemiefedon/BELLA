# -*- coding: utf-8 -*-
"""
Function arranging the outer loops of the optimisation technique LAYLA

First outer loop:
    Layerwise approach and partial lamination parameters of unknown plies
    assumed as 0 (quasi-isotropicity)

The refinement loops:
    Layerwise approach and partial lamination parameters of unknown plies
    calculated based on the orientation of plies at same locations taken from
    a previously determined stacking sequence
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA_and_LAYLA')

from src.LAYLA_V02.divide_laminate_sym import divide_laminate_sym
from src.LAYLA_V02.divide_laminate_asym import divide_laminate_asym
from src.LAYLA_V02.ply_order import calc_ply_order, calc_levels
from src.LAYLA_V02.moment_of_areas import calc_mom_of_areas
from src.LAYLA_V02.results import LAYLA_Results

from src.LAYLA_V02.outer_step_sym import outer_step_sym
from src.LAYLA_V02.outer_step_asym import outer_step_asym

from src.guidelines.one_stack import check_lay_up_rules
from src.CLA.lampam_functions import calc_lampam_from_delta_lp_matrix
from src.CLA.lampam_functions import calc_delta_lampam
#from src.divers.pretty_print import print_lampam, print_ss, print_list_ssI


def LAYLA_optimiser(parameters, constraints, targets,
                    mat_prop=None, not_constraints=None):
    """
    performs the retrieval of stacking sequences from lamination parameter
    targets

     OUTPUTS

    - LAYLA_results: results of the optimisation

    INPUTS

    - parameters: parameters of the optimiser
    - constraints: lay-up design guidelines
    - targets: target lamination parameters and ply counts
    - mat_prop: material properties
    """

    # details # do not consider
    if not_constraints is not None and not constraints.sym:
        raise Exception("""
Set of constraints that must not be satisfied only accounted for symmetric
laminates""")

    # filter the lamination parameters for orthotropy requirements
    targets.filter_target_lampams(constraints)

    if constraints.sym:
        if targets.n_plies % 2 == 1:
            middle_ply = int((targets.n_plies + 1) / 2)
        else:
            middle_ply = 0

    ply_order = calc_ply_order(constraints, targets)

    # division of the laminates in groups of plies
    if constraints.sym:
        n_plies_in_groups, pos_first_ply_groups, n_groups = \
        divide_laminate_sym(parameters, targets, step=0)
    else:
        n_plies_in_groups, pos_first_ply_groups, n_groups = \
        divide_laminate_asym(parameters, targets, step=0)

    levels_in_groups = calc_levels(ply_order, n_plies_in_groups, n_groups)

    # mom_areas: signed ply moments of areas
    # cummul_mom_areas: cummulated positive ply moments of areas
    # group_mom_areas: ply-group positive moments of areas
    mom_areas, cummul_mom_areas, group_mom_areas = calc_mom_of_areas(
        constraints, targets, ply_order, n_plies_in_groups)

    # lampam_weightings: lamination parameter weightings at each level of the
    # search (used in the objective function calculation)
    lampam_weightings = parameters.lampam_weightings_final * np.hstack((
        np.matlib.repmat(cummul_mom_areas[:, 0][:, np.newaxis], 1, 4),
        np.matlib.repmat(cummul_mom_areas[:, 1][:, np.newaxis], 1, 4),
        np.matlib.repmat(cummul_mom_areas[:, 2][:, np.newaxis], 1, 4)))
    lampam_weightings = np.array([
        lampam_weightings[ind]/sum(lampam_weightings[ind]) \
        for ind in range(lampam_weightings.shape[0])])

    # calculation of ply partial lamination parameters
    if constraints.sym:
        delta_lampams = np.empty((targets.n_plies // 2 + targets.n_plies % 2,
                                   constraints.n_set_of_angles, 12), float)
        for ind in range(delta_lampams.shape[0]):
            delta_lampams[
                ind, :, 0:4] = mom_areas[ind, 0] * constraints.cos_sin
            delta_lampams[ind, :, 4:8] = 0
            delta_lampams[
                ind, :, 8:12] = mom_areas[ind, 2] * constraints.cos_sin
    else:
        delta_lampams = np.empty((
            targets.n_plies, constraints.n_set_of_angles, 12), float)
        for ind in range(lampam_weightings.shape[0]):
            delta_lampams[
                ind, :, 0:4] = mom_areas[ind, 0] * constraints.cos_sin
            delta_lampams[
                ind, :, 4:8] = mom_areas[ind, 1] * constraints.cos_sin
            delta_lampams[
                ind, :, 8:12] = mom_areas[ind, 2] * constraints.cos_sin

    # asummed lamination parameters for ech ply groups
    lampam_assumed = np.zeros((n_groups, 12), float)

    results = LAYLA_Results(parameters, targets)

    for n_outer_step in range(parameters.n_outer_step):
        print('n_outer_step', n_outer_step)

        if constraints.sym:
            outputs = outer_step_sym(
                cummul_mom_areas=cummul_mom_areas,
                delta_lampams=delta_lampams,
                lampam_weightings=lampam_weightings,
                parameters=parameters,
                constraints=constraints,
                targets=targets,
                lampam_assumed=lampam_assumed,
                n_plies_in_groups=n_plies_in_groups,
                levels_in_groups=levels_in_groups,
                middle_ply=middle_ply,
                n_groups=n_groups,
                mat_prop=mat_prop,
                not_constraints=not_constraints)
        else:
            outputs = outer_step_asym(
                cummul_mom_areas=cummul_mom_areas,
                delta_lampams=delta_lampams,
                lampam_weightings=lampam_weightings,
                parameters=parameters,
                constraints=constraints,
                targets=targets,
                lampam_assumed=lampam_assumed,
                n_plies_in_groups=n_plies_in_groups,
                levels_in_groups=levels_in_groups,
                n_groups=n_groups,
                mat_prop=mat_prop,
                not_constraints=not_constraints)

        lampam_check = calc_lampam_from_delta_lp_matrix(
            outputs.ss_best, constraints, delta_lampams)
        check_lay_up_rules(outputs.ss_best, constraints)

        if sum(abs(outputs.lampam_best - lampam_check)) > 1e-10:
            raise Exception('Lamination parameters not matching lay-up')

        if outputs.ss_best.size != targets.n_plies:
            raise Exception('Stacking sequence with incorrect ply count')

        results.ss_tab[n_outer_step] = outputs.ss_best
        results.lampam_tab_tab[n_outer_step] = outputs.lampam_best
        results.obj_tab[n_outer_step] = outputs.obj_const
#        results.n_obj_func_calls_tab[n_outer_step] = outputs.n_obj_func_calls
        results.n_designs_last_level_tab[
            n_outer_step] = outputs.n_designs_last_level
        results.n_designs_repaired_tab[
            n_outer_step] = outputs.n_designs_repaired
        results.n_designs_repaired_unique_tab[
            n_outer_step] = outputs.n_designs_repaired_unique

        # if the stacking sequence is the same or good enough, exit the loop
        if n_outer_step == 0:
            if outputs.obj_const < 1e-10:
                break
        elif np.allclose(outputs.ss_best, results.ss_tab[n_outer_step -1]) \
        or outputs.obj_const < 1e-10:
                break

        # Repartitioning of the laminate into groups ?
        if (n_outer_step != parameters.n_outer_step - 1) \
        and parameters.group_size_max[n_outer_step] \
        != parameters.group_size_max[n_outer_step + 1]:

            # division of the laminates in groups of plies
            if constraints.sym:
                n_plies_in_groups, pos_first_ply_groups, n_groups = \
                divide_laminate_sym(parameters, targets, step=n_outer_step + 1)
            else:
                n_plies_in_groups, pos_first_ply_groups, n_groups = \
                divide_laminate_asym(
                    parameters, targets, step=n_outer_step + 1)

            mom_areas, cummul_mom_areas, group_mom_areas = calc_mom_of_areas(
                constraints, targets, ply_order, n_plies_in_groups)

            levels_in_groups = calc_levels(
                ply_order, n_plies_in_groups, n_groups)

        # Updating the table of assumptions for the next step
        lampam_assumed = np.zeros((n_groups, 12), float)
        for ind_group in range(n_groups):
            for ind_ply in range(n_plies_in_groups[ind_group]):
                lampam_assumed[ind_group] += delta_lampams[
                    levels_in_groups[ind_group][ind_ply],
                    constraints.ind_angles_dict[
                        outputs.ss_best[
                            levels_in_groups[ind_group][ind_ply]]]]

    if np.isnan(results.obj_tab).all():
        raise Exception('No successful repair during lay-up optimisation')
        return results

    ind_mini = np.nanargmin(results.obj_tab)
    results.number_of_outer_steps_performed = n_outer_step + 1
    results.n_outer_step_best_solution = ind_mini + 1
    results.objective = results.obj_tab[ind_mini]
    results.ss = results.ss_tab[ind_mini]
    results.lampam = results.lampam_tab_tab[ind_mini]
    results.completed = True

    return results