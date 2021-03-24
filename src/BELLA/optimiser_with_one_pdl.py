# -*- coding: utf-8 -*-
"""
Optimisation of a composite laminate design based on an input py drop layout
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')

from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

from src.BELLA.results import BELLA_ResultsOnePdl
from src.BELLA.objectives import calc_obj_each_panel
from src.BELLA.objectives import calc_obj_multi_panel
from src.CLA.lampam_functions import calc_lampam
from src.BELLA.moments_of_areas import calc_mom_of_areas2
from src.BELLA.lampam_matrix import calc_delta_lampams2
from src.guidelines.ipo_oopo import calc_penalty_oopo_ss
from src.guidelines.ipo_oopo import calc_penalty_ipo
from src.guidelines.contiguity import calc_penalty_contig_mp
from src.guidelines.disorientation import calc_number_violations_diso_mp
from src.guidelines.ten_percent_rule import calc_penalty_10_pc
from src.guidelines.ten_percent_rule import calc_ply_counts
from src.guidelines.ten_percent_rule import calc_penalty_10_ss
from src.guidelines.ply_drop_spacing  import calc_penalty_spacing
from src.BELLA.pruning import pruning_diso_contig_damtol
from src.RELAY.repair_mp import repair_mp
#from src.guidelines.one_stack import check_ply_drop_rules
from src.BELLA.format_pdl import extend_after_guide_based_blending
#from src.divers.pretty_print import print_list_ss

def BELLA_optimiser_one_pdl(
        multipanel, parameters, obj_func_param, constraints, ply_order,
        mom_areas_plus, delta_lampams, pdl, mat=0):
    "Ply angle optimisation with beam search (BELLA step 3)"

    results = BELLA_ResultsOnePdl()
#    print_list_ss(pdl)

    # to normalise the disorientation and contiguity constraints' penalties
    norm_diso_contig = np.array(
        [panel.n_plies for panel in multipanel.reduced.panels])

    # ply indices in the order in which plies are optimised
    indices = ply_order[-1]
    n_plies_to_optimise = indices.size

    lampam_matrix = calc_delta_lampams2(
        multipanel, constraints, delta_lampams, pdl, n_plies_to_optimise)

    (cummul_areas, cummul_first_mom_areas, cummul_sec_mom_areas) \
    = calc_mom_of_areas2(
        multipanel, constraints, mom_areas_plus, pdl, n_plies_to_optimise)

    # lampam_weightings: lamination parameter weightings at each level of the
    # search (used in the objective function calculation)
    lampam_weightings = np.zeros(
        (n_plies_to_optimise, multipanel.reduced.n_panels, 12))
    for ind_ply in range(n_plies_to_optimise):
        for ind_panel, panel in enumerate(multipanel.reduced.panels):
            lampam_weightings[ind_ply, ind_panel, 0:4] \
            = panel.lampam_weightings[0:4] \
            * cummul_areas[ind_panel, ind_ply]

            lampam_weightings[ind_ply, ind_panel, 4:8] \
            = panel.lampam_weightings[4:8] \
            * cummul_first_mom_areas[ind_panel, ind_ply]

            lampam_weightings[ind_ply, ind_panel, 8:12] \
            = panel.lampam_weightings[8:12] \
            * cummul_sec_mom_areas[ind_panel, ind_ply]

            if not np.isclose(lampam_weightings[ind_ply, ind_panel],
                              np.zeros((12,), float)).all():
                lampam_weightings[ind_ply, ind_panel] \
                /= np.sum(lampam_weightings[ind_ply, ind_panel])

    # Stacking sequences of the active nodes
    ss_bot_tab = [[np.array([], dtype=int) \
                   for ii in range(multipanel.reduced.n_panels)]]
    if not constraints.sym:
        ss_top_tab = [[np.array([], dtype=int) \
                       for ii in range(multipanel.reduced.n_panels)]]

    # Lamination parameters of the active nodes
    lampam_tab = np.zeros((1, multipanel.reduced.n_panels, 12), float)

    # Estimate function values of the active nodes
    obj_constraints_tab = np.zeros((1,), dtype=float)
    obj_no_constraints_tab = np.zeros((1, multipanel.reduced.n_panels,), float)

    # Ply counts in each fibre orientation of the active nodes
    n_plies_per_angle_tab = np.zeros((
        1, multipanel.reduced.n_panels, constraints.n_set_of_angles),
        dtype='float16')

    n_obj_func_calls = 0 # Number of objective function calls

    ss_final = []
    sst_final = []
    pdl_final = []
    penalty_diso_final = []
    penalty_contig_final = []
    penalty_bal_ipo_final = []
    penalty_10_final = []

    for level in indices:
# =============================================================================
# preparation of the exploration of a level during beam search #
# =============================================================================
        last_level = level == indices[-1]
        # print('level in the search tree', level)

        n_nodes = obj_constraints_tab.size

        for node in range(n_nodes):
# =============================================================================
# selection of first active node to be branched #
# =============================================================================
            mother_ss_bot = ss_bot_tab.pop(0)
            if not constraints.sym:
                mother_ss_top = ss_top_tab.pop(0)

            mother_lampam = lampam_tab[0]

            mother_n_plies_per_angle = n_plies_per_angle_tab[0]

            lampam_tab = np.delete(lampam_tab, np.s_[0], axis=0)

            n_plies_per_angle_tab = np.delete(
                n_plies_per_angle_tab, np.s_[0], axis=0)

            obj_constraints_tab = np.delete(obj_constraints_tab, np.s_[0])
            obj_no_constraints_tab = np.delete(
                obj_no_constraints_tab, np.s_[0], axis=0)

#            print('at the beginning')
#            print('mother_ss_bot', mother_ss_bot)
#            print()
#            print('mother_ss_top', mother_ss_top)

# =============================================================================
# branching #
# =============================================================================
            child_ss = np.copy(constraints.set_of_angles)

# =============================================================================
# pruning for damage tolerance, disorientation and contiguity in guide panel
# =============================================================================
            if constraints.sym:
                ss = pruning_diso_contig_damtol(
                    child_ss=child_ss,
                    mother_ss_bot=mother_ss_bot[multipanel.reduced.ind_thick],
                    level=level,
                    constraints=constraints,
                    has_middle_ply=multipanel.has_middle_ply,
                    n_plies_to_optimise=n_plies_to_optimise)
            else:
                ss = pruning_diso_contig_damtol(
                    child_ss=child_ss,
                    mother_ss_top=mother_ss_top[multipanel.reduced.ind_thick],
                    mother_ss_bot=mother_ss_bot[multipanel.reduced.ind_thick],
                    level=level,
                    constraints=constraints,
                    n_plies_to_optimise=n_plies_to_optimise)
            if ss is None: # branch totally pruned
                continue # go to next step
#            print('ss after pruning for design guidelines', ss.T)

# =============================================================================
# calculation of ply counts #
# =============================================================================
            n_plies_per_angle = np.zeros((
                ss.size, multipanel.reduced.n_panels,
                constraints.n_set_of_angles), dtype=float)
            for ind_angle in range(ss.size):
                n_plies_per_angle[ind_angle] = np.copy(
                    mother_n_plies_per_angle)

                for ind_panel in range(multipanel.reduced.n_panels):
                    if pdl[ind_panel, level] >= 0:
                        index = constraints.ind_angles_dict[ss[ind_angle]]
                        if last_level and multipanel.has_middle_ply:
                            n_plies_per_angle[
                                ind_angle, ind_panel, index] += 1/2
                        else:
                            n_plies_per_angle[
                                ind_angle, ind_panel, index] += 1

            n_plies_per_angle_tab = np.vstack((
                n_plies_per_angle_tab, n_plies_per_angle))
#            print('n_plies_per_angle', n_plies_per_angle)

# =============================================================================
# calculation of lamination parameters #
# =============================================================================
            for ind_angle in range(ss.size):
                lampam = np.copy(mother_lampam)
                for ind_panel in range(multipanel.reduced.n_panels):
                    lampam[ind_panel] += lampam_matrix[
                        ind_panel,
                        constraints.ind_angles_dict[ss[ind_angle]], level]
                lampam_tab = np.vstack(
                    (lampam_tab,
                     lampam.reshape(1, multipanel.reduced.n_panels, 12)))

                n_obj_func_calls += 1

# =============================================================================
# computation of stacking sequences #
# =============================================================================
            if constraints.sym:
                for ind_angle in range(ss.size):
                    ss_bot = list(mother_ss_bot)
                    for ind_panel in range(multipanel.reduced.n_panels):
                        if pdl[ind_panel, level] >= 0:
                            ss_bot[ind_panel] = np.hstack((
                                ss_bot[ind_panel], ss[ind_angle]))
                    ss_bot_tab.append(ss_bot)

            else:
                for ind_angle in range(ss.size):

                    ss_top = list(mother_ss_top)
                    ss_bot = list(mother_ss_bot)

                    if level % 2 == 0:
                        for ind_panel in range(multipanel.reduced.n_panels):
                            if pdl[ind_panel, level] >= 0:
                                ss_bot[ind_panel] = np.hstack((
                                    ss_bot[ind_panel], ss[ind_angle]))
                    else:
                        for ind_panel in range(multipanel.reduced.n_panels):
                            if pdl[ind_panel, level] >= 0:
                                ss_top[ind_panel] = np.hstack((
                                    ss[ind_angle], ss_top[ind_panel]))

                    ss_top_tab.append(ss_top)
                    ss_bot_tab.append(ss_bot)

#            print('ss_bot_tab', ss_bot_tab)
#            print()

# =============================================================================
# computation of penalties and objective function values
# =============================================================================
            for ind_angle, ind_angle2 \
            in zip(range(-ss.size, 0, 1), range(ss.size)):

                # calulation of the objective/estimate function values
                obj_no_constraints = calc_obj_each_panel(
                    multipanel=multipanel,
                    lampam=lampam_tab[ind_angle],
                    mat=mat,
                    obj_func_param=obj_func_param,
                    lampam_weightings=lampam_weightings[level])
                obj_no_constraints_tab = np.vstack((
                    obj_no_constraints_tab, obj_no_constraints))


                # calulation of the penalties for balance
                penalty_bal_ipo = np.zeros((
                    multipanel.reduced.n_panels,), float)

                # calulation of the penalties for the 10% rule
                penalty_10 = np.zeros((
                    multipanel.reduced.n_panels,), float)

                # calulation of the penalties for the disorientation rule
                penalty_diso = np.zeros((
                    multipanel.reduced.n_panels,), float)

                # calulation of the penalties for the contiguity rule
                penalty_contig = np.zeros((
                    multipanel.reduced.n_panels,), float)

                # calulation of the estimate/objective functions
                obj_constraints = calc_obj_multi_panel(
                    objective=obj_no_constraints,
                    actual_panel_weightings=\
                    multipanel.reduced.actual_panel_weightings,
                    penalty_diso=penalty_diso,
                    penalty_contig=penalty_contig,
                    penalty_10=penalty_10,
                    penalty_bal_ipo=penalty_bal_ipo,
                    coeff_diso=obj_func_param.coeff_diso,
                    coeff_contig=obj_func_param.coeff_contig,
                    coeff_10=obj_func_param.coeff_10,
                    coeff_bal_ipo=obj_func_param.coeff_bal_ipo)
                obj_constraints_tab = np.hstack((
                    obj_constraints_tab, obj_constraints))

#                print('mother_ss_top', mother_ss_top)
#                print('mother_ss_bot', mother_ss_bot)
#                print('obj_no_constraints', obj_no_constraints)
#                print('obj_constraints_tab', obj_constraints_tab)

# =============================================================================
# local pruning
# =============================================================================
            n_local_nodes = ss.size
            if last_level:
                node_limit = parameters.local_node_limit_final
            else:
                node_limit = parameters.local_node_limit
            n_excess_nodes = n_local_nodes - node_limit
            if level != n_plies_to_optimise - 1 and n_excess_nodes > 0:
                obj_constraints_tab_to_del = np.copy(
                    obj_constraints_tab)[-n_local_nodes:]
                to_del = []
                for counter in range(n_excess_nodes):
                    ind_max = np.argmax(obj_constraints_tab_to_del)
                    obj_constraints_tab_to_del[ind_max] = -6666
                    to_del.append(
                        ind_max + obj_constraints_tab.size - n_local_nodes)
                for ind_del in sorted(to_del, reverse=True):
                    del ss_bot_tab[ind_del]
                    if not constraints.sym:
                        del ss_top_tab[ind_del]
                lampam_tab = np.delete(lampam_tab, np.s_[to_del], axis=0)
                n_plies_per_angle_tab = np.delete(
                    n_plies_per_angle_tab, np.s_[to_del], axis=0)
                obj_constraints_tab = np.delete(
                    obj_constraints_tab, np.s_[to_del])
                obj_no_constraints_tab = np.delete(
                    obj_no_constraints_tab, np.s_[to_del], axis=0)

#            print('mother_ss_top', mother_ss_top)
#            print('mother_ss_bot', mother_ss_bot)
#            print('len(ss_top_tab)', len(ss_top_tab))
#            print('ss after local pruning', ss.T)
#            print('obj_constraints_tab', obj_constraints_tab[-ss.size:])
#            print('lampam_tab.shape', lampam_tab.shape)
# =============================================================================
# global pruning
# =============================================================================
        if obj_constraints_tab.size == 0:
            raise Exception("""
beam search with no solutions !
    -Constraints too strict or optimiser parameters not allowing for
    sufficient design space exploration""")
        if last_level:
            node_limit = parameters.global_node_limit_final
        else:
            node_limit = parameters.global_node_limit
        if obj_constraints_tab.size > node_limit:
            obj_constraints_tab_to_del = np.copy(obj_constraints_tab)
            for counter in range(node_limit):
                ind_min = np.argmin(obj_constraints_tab_to_del)
                obj_constraints_tab_to_del[ind_min] = 6666
            to_keep = obj_constraints_tab_to_del == 6666
            to_del = np.invert(to_keep).astype(int)
            to_del = [i for i, x in enumerate(to_del) if x]
            for ind_del in sorted(to_del, reverse=True):
                del ss_bot_tab[ind_del]
                if not constraints.sym:
                    del ss_top_tab[ind_del]
            lampam_tab = np.delete(lampam_tab, np.s_[to_del], axis=0)
            n_plies_per_angle_tab = np.delete(
                n_plies_per_angle_tab, np.s_[to_del], axis=0)
            obj_constraints_tab = np.delete(
                obj_constraints_tab, np.s_[to_del])
            obj_no_constraints_tab = np.delete(
                obj_no_constraints_tab, np.s_[to_del], axis=0)
#        print('size after global pruning', obj_constraints_tab.size)


# =============================================================================
# repair with RELAY
# =============================================================================
    for ind_angle in range(obj_constraints_tab.size):

        # generate full stacking
        stack = []
        if constraints.sym:
            for ind_panel, panel \
            in enumerate(multipanel.reduced.panels):
                stack.append(np.array((), dtype=int))
                stack[ind_panel] = np.hstack((
                    ss_bot_tab[ind_angle][ind_panel],
                    np.flip(ss_bot_tab[ind_angle][ind_panel],
                            axis=0)
                    )).astype('int16')

                if panel.has_middle_ply:
                    stack[ind_panel] = np.delete(
                        stack[ind_panel],
                        np.s_[panel.middle_ply_index],
                        axis=0)

                if stack[ind_panel].size != panel.n_plies:
                    raise Exception("Wrong ply count")
        else:
            for ind_panel, panel \
            in enumerate(multipanel.reduced.panels):
                stack.append(np.array((), dtype=int))
                stack[ind_panel] = np.hstack((
                    ss_bot_tab[ind_angle][ind_panel],
                    ss_top_tab[ind_angle][ind_panel]
                    )).astype('int16')
                if stack[ind_panel].size != panel.n_plies:
                    raise Exception("Wrong ply count")
    #                    print('stack', stack)

        # laminate repair with RELAY
        results_repair = repair_mp(
            multipanel, stack, constraints, parameters,
            obj_func_param, reduced_pdl=pdl, mat=mat)
        results.n_designs_last_level += 1
        # successful repair
        if results_repair[0]:
            results.n_designs_after_ss_ref_repair += 1
            results.n_designs_after_thick_to_thin += 1
            results.n_designs_after_thin_to_thick += 1
        # repair reference panel unsuccessful
        elif results_repair[4] == 1:
            pass
        # repair thick-to-thin unsuccessful
        elif results_repair[4] == 2:
            results.n_designs_after_ss_ref_repair += 1
        # repair thin-to-thick unsuccessful
        else:
            results.n_designs_after_ss_ref_repair += 1
            results.n_designs_after_thick_to_thin += 1

        # possible changes during repair
        lampam_tab[ind_angle] = calc_lampam(
            results_repair[1], constraints)
        ss_final.append(results_repair[1])
        sst_final.append(results_repair[2])
        pdl_final.append(results_repair[3])
        ply_counts = calc_ply_counts(
            multipanel, results_repair[1], constraints)
        n_plies_per_angle_tab[ind_angle] = ply_counts

        # penalty for disorientation
        n_diso = calc_number_violations_diso_mp(results_repair[1], constraints)
        if constraints.diso and n_diso.any():
            penalty_diso = n_diso / norm_diso_contig
        else:
            penalty_diso = np.zeros((multipanel.reduced.n_panels,), float)

        # penalty for contiguity
        n_contig = calc_penalty_contig_mp(results_repair[1], constraints)
        if constraints.contig and n_contig.any():
            penalty_contig = n_contig / norm_diso_contig
        else:
            penalty_contig = np.zeros((multipanel.reduced.n_panels,), float)

        # penalty for 10% rule
        if constraints.rule_10_percent and constraints.rule_10_Abdalla:
            penalty_10 = calc_penalty_10_ss(results_repair[1],
                                            constraints,
                                            LPs=lampam_tab[ind_angle], mp=True)
        else:
            penalty_10 = calc_penalty_10_pc(ply_counts, constraints)

        # penalty for balance
        penalty_bal_ipo = np.zeros((multipanel.reduced.n_panels,), float)

        penalty_diso_final.append(penalty_diso)
        penalty_contig_final.append(penalty_contig)
        penalty_bal_ipo_final.append(penalty_bal_ipo)
        penalty_10_final.append(penalty_10)

        # calulation of the objective/estimate function values
        if not results_repair[0]:
                obj_no_constraints = 1e10 * np.ones((
                    multipanel.reduced.n_panels,), float)
        else:
            obj_no_constraints = calc_obj_each_panel(
                multipanel=multipanel,
                lampam=lampam_tab[ind_angle],
                mat=mat,
                obj_func_param=obj_func_param,
                lampam_weightings=lampam_weightings[level])
        obj_no_constraints_tab[ind_angle] = obj_no_constraints

        if not (obj_no_constraints != 1e10).all():
            obj_constraints = 1e10
        else:
            obj_constraints = calc_obj_multi_panel(
                objective=obj_no_constraints,
                actual_panel_weightings=\
                multipanel.reduced.actual_panel_weightings,
                penalty_diso=penalty_diso,
                penalty_contig=penalty_contig,
                penalty_10=penalty_10,
                penalty_bal_ipo=penalty_bal_ipo,
                coeff_diso=obj_func_param.coeff_diso,
                coeff_contig=obj_func_param.coeff_contig,
                coeff_10=obj_func_param.coeff_10,
                coeff_bal_ipo=obj_func_param.coeff_bal_ipo)
        obj_constraints_tab[ind_angle] = obj_constraints

# =============================================================================
# select best result
# =============================================================================
    ss_bot_tab = np.copy(ss_final)

    if parameters.save_success_rate:
        to_keep = np.array([ind for ind in range(len(obj_constraints_tab))\
                            if obj_constraints_tab[ind] < 1e10])
        to_del = np.array([ind for ind in range(len(obj_constraints_tab))\
                           if ind not in to_keep])
        if to_keep.size:
            ss_final = [ss_final[ind] for ind in to_keep]
            sst_final = [sst_final[ind] for ind in to_keep]
            pdl_final = [pdl_final[ind] for ind in to_keep]
            penalty_diso_final = [penalty_diso_final[ind] for ind in to_keep]
            penalty_contig_final = [
                penalty_contig_final[ind] for ind in to_keep]
            penalty_bal_ipo_final = [
                penalty_bal_ipo_final[ind] for ind in to_keep]
            penalty_10_final = [penalty_10_final[ind] for ind in to_keep]

        for ind in to_del[::-1]:
            lampam_tab = np.delete(lampam_tab, np.s_[ind], axis=0)
            obj_constraints_tab = np.delete(
                obj_constraints_tab, np.s_[ind], axis=0)
            obj_no_constraints_tab = np.delete(
                obj_no_constraints_tab, np.s_[ind], axis=0)

        results.n_designs_repaired_unique = np.unique(
            sst_final, axis=0).shape[0]

    if len(obj_constraints_tab) == 0:
        print('No laminate design solution with this initial ply drop layout')
        return None

    index = np.argmin(obj_constraints_tab)
    ss = ss_bot_tab[index]
    lampam = lampam_tab[index]
    n_plies_per_angle = n_plies_per_angle_tab[index]
    obj_no_constraints = obj_no_constraints_tab[index]
    obj_constraints = obj_constraints_tab[index]
# =============================================================================
# Check the results
# =============================================================================
    # test for the ply counts
    for ind_panel, panel in enumerate(multipanel.reduced.panels):
        if ss[ind_panel].size != panel.n_plies:
            raise Exception("""
Wrong ply counts in the laminate.""")

    # test for the partial lamination parameters
    ss = list(ss)
    lampam_test = calc_lampam(ss, constraints)
    if not abs(lampam - lampam_test).all() < 1e-13:
        print_lampam(lampam[0], lampam_test[0])
        print_ss(ss[0])
        raise Exception("""
Lamination parameters not matching the stacking sequences.""")

    # test for the objective function value
    obj_no_constraints_test = calc_obj_each_panel(
        multipanel=multipanel,
        lampam=lampam,
        lampam_weightings=lampam_weightings[level],
        mat=mat,
        obj_func_param=obj_func_param)

    if abs(obj_no_constraints_test - obj_no_constraints > 1e-10).any():
        print('obj_no_constraints_test', obj_no_constraints_test)
        print('obj_no_constraints', obj_no_constraints)
        raise Exception("""
Objective function value not matching the stacking sequences.""")


    # disorientaion - penalty used in blending steps 4.2 and 4.3
    n_diso = calc_number_violations_diso_mp(ss, constraints)
    if constraints.diso and n_diso.any():
        penalty_diso = n_diso / norm_diso_contig
    else:
        penalty_diso = np.zeros((multipanel.reduced.n_panels,))

    # contiguity - penalty used in blending steps 4.2 and 4.3
    n_contig = calc_penalty_contig_mp(ss, constraints)
    if constraints.contig and n_contig.any():
        penalty_contig = n_contig / norm_diso_contig
    else:
        penalty_contig = np.zeros((multipanel.reduced.n_panels,))

    # 10% rule - no penalty used in blending steps 4.2 and 4.3
    if constraints.rule_10_percent and constraints.rule_10_Abdalla:
        penalty_10 = calc_penalty_10_ss(ss, constraints, lampam, mp=True)
    else:
        penalty_10 = calc_penalty_10_pc(n_plies_per_angle, constraints)

    # balance
    penalty_bal_ipo = np.zeros((multipanel.reduced.n_panels,))

#    print()
#    print_ss(ss[2])
#    print('obj_no_constraints', obj_no_constraints)
#    print(penalty_diso)
#    print(penalty_contig)
#    print(penalty_bal_ipo)
#    print(penalty_10)

    obj_constraints_test = calc_obj_multi_panel(
        objective=obj_no_constraints,
        actual_panel_weightings=multipanel.reduced.actual_panel_weightings,
        penalty_diso=penalty_diso,
        penalty_contig=penalty_contig,
        penalty_10=penalty_10,
        penalty_bal_ipo=penalty_bal_ipo,
        coeff_diso=obj_func_param.coeff_diso,
        coeff_contig=obj_func_param.coeff_contig,
        coeff_10=obj_func_param.coeff_10,
        coeff_bal_ipo=obj_func_param.coeff_bal_ipo)

    if abs(obj_constraints_test - obj_constraints > 1e-10).any():
        print('obj_constraints_test', obj_constraints_test)
        print('obj_constraints', obj_constraints)
        pass
        raise Exception("""
Objective function value not matching the stacking sequences.""")


    ### calculates penalties for all panels

    # theses penalties are not used in caculations, just used to show the
    # degrees of feasibility of the designed laminates in the result report

    # balance
    penalty_bal_ipo = calc_penalty_ipo(lampam)

    # out-of-plane orthotropy
    penalty_oopo = calc_penalty_oopo_ss(lampam, constraints=constraints)

    # penalty_spacing
    penalty_spacing = calc_penalty_spacing(
        pdl=pdl_final[index],
        multipanel=multipanel,
        constraints=constraints,
        on_blending_strip=True)

# =============================================================================
# return the results
# =============================================================================
    results.ss = extend_after_guide_based_blending(multipanel, ss)
    results.lampam = extend_after_guide_based_blending(multipanel, lampam)
    results.n_plies_per_angle = extend_after_guide_based_blending(
        multipanel, n_plies_per_angle)
    results.n_obj_func_calls = n_obj_func_calls
    results.obj_constraints = obj_constraints
    results.obj_no_constraints = extend_after_guide_based_blending(
        multipanel, obj_no_constraints)
    results.penalty_diso = extend_after_guide_based_blending(
        multipanel, penalty_diso)
    results.penalty_contig = extend_after_guide_based_blending(
        multipanel, penalty_contig)
    results.penalty_10 = extend_after_guide_based_blending(
        multipanel, penalty_10)
    results.penalty_bal_ipo = extend_after_guide_based_blending(
        multipanel, penalty_bal_ipo)
    results.penalty_oopo = extend_after_guide_based_blending(
        multipanel, penalty_oopo)
    results.penalty_spacing = penalty_spacing
    results.n_diso = extend_after_guide_based_blending(multipanel, n_diso)
    results.n_contig = extend_after_guide_based_blending(multipanel, n_contig)
    results.sst = extend_after_guide_based_blending(
        multipanel, sst_final[index])
    results.pdl = extend_after_guide_based_blending(
        multipanel, pdl_final[index])
    return results
