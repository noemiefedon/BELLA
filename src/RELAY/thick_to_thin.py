# -*- coding: utf-8 -*-
"""
Functions used for the repair of panels with the thick-to-thin methodology:
    plies of a reference stacking sequence are left unchanged and the
    positions of the ply drops necessary to design the thinner panels are
    re-optimised with the objective to better match panel lamination
    parameter targets while also satisfying design and manufacturing
    constraints

- repair_thick_to_thin
    performs repair of multi-panel structure by modifying the ply drop layout
    with the thick-to-thin methodology

- calc_panel_weigthings
    returns the panel weightings for the objective function at each step of the
    thick-to-thin or thin-to-thick repair

- calc_parameters
    caluclates values usefull during the beam search for thick-to-thin or
    thin-to-thick repair
"""
import sys
import numpy as np
from copy import deepcopy

sys.path.append(r'C:\BELLA')
from src.guidelines.ply_drop_spacing import calc_penalty_spacing
from src.guidelines.ply_drop_spacing import is_same_pdl
from src.guidelines.contiguity import calc_penalty_contig_ss
from src.guidelines.disorientation import calc_n_penalty_diso_ss
# from src.guidelines.ipo_oopo import calc_penalty_oopo_ss
from src.guidelines.ten_percent_rule import calc_penalty_10_pc
from src.guidelines.ten_percent_rule import calc_penalty_10_ss
from src.BELLA.objectives import calc_obj_one_panel
from src.BELLA.objectives import calc_obj_multi_panel
from src.CLA.lampam_functions import calc_lampam
from src.divers.pretty_print import print_list_ss

# if apply_balance_1_by_1 == True
#     - if an angled ply is added/removed from a balance panel, the next ply
#        to be removed/added rectify balance
# otherwise
#     - all panels of the blending strip are enforced to be balance. This can
#       cause issues if the blending strip contains many panels with
#       incremental ply counts
apply_balance_1_by_1 = True

def repair_thick_to_thin(
        reduced_lampam, reduced_sst, reduced_ss,
        multipanel, parameters, obj_func_param, constraints, mat=0):
    """
    performs repair of multi-panel structure by modifying the ply drop layout
    with the thick-to-thin methodology:
        plies of a reference stacking sequence are left unchanged and the
        positions of the ply drops necessary to design the thinner panels are
        re-optimised with the objective to better match panel lamination
        parameter targets while also satisfying design and manufacturing
        constraints
    """
    ### initialisation
    ss_ref = np.copy(reduced_ss[multipanel.reduced.ind_ref])
    pdl_ref = np.copy(reduced_sst[multipanel.reduced.ind_ref])
    n_plies_max = multipanel.reduced.n_plies_in_panels[-1]
#    print_list_ss(reduced_sst)

    n_steps = multipanel.reduced.n_steps_thin
    ind_panel_tab = multipanel.reduced.ind_panel_thin_tab
    new_boundary_tab = multipanel.reduced.new_boundary_thin_tab
#    print('n_steps', n_steps)
#    print('ind_panel_tab', ind_panel_tab)
#    print('new_boundary_tab', new_boundary_tab)

    if n_steps == 0: # no change
        return True, reduced_sst, reduced_lampam, reduced_ss

    ## list of ply drop position indices
    if constraints.sym:
        all_pdl_indices_tab = [[ind for ind, elem in enumerate(reduced_sst[
            multipanel.reduced.ind_ref][:n_plies_max // 2]) if elem == -1]]
    else:
        all_pdl_indices_tab = [[ind for ind, elem in enumerate(reduced_sst[
            multipanel.reduced.ind_ref][:n_plies_max]) if elem == -1]]
#    print('all_pdl_indices_tab', all_pdl_indices_tab)
    last_pdl_index_tab = [None]

    ## ply-drop layouts
    initial_pdl = [None for ind in range(multipanel.reduced.n_panels)]
    initial_pdl[
        multipanel.reduced.ind_ref] = reduced_sst[multipanel.reduced.ind_ref]
    pdls_tab = [initial_pdl]
#    print('initial_pdl', initial_pdl)

    ## lamination parameters
    initial_lampam = [None for ind in range(multipanel.reduced.n_panels)]
    initial_lampam[multipanel.reduced.ind_ref] \
    = reduced_lampam[multipanel.reduced.ind_ref]
    lampam_tab = [initial_lampam]

    ## number of plies in each direction
    n_plies_per_angle_ref = np.zeros(
        (constraints.n_set_of_angles), dtype='float16')
    for index in range(ss_ref.size):
        index = constraints.ind_angles_dict[ss_ref[index]]
        n_plies_per_angle_ref[index] += 1
    initial_n_plies_per_angle = [
        None for ind in range(multipanel.reduced.n_panels)]
    initial_n_plies_per_angle[multipanel.reduced.ind_ref] \
    = np.copy(n_plies_per_angle_ref)
    n_plies_per_angle_tab = [initial_n_plies_per_angle]

    ## penalties
    initial_pdl_diso = [None for ind in range(multipanel.reduced.n_panels)]
    initial_pdl_diso[multipanel.reduced.ind_ref] = 0
    penalty_diso_tab = [initial_pdl_diso]

    # initial_pdl_contig = [None for ind in range(multipanel.reduced.n_panels)]
    # initial_pdl_contig[multipanel.reduced.ind_ref] = 0
    # penalty_contig_tab = [initial_pdl_contig]

    # initial_pdl_oopo = [None for ind in range(multipanel.reduced.n_panels)]
    # initial_pdl_oopo[multipanel.reduced.ind_ref] = 0
    # penalty_oopo_tab = [initial_pdl_oopo]

    # initial_pdl_10 = [None for ind in range(multipanel.reduced.n_panels)]
    # initial_pdl_10[multipanel.reduced.ind_ref] = 0
    # penalty_10_tab = [initial_pdl_10]

    penalty_spacing_tab = [None]

    angle_queue_tab = [np.array((), dtype=int)]

    ## objectives
    obj_no_constraints_tab = [[
        None for ind in range(multipanel.reduced.n_panels)]]
    obj_constraints_tab = np.zeros((1,), dtype=float)

    n_obj_func_calls = 0

    for ind_step in range(n_steps):
#        print('ind_step', ind_step)

        if len(obj_constraints_tab) == 0:
            return False, reduced_sst, reduced_lampam, reduced_ss

#        print('len(obj_constraints_tab)',
#              len(obj_constraints_tab))

        for node in range(len(obj_constraints_tab)):

            ### selection of node to be branched (first in the list)
            mother_pdl = pdls_tab.pop(0)
            mother_all_pdl_indices = all_pdl_indices_tab.pop(0)
            mother_lampam = lampam_tab.pop(0)
            mother_n_plies_per_angle = n_plies_per_angle_tab.pop(0)
            mother_obj_no_constraints = obj_no_constraints_tab.pop(0)
            mother_penalty_diso = penalty_diso_tab.pop(0)
            # mother_penalty_contig = penalty_contig_tab.pop(0)
            mother_angle_queue = angle_queue_tab.pop(0)
            # mother_penalty_oopo = penalty_oopo_tab.pop(0)
            # mother_penalty_10 = penalty_10_tab.pop(0)
            del penalty_spacing_tab[0]
            del last_pdl_index_tab[0]
            obj_constraints_tab = np.delete(obj_constraints_tab, np.s_[0])

            ### branching + pruning for damage tolerance rule and covering rule
            # pd : index of ply deleted from reference stacking sequence
            if constraints.covering:
                if constraints.n_covering == 1:
                    if constraints.sym:
                        child_pd_indices = np.arange(1, n_plies_max // 2)
                    else:
                        child_pd_indices = np.arange(1, n_plies_max - 1)

                elif constraints.n_covering == 2:
                    if constraints.sym:
                        child_pd_indices = np.arange(2, n_plies_max // 2)
                    else:
                        child_pd_indices = np.arange(2, n_plies_max - 2)
            else:
                if constraints.sym:
                    child_pd_indices = np.arange(n_plies_max // 2)
                else:
                    child_pd_indices = np.arange(n_plies_max)

            ### remove duplicates
            child_pd_indices = np.setdiff1d(child_pd_indices,
                                            mother_all_pdl_indices)

            n_tab_nodes = 0
            tab_child_pdl = []
            tab_all_pdl_indices = []
            tab_penalty_spacing = []
            tab_angle_queue = []
            tab_last_pdl_index = []

            for one_pd_index in child_pd_indices:
#                print('one_pd_index', one_pd_index)

                ### pruning for balance
                if constraints.bal:
                    angle = pdl_ref[one_pd_index]

                    if apply_balance_1_by_1:
                        if mother_angle_queue \
                        and angle != -mother_angle_queue[0]:
                            continue
                    else:
                        if mother_angle_queue:
                            if angle != - mother_angle_queue[0]:
                                continue
                        elif (ind_step == n_steps - 1 \
                              or new_boundary_tab[ind_step + 1]) \
                        and angle not in (0, 90):
                            continue

                if constraints.bal:
                    if mother_angle_queue:
                        tab_angle_queue.append([])
                    elif angle not in (0, 90):
                        tab_angle_queue.append([angle])
                    else:
                        tab_angle_queue.append([])
                else:
                    tab_angle_queue.append([])

                ### ply drop layout
                child_pdl = deepcopy(mother_pdl)
                if new_boundary_tab[ind_step]:
                    child_pdl[ind_panel_tab[ind_step]] = np.copy(
                        child_pdl[ind_panel_tab[ind_step] + 1])
                child_pdl[ind_panel_tab[ind_step]][one_pd_index] = - 1
                if constraints.sym:
                    child_pdl[ind_panel_tab[ind_step]][
                        n_plies_max - one_pd_index - 1] = - 1

                ### penalties for the ply-drop layout rule
                penalty_spacing = calc_penalty_spacing(
                    pdl=child_pdl,
                    multipanel=multipanel,
                    obj_func_param=obj_func_param,
                    constraints=constraints,
                    on_blending_strip=True)

                tab_child_pdl.append(child_pdl[:])
                tab_penalty_spacing.append(penalty_spacing)
                new_list = list(mother_all_pdl_indices)
                new_list.append(one_pd_index)
                tab_all_pdl_indices.append(new_list)
                tab_last_pdl_index.append(one_pd_index)

#                print('child_pdl', child_pdl)

            ### local pruning for the ply-drop layout rules
            indices_to_keep = []
            tab_penalty_spacing_for_pruning = np.copy(tab_penalty_spacing)
            if len(tab_penalty_spacing_for_pruning) \
            > parameters.local_node_limit2:

                while len(indices_to_keep) < parameters.local_node_limit2:

                    min_value = min(tab_penalty_spacing_for_pruning)
                    indices_to_add = np.where(
                        tab_penalty_spacing_for_pruning == min_value)[0]
                    for elem in indices_to_add:
                        indices_to_keep.append(elem)
                        tab_penalty_spacing_for_pruning[elem] = 1000

                indices_to_keep.sort()
#                print('indices_to_keep', indices_to_keep)
                tab_child_pdl = [tab_child_pdl[index] \
                                   for index in indices_to_keep]
                tab_penalty_spacing = [tab_penalty_spacing[index] \
                                     for index in indices_to_keep]
                tab_angle_queue = [tab_angle_queue[index] \
                                     for index in indices_to_keep]
                tab_last_pdl_index = [tab_last_pdl_index[index] \
                                     for index in indices_to_keep]
                tab_all_pdl_indices = [tab_all_pdl_indices[index] \
                                      for index in indices_to_keep]

            ### calculations of lay-up penalties and multi-panel objective
            # function values
            tab_child_n_plies_per_angle = []
            tab_child_lampam = []
            # tab_child_penalty_oopo = []
            tab_child_penalty_diso = []
            # tab_child_penalty_contig = []
            tab_child_obj_no_constraints = []
            tab_child_obj_constraints = []
            # tab_child_penalty_10 = []

            for ind_pd in range(len(tab_child_pdl))[::-1]:
                ### calculation of the stacking sequence in the currently
                # optimised panel
                child_ss = np.copy(tab_child_pdl[
                        ind_pd][ind_panel_tab[ind_step]])
                child_ss = child_ss[child_ss != -1]
#                print('child_ss', child_ss.size)
#                print_ss(child_ss[:child_ss.size // 2], 30)

                ### calculation of the number of plies in each direction
                # ~ child_n_plies_per_angle = np.copy(mother_n_plies_per_angle)
                child_n_plies_per_angle = deepcopy(mother_n_plies_per_angle)
                if new_boundary_tab[ind_step]:
                    child_n_plies_per_angle[ind_panel_tab[ind_step]] = np.copy(
                        child_n_plies_per_angle[ind_panel_tab[ind_step] + 1])
                index_pd = tab_last_pdl_index[ind_pd]
                my_angle = pdl_ref[index_pd]
                index = constraints.ind_angles_dict[my_angle]

                child_n_plies_per_angle[ind_panel_tab[ind_step]][index] -= 1
                if constraints.sym and n_plies_max - index_pd - 1 != index_pd:
                    child_n_plies_per_angle[
                        ind_panel_tab[ind_step]][index] -= 1

#                print('tab_all_pdl_indices', tab_all_pdl_indices[ind_pd])
#                print('index_pd', index_pd, 'my_angle', my_angle)
#                print('child_n_plies_per_angle')
#                print(child_n_plies_per_angle)

                ### calculation of penalties for the disorientation constraint
                if constraints.diso:
                    penalty_diso = calc_n_penalty_diso_ss(
                        child_ss, constraints)
                else:
                    penalty_diso = 0
                child_penalty_diso = deepcopy(mother_penalty_diso)
                child_penalty_diso[ind_panel_tab[ind_step]] = penalty_diso
#                print('child_penalty_diso', child_penalty_diso)

                ### calculation of penalties for the contiguity constraint
                if constraints.contig:
                    penalty_contig = calc_penalty_contig_ss(
                        child_ss, constraints)
                    # pruning for contiguity
                    if penalty_contig != 0:
                        # print('contig')
                        del tab_child_pdl[ind_pd]
                        del tab_penalty_spacing[ind_pd]
                        del tab_angle_queue[ind_pd]
                        del tab_all_pdl_indices[ind_pd]
                        del tab_last_pdl_index[ind_pd]
                        continue
                else:
                    penalty_contig = 0
                # child_penalty_contig = deepcopy(mother_penalty_contig)
                # child_penalty_contig[ind_panel_tab[ind_step]] = penalty_contig
#                print('child_penalty_contig', child_penalty_contig)

                ### calculation of lamination parameters
                child_lampam = deepcopy(mother_lampam)
                child_lampam[ind_panel_tab[ind_step]] \
                = calc_lampam(child_ss, constraints)
#                print('child_lampam', child_lampam)

                ### 10% rule
                if constraints.rule_10_percent:
                    if constraints.rule_10_Abdalla:
                        penalty_10 = calc_penalty_10_ss(
                            child_ss,
                            constraints,
                            LPs=child_lampam[ind_panel_tab[ind_step]],
                            mp=False)
                    else:
                        penalty_10 = calc_penalty_10_pc(
                            child_n_plies_per_angle[ind_panel_tab[ind_step]],
                            constraints)
                    # pruning 10% rule
                    if penalty_10 != 0:
                        # print('10')
                        del tab_child_pdl[ind_pd]
                        del tab_penalty_spacing[ind_pd]
                        del tab_angle_queue[ind_pd]
                        del tab_all_pdl_indices[ind_pd]
                        del tab_last_pdl_index[ind_pd]
                        continue
                else:
                    penalty_10 = 0
                # child_penalty_10 = deepcopy(mother_penalty_10)
                # child_penalty_10[ind_panel_tab[ind_step]] = penalty_10
#                print('child_penalty_10', child_penalty_10)

                ### calculation of objective function values
                obj_no_constraints = calc_obj_one_panel(
                    lampam=child_lampam[ind_panel_tab[ind_step]],
                    lampam_target=multipanel.reduced.panels[
                        ind_panel_tab[ind_step]].lampam_target,
                    lampam_weightings=multipanel.reduced.panels[
                        ind_panel_tab[ind_step]].lampam_weightings)
                child_obj_no_constraints = deepcopy(mother_obj_no_constraints)
                child_obj_no_constraints[
                    ind_panel_tab[ind_step]] = obj_no_constraints
#                print('child_obj_no_constraints', child_obj_no_constraints)

                child_obj_constraints = calc_obj_multi_panel(
                    objective=child_obj_no_constraints,
                    actual_panel_weightings=multipanel.reduced.actual_panel_weightings,
                    penalty_diso=child_penalty_diso,
                    penalty_contig=None,
                    penalty_oopo=None,
                    penalty_10=None,
                    penalty_bal_ipo=None,
                    penalty_weight=None,
                    with_Nones=True)
#                print('child_obj_constraints', child_obj_constraints)

                ### saving
                tab_child_n_plies_per_angle.append(child_n_plies_per_angle)
                tab_child_lampam.append(child_lampam)
                # tab_child_penalty_oopo.append(child_penalty_oopo)
                tab_child_penalty_diso.append(child_penalty_diso)
                # tab_child_penalty_contig.append(child_penalty_contig)
                # tab_child_penalty_10.append(child_penalty_10)
                tab_child_obj_no_constraints.append(child_obj_no_constraints)
                tab_child_obj_constraints.append(child_obj_constraints)

                n_obj_func_calls += 1
                n_tab_nodes += 1


            ### local pruning for the other guidelines and stiffness optimality
            indices_to_keep = []
            tab_child_obj_constraints_for_pruning \
            = np.copy(tab_child_obj_constraints)
            if ind_step != n_steps - 1 \
            and len(tab_child_obj_constraints_for_pruning) \
            > parameters.local_node_limit2:

                while len(indices_to_keep) < parameters.local_node_limit2:

                    min_value = min(tab_child_obj_constraints_for_pruning)
                    index_to_add = np.where(
                    tab_child_obj_constraints_for_pruning == min_value)[0][0]
                    indices_to_keep.append(index_to_add)
                    tab_child_obj_constraints_for_pruning[index_to_add] = 1000

                indices_to_keep.sort()
#                print('indices_to_keep', indices_to_keep)
                tab_child_pdl = [tab_child_pdl[index] \
                                 for index in indices_to_keep]
                tab_all_pdl_indices = [tab_all_pdl_indices[index] \
                                       for index in indices_to_keep]
                tab_last_pdl_index = [tab_last_pdl_index[index] \
                                      for index in indices_to_keep]
                tab_penalty_spacing = [tab_penalty_spacing[index] \
                                   for index in indices_to_keep]
                tab_angle_queue = [tab_angle_queue[index] \
                                   for index in indices_to_keep]
                tab_child_n_plies_per_angle = [
                    tab_child_n_plies_per_angle[index] \
                    for index in indices_to_keep]
                tab_child_lampam = [tab_child_lampam[index] \
                                    for index in indices_to_keep]
                # tab_child_penalty_oopo = [tab_child_penalty_oopo[index] \
                #                           for index in indices_to_keep]
                tab_child_penalty_diso = [tab_child_penalty_diso[index] \
                                          for index in indices_to_keep]
                # tab_child_penalty_contig = [tab_child_penalty_contig[index] \
                #                             for index in indices_to_keep]
                # tab_child_penalty_10 = [tab_child_penalty_10[index] \
                #                         for index in indices_to_keep]
                tab_child_obj_no_constraints = [
                    tab_child_obj_no_constraints[index] \
                    for index in indices_to_keep]
                tab_child_obj_constraints = [
                    tab_child_obj_constraints[index] \
                    for index in indices_to_keep]


            ### save local solutions as global solutions
            for ind in range(len(tab_child_obj_constraints)):

                pdls_tab.append(tab_child_pdl[ind])
                all_pdl_indices_tab.append(tab_all_pdl_indices[ind])
                last_pdl_index_tab.append(tab_last_pdl_index[ind])
                penalty_spacing_tab.append(tab_penalty_spacing[ind])
                angle_queue_tab.append(tab_angle_queue[ind])
                n_plies_per_angle_tab.append(tab_child_n_plies_per_angle[ind])
                lampam_tab.append(tab_child_lampam[ind])
                penalty_diso_tab.append(tab_child_penalty_diso[ind])
                # penalty_contig_tab.append(tab_child_penalty_contig[ind])
                # penalty_oopo_tab.append(tab_child_penalty_oopo[ind])
                # penalty_10_tab.append(tab_child_penalty_10[ind])
                obj_constraints_tab = np.hstack((
                    obj_constraints_tab, tab_child_obj_constraints[ind]))
                obj_no_constraints_tab.append(
                    tab_child_obj_no_constraints[ind])


        ### remove duplicates
        to_del = []
        for ind_pdl_1 in range(len(pdls_tab)):
            for ind_pdl_2 in range(ind_pdl_1 + 1, len(pdls_tab)):
                if is_same_pdl(pdls_tab[ind_pdl_1],
                               pdls_tab[ind_pdl_2],
                               thick_to_thin=True,
                               ind_ref=multipanel.reduced.ind_ref):
                    to_del.append(ind_pdl_1)
                    break
        to_del.sort(reverse=True)
        for ind_to_del in to_del:
            del pdls_tab[ind_to_del]
            del all_pdl_indices_tab[ind_to_del]
            del last_pdl_index_tab[ind_to_del]
            del penalty_spacing_tab[ind_to_del]
            del angle_queue_tab[ind_to_del]
            del n_plies_per_angle_tab[ind_to_del]
            del lampam_tab[ind_to_del]
            del penalty_diso_tab[ind_to_del]
            # del penalty_contig_tab[ind_to_del]
            # del penalty_oopo_tab[ind_to_del]
            # del penalty_10_tab[ind_to_del]
            del obj_no_constraints_tab[ind_to_del]
            obj_constraints_tab = np.delete(obj_constraints_tab,
                                            np.s_[ind_to_del])


        #### global pruning for ply-drop layout rules
        indices_to_keep = []
        penalty_spacing_tab_for_pruning = np.copy(penalty_spacing_tab)

        if ind_step != n_steps - 1 \
        and len(penalty_spacing_tab_for_pruning) \
        > parameters.global_node_limit2:

            while len(indices_to_keep) < parameters.global_node_limit2:

                min_value = min(penalty_spacing_tab_for_pruning)
                indices_to_add = np.where(
                penalty_spacing_tab_for_pruning == min_value)[0]
                for elem in indices_to_add:
                    indices_to_keep.append(elem)
                    penalty_spacing_tab_for_pruning[elem] = 1000

            indices_to_keep.sort()
            pdls_tab = [pdls_tab[index] for index in indices_to_keep]
            all_pdl_indices_tab = [all_pdl_indices_tab[index] \
                                   for index in indices_to_keep]
            last_pdl_index_tab = [last_pdl_index_tab[index] \
                                  for index in indices_to_keep]
            penalty_spacing_tab = [penalty_spacing_tab[index] \
                               for index in indices_to_keep]
            angle_queue_tab = [angle_queue_tab[index] \
                               for index in indices_to_keep]
            n_plies_per_angle_tab = [n_plies_per_angle_tab[index] \
                                     for index in indices_to_keep]
            lampam_tab = [lampam_tab[index] \
                          for index in indices_to_keep]
            penalty_diso_tab = [penalty_diso_tab[index] \
                                for index in indices_to_keep]
            # penalty_contig_tab = [penalty_contig_tab[index] \
            #                       for index in indices_to_keep]
            # penalty_oopo_tab = [penalty_oopo_tab[index] \
            #                     for index in indices_to_keep]
            # penalty_10_tab = [penalty_10_tab[index] \
            #                   for index in indices_to_keep]
            obj_constraints_tab = [obj_constraints_tab[index] \
                                   for index in indices_to_keep]
            obj_no_constraints_tab = [obj_no_constraints_tab[index] \
                                      for index in indices_to_keep]

#        print('len(obj_constraints_tab) before global pruning stiffness',
#              len(obj_constraints_tab))

        #### global pruning for the other guidelines and stiffness optimality
        indices_to_keep = []
        tab_child_obj_constraints_for_pruning \
        = np.copy(obj_constraints_tab)

        if ind_step != n_steps - 1 \
        and len(tab_child_obj_constraints_for_pruning) \
        > parameters.global_node_limit2:

            while len(indices_to_keep) < parameters.global_node_limit2:

                min_value = min(tab_child_obj_constraints_for_pruning)
                index_to_add = np.where(
                tab_child_obj_constraints_for_pruning == min_value)[0][0]
                indices_to_keep.append(index_to_add)
                tab_child_obj_constraints_for_pruning[index_to_add] = 1000

            indices_to_keep.sort()
            pdls_tab = [pdls_tab[index] for index in indices_to_keep]
            all_pdl_indices_tab = [all_pdl_indices_tab[index] \
                                   for index in indices_to_keep]
            last_pdl_index_tab = [last_pdl_index_tab[index] \
                                  for index in indices_to_keep]
            penalty_spacing_tab = [penalty_spacing_tab[index] \
                               for index in indices_to_keep]
            angle_queue_tab = [angle_queue_tab[index] \
                               for index in indices_to_keep]
            n_plies_per_angle_tab = [n_plies_per_angle_tab[index] \
                                     for index in indices_to_keep]
            lampam_tab = [lampam_tab[index] \
                          for index in indices_to_keep]
            penalty_diso_tab = [penalty_diso_tab[index] \
                                for index in indices_to_keep]
            # penalty_contig_tab = [penalty_contig_tab[index] \
            #                       for index in indices_to_keep]
            # penalty_oopo_tab = [penalty_oopo_tab[index] \
            #                     for index in indices_to_keep]
            # penalty_10_tab = [penalty_10_tab[index] \
            #                   for index in indices_to_keep]
            obj_constraints_tab = [obj_constraints_tab[index] \
                                   for index in indices_to_keep]
            obj_no_constraints_tab = [obj_no_constraints_tab[index] \
                                      for index in indices_to_keep]

#        print('len(obj_constraints_tab) after global pruning',
#              len(obj_constraints_tab))

    if len(obj_constraints_tab) == 0:
        return False, reduced_sst, reduced_lampam, reduced_ss

    # select best repaired solution
    index = np.argmin(obj_constraints_tab)

    for ind_panel in range(multipanel.reduced.ind_ref):
        reduced_sst_panel = pdls_tab[index][ind_panel]
        reduced_sst[ind_panel] = reduced_sst_panel
        reduced_lampam[ind_panel] = lampam_tab[index][ind_panel]
        reduced_ss[ind_panel] = reduced_sst_panel[reduced_sst_panel != -1]

    ## check for symmetry
    if constraints.sym:
        for elem in reduced_ss[0: multipanel.reduced.ind_ref]:
            for ind in range(elem.size // 2):
                if elem[ind] != elem[- ind - 1]:
                    raise Exception('reduced_ss not symmetric')
        for elem in reduced_sst[0: multipanel.reduced.ind_ref]:
            for ind in range(n_plies_max // 2):
                if elem[ind] != elem[- ind - 1]:
                    raise Exception('reduced_sst not symmetric')

    ## test for the partial lamination parameters
    reduced_lampam_test = calc_lampam(
        reduced_ss[0: multipanel.reduced.ind_ref], constraints)
    if not abs(reduced_lampam[0: multipanel.reduced.ind_ref] \
               - reduced_lampam_test).all() < 1e-13:
        raise Exception("""
beam search does not return group lamination parameters matching
the group stacking sequences.""")

    ## test for the ply counts
    for ind_panel in range(multipanel.reduced.ind_ref):
        if reduced_ss[ind_panel].size \
        != multipanel.reduced.n_plies_in_panels[ind_panel]:
            raise Exception("""
Wrong ply counts in the laminate. This should not happen.""")

    return True, reduced_sst, reduced_lampam, reduced_ss
