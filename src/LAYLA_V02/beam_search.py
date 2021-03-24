#  - * -  coding: utf - 8  - * -
"""
Beam search for a ply group search
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import math as ma
import numpy as np

sys.path.append(r'C:\BELLA_and_LAYLA')
from src.divers.pretty_print import print_lampam, print_ss, print_list_ss
from src.guidelines.ten_percent_rule import calc_penalty_10_pc
from src.guidelines.balance import calc_penalty_bal
from src.guidelines.ipo_oopo import calc_penalty_ipo_oopo_ss
from src.guidelines.one_stack import check_lay_up_rules
from src.CLA.lampam_functions import calc_lampam_from_delta_lp_matrix
from src.LAYLA_V02.objectives import calc_obj_multi_ss, objectives
from src.LAYLA_V02.pruning import pruning_diso_contig_damtol
from src.RELAY.repair_ss import repair_ss
from src.divers.arrays import max_arrays

# do you want to save the success rate of the repair strategy?
SAVE_SUCCESS_RATE = True

def beam_search(
        levels,
        lampam_current,
        lampam_weightings,
        group_size,
        targets,
        parameters,
        constraints,
        n_plies_per_angle,
        cummul_mom_areas,
        delta_lampams,
        last_group=False,
        mat_prop=None,
        not_constraints=None,
        random_for_10=0,
        middle_ply=0,
        ss_top=None,
        ss_bot=None):
    '''
    beam search for a ply group search

    INPUTS

    - levels: ply indices
    - n_plies_in_groups: number of plies in each group of plies
    - lampam_current: lamination parameters of all plies in the laminate but
    not the plies in the group under optimisation
    - n_plies_per_angle: ply counts in each fibre direction of all plies in the
    laminate but not the plies in the group under optimisation
    - lampam_weightings: lamination parameter weightings at each search level
    - group_size: size of the ply group under determination
    - parameters: parameters of the optimiser
    - constraints: lay-up design guidelines
    - targets: target lamination parameters and ply counts
    - mat_prop: material properties
    - cummul_mom_areas: cummulated ply moments of areas
    - last_group: flag indicating of the last ply group is to be optimised
    - not_constraints: design guidelines that should not be satisfied
    - random_for_10: number to decide a ply orientation in which the 10% rule
    must not be satisfied
    - delta_lampams: ply partial lamination parameters
    - middle_ply: middle ply position of symmetric laminates, 0 otherwise
    - ss_top_ini: lay-up of plies at the top of the laminate
    - ss_bot_ini:lay-up of plies at the bottom of the laminate
    '''
    results = beamsearchResults()

    # details for when the 10% rule should be violated
    if not_constraints is not None and not_constraints.rule_10_percent:
        n_plies_0_lim = ma.ceil(not_constraints.percent_0 * targets.n_plies)
        n_plies_90_lim = ma.ceil(
            not_constraints.percent_90 * targets.n_plies)
        n_plies_45_lim = ma.ceil(
            not_constraints.percent_45 * targets.n_plies)
        n_plies_135_lim = ma.ceil(
            not_constraints.percent_135 * targets.n_plies)
        n_plies_45_135_lim = ma.ceil(
            not_constraints.percent_45_135 * targets.n_plies)
        if constraints.sym:
            n_plies_0_lim = ma.ceil(n_plies_0_lim / 2)
            n_plies_90_lim = ma.ceil(n_plies_90_lim / 2)
            n_plies_45_lim = ma.ceil(n_plies_45_lim / 2)
            n_plies_135_lim = ma.ceil(n_plies_135_lim / 2)
            n_plies_45_135_lim = ma.ceil(n_plies_45_135_lim / 2)

    # simplifies ss_top and ss_bot for faster feasibility checks?
    ss_bot_simp = np.copy(ss_bot)
    if not constraints.sym:
        ss_top_simp = np.copy(ss_top)

    if constraints.contig or constraints.diso:
        if ss_bot.size > constraints.n_contig:
            ss_bot_simp = ss_bot[-constraints.n_contig:]

        if not constraints.sym:
            if ss_top.size > constraints.n_contig:
                ss_top_simp = ss_top[:constraints.n_contig]

    # nodal stacking sequences
    ss_bot_tab = [np.array((), dtype='int16')]
    if not constraints.sym:
        ss_top_tab = [np.array((), dtype='int16')]
    # nodal lamination parameters
    lampam_tab = lampam_current.reshape((1, 12))
    # estimate functions
    obj_const_tab = np.zeros((1,), float)
    # nodal ply counts in each fibre orientations
    ply_counts_tab = n_plies_per_angle.reshape((1, constraints.n_set_of_angles))

    if last_group:
        ss_final = np.array((), dtype='int16').reshape((0, targets.n_plies))

    for local_level in range(group_size):

        level = levels[local_level]
#        print('********************************')
#        print('level', level)

        for el in range(obj_const_tab.size):

            try:
                mother_lampam = lampam_tab[0]
            except IndexError:
                raise Exception(
                    'Infeasible beam-search, increase the branching limits')

            mother_ss_bot = ss_bot_tab.pop(0)
            if not constraints.sym:
                mother_ss_top = ss_top_tab.pop(0)
            mother_ply_counts = ply_counts_tab[0]


            # deletion of the mother lamination parameters/ply counts/estimate
            # function values from the queues
            lampam_tab = np.delete(lampam_tab, np.s_[0], axis=0)
            ply_counts_tab = np.delete(ply_counts_tab, np.s_[0], axis=0)
            obj_const_tab = np.delete(obj_const_tab, np.s_[0])

            # branching
            child_ss = np.copy(constraints.set_of_angles)

            # pruning
            if constraints.sym:
                child_ss = pruning_diso_contig_damtol(
                    child_ss=child_ss,
                    mother_ss_bot=mother_ss_bot,
                    ss_bot_simp=ss_bot_simp,
                    level=level,
                    constraints=constraints,
                    targets=targets)
            else:
                child_ss = pruning_diso_contig_damtol(
                    child_ss=child_ss,
                    mother_ss_bot=mother_ss_bot,
                    mother_ss_top=mother_ss_top,
                    ss_bot_simp=ss_bot_simp,
                    ss_top_simp=ss_top_simp,
                    level=level,
                    constraints=constraints,
                    targets=targets)

#            print(child_ss)

            if child_ss is None:
                continue

#            print('mother_ss', mother_ss)
#            print('child_ss after pruning diso/contig/damtol', child_ss.T)

            # caluclate the lamination parameters
            child_lampam = np.matlib.repmat(mother_lampam, child_ss.size, 1)
            for myindex in range(child_ss.size):
                child_lampam[myindex] += delta_lampams[
                    level, constraints.ind_angles_dict[child_ss[myindex]]]

            # calculate the ply counts
            # & compute the stacking sequences
            # & calculate the penalties for the 10% rule
            penalty_10 = np.array((), float)
            n_solution = 0
            for indd in range(child_ss.size)[::-1]:

                ply_counts = np.copy(mother_ply_counts)
                index = constraints.ind_angles_dict[child_ss[indd]]
                if middle_ply != 0 and local_level == group_size - 1:
                    ply_counts[index] += 1/2
                else:
                    ply_counts[index] += 1

                # pruning for not_constraints
                if not_constraints is not None \
                and not_constraints.rule_10_percent:
                    if random_for_10 == 0 \
                    and ply_counts[constraints.index0] + 1 >= n_plies_0_lim:
                        continue
                    if random_for_10 == 1\
                    and ply_counts[constraints.index90] + 1 >= n_plies_90_lim:
                        continue
                    if random_for_10 == 2 \
                    and ply_counts[constraints.index45] + 1 >= n_plies_45_lim:
                        continue
                    if random_for_10 == 3 \
                    and ply_counts[constraints.index135] + 1 >= n_plies_135_lim:
                        continue
                    if random_for_10 == 4 \
                    and ply_counts[constraints.index45] + \
                    ply_counts[constraints.index135] + 1 >= n_plies_45_135_lim:
                        continue

                # penalty for the 10% rule
                if local_level == group_size - 1 \
                and constraints.rule_10_percent \
                and parameters.penalty_10_pc_switch:
                    penalty_10 = np.hstack((
                        penalty_10,
                        calc_penalty_10_pc(ply_counts, constraints)))

                n_solution += 1

                if constraints.sym:
                    new_stack_bot = np.copy(mother_ss_bot)
                    new_stack_bot = np.hstack((new_stack_bot, child_ss[indd]))
                    ss_bot_tab.append(new_stack_bot)
                else:
                    new_stack_bot = np.copy(mother_ss_bot)
                    new_stack_top = np.copy(mother_ss_top)

                    if local_level % 2:
                        new_stack_top = np.hstack((
                            child_ss[indd], new_stack_top))
                    else:
                        new_stack_bot = np.hstack((
                            new_stack_bot, child_ss[indd]))
                    ss_bot_tab.append(new_stack_bot)
                    ss_top_tab.append(new_stack_top)

#                    print('new_stack_bot', new_stack_bot)
#                    print('new_stack_top', new_stack_top)

                lampam_tab = np.vstack((lampam_tab, child_lampam[indd]))
                ply_counts_tab = np.vstack((ply_counts_tab, ply_counts))

            if penalty_10.size == 0:
                penalty_10 = 0

            results.n_nodes += n_solution
            if n_solution == 0:
                continue # go to next branching

            # estimate function values with no constraints
#            print('no - constraints')
#            for el in lampam_tab[lampam_tab.shape[0] - n_solution:]:
#                print('lampam', el[0:4])
            obj_no_constraints = objectives(
                lampam=lampam_tab[lampam_tab.shape[0] - n_solution:],
                targets=targets,
                lampam_weightings=lampam_weightings[level],
                constraints=constraints,
                parameters=parameters,
                mat_prop=mat_prop)

            if last_group and local_level == group_size - 1:

                # full stacking sequences
                for indd in range(n_solution)[::-1]:
                    indddd = len(ss_bot_tab) - indd - 1
#                    print('child_ss_tab[indddd]', child_ss_tab[indddd])

                    if constraints.sym:
                        ss = np.hstack((
                            ss_bot,
                            ss_bot_tab[indddd],
                            np.flip(ss_bot_tab[indddd], axis=0),
                            np.flip(ss_bot, axis=0))).astype('int16')
                        if middle_ply != 0:
                            ss = np.delete(ss, np.s_[middle_ply], axis=0)
                    else:
                        ss = np.hstack((
                            ss_bot,
                            ss_bot_tab[indddd],
                            ss_top_tab[indddd],
                            ss_top)).astype('int16')
                    if ss.size != targets.n_plies:
                        print('ss.size', ss.size)
                        print('targets.n_plies', targets.n_plies)
                        raise Exception("This should not happen")

                    # repair
                    results.n_designs_last_level += 1
#                    print('before repair')
#                    print_ss(ss)
#                    print('obj_no_constraints ', obj_no_constraints )
                    ss, success_repair, n_obj_func_D_calls = repair_ss(
                        ss, constraints, parameters, targets.lampam, True)
#                    results.n_obj_func_calls += n_obj_func_D_calls
#                    print('repair successful?', success_repair)
#                    print('after repair')
#                    print_ss(ss)
#                    print('obj_no_constraints ', obj_no_constraints )
                    if success_repair:
                        results.n_designs_repaired += 1
#                        print('no - constraints - after success full repair')
                        obj_no_constraints[n_solution - indd - 1] = objectives(
                            lampam=calc_lampam_from_delta_lp_matrix(
                                ss, constraints, delta_lampams),
                            targets=targets,
                            lampam_weightings=lampam_weightings[level],
                            constraints=constraints,
                            parameters=parameters,
                            mat_prop=mat_prop)
#                        results.n_obj_func_calls += 1
                        check_lay_up_rules(ss, constraints)
                        lampam_tab[indddd] = calc_lampam_from_delta_lp_matrix(
                            ss, constraints, delta_lampams)
                    else:
#                        print('unsuccess full repair')
                        obj_no_constraints[n_solution - indd - 1] = 1e10
                    ss_final = np.vstack((ss_final, ss))


            # calculation the penalties for the in-plane and out-of-plane
            # orthotropy requirements based on lamination parameters
            penalty_ipo_lampam, penalty_oopo = calc_penalty_ipo_oopo_ss(
                lampam_tab[lampam_tab.shape[0] - n_solution:],
                constraints=constraints,
                parameters=parameters,
                cummul_areas=cummul_mom_areas[level, 0],
                cummul_sec_mom_areas=cummul_mom_areas[level, 2])
#            print('penalty_ipo_lampam', penalty_ipo_lampam.shape)
#            print('penalty_oopo', penalty_oopo.shape)

            # calculation the penalties for the in-plane orthotropy
            # requirements based on ply counts
            penalty_ipo_pc = np.zeros((n_solution,))
            if constraints.ipo and parameters.penalty_bal_switch:
                penalty_ipo_pc = calc_penalty_bal(
                    ply_counts,
                    constraints=constraints,
                    cummul_areas=cummul_mom_areas[level, 0])
#                print('penalty_ipo_pc', penalty_ipo_pc.shape)

            penalty_bal_ipo = max_arrays(penalty_ipo_pc, penalty_ipo_lampam)

            # calculation of the bounds
            obj_const = calc_obj_multi_ss(
                objective=obj_no_constraints,
                penalty_10=penalty_10,
                penalty_bal_ipo=penalty_bal_ipo,
                penalty_oopo=penalty_oopo,
                coeff_10=parameters.coeff_10,
                coeff_bal_ipo=parameters.coeff_bal_ipo,
                coeff_oopo=parameters.coeff_oopo)
            obj_const_tab = np.hstack((obj_const_tab, obj_const))

#            print('')
#            print('lampam_tab', lampam_tab.shape)
#            print('obj_const_tab', obj_const_tab.shape)
#            print('mother_ss', mother_ss)
#            print('mother_lampam', mother_lampam[0:4])
#            print('child_ss', child_ss.T)
#            print('lampam_tab', lampam_tab[-1][0:4])
#            print('lampam_tab', lampam_tab[-2][0:4])
#            print('lampam_tab', lampam_tab[-3][0:4])
#            print('obj_const', obj_const, '\n\n')

            # local pruning
            n_local_nodes = obj_const.size
            if last_group and local_level != group_size - 1:
                node_limit = parameters.local_node_limit_final
            else:
                node_limit = parameters.local_node_limit
            n_excess_nodes = n_local_nodes - node_limit
            if local_level != group_size - 1 and n_excess_nodes > 0:
                obj_const_tab_to_del = np.copy(obj_const)
                to_del = []
                for counter in range(n_excess_nodes):
                    ind_max = np.argmax(obj_const_tab_to_del)
                    obj_const_tab_to_del[ind_max] = -6666
                    to_del.append(ind_max + obj_const_tab.size - n_local_nodes)
                ss_bot_tab = [elem for ind, elem in enumerate(ss_bot_tab) \
                              if ind not in to_del]
                if not constraints.sym:
                    ss_top_tab = [elem for ind, elem in enumerate(ss_top_tab) \
                                  if ind not in to_del]
                lampam_tab = np.delete(lampam_tab, np.s_[to_del], axis=0)
                ply_counts_tab = np.delete(
                    ply_counts_tab, np.s_[to_del], axis=0)
                obj_const_tab = np.delete(obj_const_tab, np.s_[to_del])

        if obj_const_tab.size == 0:
            raise Exception(
                'Infeasible beam-search, increase the branching limits')

        # global pruning
        if last_group and local_level == group_size - 2:
            if obj_const_tab.size > parameters.global_node_limit_p:
                obj_const_to_del = np.copy(obj_const_tab)
                for counter in range(parameters.global_node_limit_p):
                    my_min = np.argmin(obj_const_to_del)
                    obj_const_to_del[my_min] = 6666
                to_keep = obj_const_to_del == 6666
                to_del = np.invert(to_keep).astype(int)
                to_del = [i for i, x in enumerate(to_del) if x]
                ss_bot_tab = [elem for ind, elem in enumerate(ss_bot_tab) \
                              if ind not in to_del]
                if not constraints.sym:
                    ss_top_tab = [elem for ind, elem in enumerate(ss_top_tab) \
                                  if ind not in to_del]
                lampam_tab = np.delete(lampam_tab, np.s_[to_del], axis=0)
                ply_counts_tab = np.delete(
                    ply_counts_tab, np.s_[to_del], axis=0)
                obj_const_tab = np.delete(obj_const_tab, np.s_[to_del])

        elif local_level != group_size - 1:
            if obj_const_tab.size > parameters.global_node_limit:
                obj_const_to_del = np.copy(obj_const_tab)
                for counter in range(parameters.global_node_limit):
                    my_min = np.argmin(obj_const_to_del)
                    obj_const_to_del[my_min] = 6666
                to_keep = obj_const_to_del == 6666
                to_del = np.invert(to_keep).astype(int)
                to_del = [i for i, x in enumerate(to_del) if x]
                ss_bot_tab = [elem for ind, elem in enumerate(ss_bot_tab) \
                              if ind not in to_del]
                if not constraints.sym:
                    ss_top_tab = [elem for ind, elem in enumerate(ss_top_tab) \
                                  if ind not in to_del]
                lampam_tab = np.delete(lampam_tab, np.s_[to_del], axis=0)
                ply_counts_tab = np.delete(
                    ply_counts_tab, np.s_[to_del], axis=0)
                obj_const_tab = np.delete(obj_const_tab, np.s_[to_del])

    ## return results
    if last_group:
        ss_bot_tab = np.copy(ss_final)
        if SAVE_SUCCESS_RATE:
            ss_final = np.array([ss_final[ind] \
                                 for ind in range(len(obj_const_tab)) \
                                 if obj_const_tab[ind] < 1e10])
            if ss_final.size == 0:
                results.n_designs_repaired_unique = 0
            else:
                results.n_designs_repaired_unique = np.unique(
                    ss_final, axis=0).shape[0]

        # identification of the best leaf
        ind_baby = np.argmin(obj_const_tab)
        if obj_const_tab[ind_baby] > 0.99*1e10:
            raise Exception("""
No successfull repair during beam search, increase the branching limits""")
        results.ss_best = ss_bot_tab[ind_baby]
        results.lampam_best = lampam_tab[ind_baby]
        results.ply_counts = ply_counts_tab[ind_baby]

    else:
        # identification of the best leaf
        ind_baby = np.argmin(obj_const_tab)
        if obj_const_tab[ind_baby] > 0.99*1e10:
            raise Exception("""
No successfull repair during beam search, increase the branching limits""")
        results.ss_bot_best = ss_bot_tab[ind_baby]
        if not constraints.sym:
            results.ss_top_best = ss_top_tab[ind_baby]
        results.lampam_best = lampam_tab[ind_baby]
        results.ply_counts = ply_counts_tab[ind_baby]

#    t_end = time.time()
#    print('time beam search', t_end - t_beg)

    return results


class beamsearchResults():
    " An object for storing the results of a ply group search in LAYLA"
    def __init__(self):
        "Initialise the results of a ply group search in LAYLA"
        # solution stacking sequence
        self.ss_bot_best = None
        self.ss_top_best = None
        self.ss_best = None
        # solution lamination parameters
        self.lampam_best = None
        # solution lamination parameters at each outer step
        self.ply_counts = None
        # number of nodes reached in the search tree
        self.n_nodes = 0
#        # number of objective function calls
#        self.n_obj_func_calls = 0
        # number of nodes reached at the last level of the search tree
        self.n_designs_last_level = 0
        # number of repaired nodes reached at the last level of the search tree
        self.n_designs_repaired = 0

    def __repr__(self):
        " Display object "

        return f'''
Results with LAYLA:

    Stacking sequence: {self.ss_best}
    Lamination parameters 1-4: {self.lampam_best[:4]}
    Lamination parameters 5-8: {self.lampam_best[4:8]}
    Lamination parameters 9-12: {self.lampam_best[8:]}
'''
