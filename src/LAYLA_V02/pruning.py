# -*- coding: utf-8 -*-
"""
Pruning during the beam search
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA_and_LAYLA')
#from src.divers.pretty_print import print_lampam, print_ss, print_list_ss
from src.guidelines.external_contig import external_contig
from src.guidelines.internal_contig import internal_contig2
from src.guidelines.disorientation import is_diso

def pruning_diso_contig_damtol(
        child_ss,
        mother_ss_bot,
        ss_bot_simp,
        level,
        constraints,
        targets,
        mother_ss_top=None,
        ss_top_simp=None):
    '''
    performs the pruning for disorientation, damage tolerance, contiguity and
    middle ply symmetry design guidelines during beam search

    INPUTS

    - level: level of the ply determination in the complete lay-up search tree
    - child_ss: possible fibre angles for new level
    - ss_bot_simp: partial lay-up before the ply group being optimised
    - ss_top_simp: partial lay-up after the ply group being optimised
    - mother_ss_bot: beginning of the partial lay-up of the ply group being
    optimised
    - mother_ss_bot: end of the partial lay-up of the ply group being optimised
    - constraints: lay-up design guidelines
    - targets.n_plies: laminate target ply counts
    '''
# =============================================================================
#       pruning for middle ply symmetry
# =============================================================================
    if constraints.sym and targets.n_plies % 2 \
    and level == targets.n_plies // 2:
        child_ss = np.array([0, 90], int)
# =============================================================================
#      pruning for damage tolerance
# =============================================================================
    my_set = set([45, -45])

    if constraints.dam_tol:

        if level == 0:
            for ind in range(child_ss.size)[:: -1]:
                if child_ss[ind] not in my_set:
                    child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            if child_ss.size == 0:
                return None
            return child_ss

        elif not constraints.sym and level == targets.n_plies - 1:
            for ind in range(child_ss.size)[:: -1]:
                if child_ss[ind] not in my_set:
                    child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            if child_ss.size == 0:
                return None
            return child_ss

        if constraints.dam_tol_rule in [2, 3]:

            if level == 1:

                if constraints.dam_tol_rule == 2:
                    for ind in range(child_ss.size)[:: -1]:
                        if child_ss[ind] != - mother_ss_bot[0]:
                            child_ss = np.delete(child_ss, np.s_[ind], axis=0)

                elif constraints.dam_tol_rule == 3:
                    for ind in range(child_ss.size)[:: -1]:
                        if child_ss[ind] not in my_set:
                            child_ss = np.delete(child_ss, np.s_[ind], axis=0)
                            continue
                        # diso
                        if constraints.diso \
                        and not is_diso(-45, 45, constraints.delta_angle):
                            if child_ss[ind] != - mother_ss_bot[0]:
                                child_ss = np.delete(
                                    child_ss, np.s_[ind], axis=0)

                if child_ss.size == 0:
                    return None
                return child_ss

            if level == targets.n_plies - 2 and not constraints.sym:

                if constraints.dam_tol_rule == 2:
                    for ind in range(child_ss.size)[:: -1]:
                        if child_ss[ind] != - mother_ss_top[-1]:
                            child_ss = np.delete(child_ss, np.s_[ind], axis=0)

                elif constraints.dam_tol_rule == 3:
                    for ind in range(child_ss.size)[:: -1]:
                        if child_ss[ind] not in my_set:
                            child_ss = np.delete(child_ss, np.s_[ind], axis=0)
                            continue
                        # diso
                        if constraints.diso \
                        and not is_diso(-45, 45, constraints.delta_angle):
                            if child_ss[ind] != - mother_ss_top[-1]:
                                child_ss = np.delete(
                                    child_ss, np.s_[ind], axis=0)

                if child_ss.size == 0:
                    return None
                return child_ss

# =============================================================================
#       pruning for disorientation
# =============================================================================
    if constraints.diso:
        if constraints.sym or level % 2: # plies at the laminate bottom part
            #   externally with ss_bot_simp
            if ss_bot_simp.size > 0 and mother_ss_bot.size == 0:
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], ss_bot_simp[-1],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            #   internally
            elif mother_ss_bot.size > 0:
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_bot[-1],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        else: # asymetric laminate top part
            #   externally with ss_top_simp
            if ss_top_simp.size > 0 and mother_ss_top.size == 0:
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], ss_top_simp[0],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            #   internally
            if mother_ss_top.size > 0:
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_top[0],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        if not constraints.sym and level == targets.n_plies - 1:
            # last ply asymmetric laminates
            if not level % 2: # check compatibility with laminate bottom part
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_bot[-1],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            else: # check compatibility with laminate top part
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_top[0],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        if child_ss.size == 0:
            return None
# =============================================================================
#     #  pruning for the contiguity constraint
# =============================================================================
    if constraints.contig:

        # laminate bottom part but not last ply
        if (constraints.sym \
            and not level == targets.n_plies // 2 + targets.n_plies % 2 - 1)\
        or (not constraints.sym and not level == targets.n_plies - 1):

            for ind in range(child_ss.size)[:: -1]:
                #   externally with ss_bot_simp
                test, _ = external_contig(
                    angle=np.array((child_ss[ind],)),
                    n_plies_group=1,
                    constraints=constraints,
                    ss_before=np.hstack((ss_bot_simp, mother_ss_bot)))
                if test.size == 0:
                    child_ss = np.delete(child_ss, np.s_[ind], axis=0)
                    continue

        # plies at the laminate bottom part but not last ply
        elif not constraints.sym and not level % 2 \
        and not level == targets.n_plies - 1:

            for ind in range(child_ss.size)[:: -1]:

                #   externally with ss_top_simp
                test, _ = external_contig(
                    angle=np.array((child_ss[ind],)),
                    n_plies_group=1,
                    constraints=constraints,
                    ss_before=np.flip(
                        np.hstack((mother_ss_top, ss_top_simp)), axis=0))
                if test.size == 0:
                    child_ss = np.delete(child_ss, np.s_[ind], axis=0)
                    continue

        # last ply of symmetric laminate
        elif constraints.sym \
        and level == targets.n_plies // 2 + targets.n_plies % 2 - 1:

            ss_before = mother_ss_bot[
                mother_ss_bot.size - constraints.n_contig:]
            if ss_before.size < constraints.n_contig:
                ss_before = np.hstack((
                    ss_bot_simp[ss_bot_simp.size \
                                - constraints.n_contig + ss_before.size:],
                    ss_before))

            if targets.n_plies % 2 == 0: # no middle ply

                for ind in range(child_ss.size)[:: -1]:
                    new_stack = np.hstack((
                        ss_before,
                        child_ss[ind],
                        child_ss[ind],
                        np.flip(ss_before, axis=0)))
                    if not internal_contig2(new_stack, constraints):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

            else: # a middle ply
                for ind in range(child_ss.size)[:: -1]:
                    new_stack = np.hstack((
                        ss_before,
                        child_ss[ind],
                        np.flip(ss_before, axis=0)))
                    if not internal_contig2(new_stack, constraints):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        # plies at the laminate bottom part and not last ply
        elif not constraints.sym and level == targets.n_plies - 1:

            ss_before = mother_ss_bot[
                mother_ss_bot.size - constraints.n_contig:]
            if ss_before.size < constraints.n_contig:
                ss_before = np.hstack((
                    ss_bot_simp[ss_bot_simp.size \
                                - constraints.n_contig + ss_before.size:],
                    ss_before))

            ss_after = mother_ss_top[:constraints.n_contig]
            if ss_after.size < constraints.n_contig:
                ss_after = np.hstack((
                    ss_after,
                    ss_top_simp[:constraints.n_contig - ss_after.size:]))

            for ind in range(child_ss.size)[:: -1]:
                if not internal_contig2(
                        new_stack=np.hstack((
                            ss_before, child_ss[ind], ss_after)),
                        constraints=constraints):
                    child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        if child_ss.size == 0:
            return None

    return child_ss
