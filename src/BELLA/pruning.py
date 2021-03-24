# -*- coding: utf-8 -*-
"""
Pruning during guide laminate lay-up optimisation
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
#from src.divers.pretty_print import print_lampam, print_ss, print_list_ss
from src.guidelines.external_contig import external_contig
from src.guidelines.internal_contig import internal_contig2
from src.guidelines.disorientation import is_diso

def pruning_diso_contig_damtol(
        child_ss,
        mother_ss_bot,
        level,
        n_plies_to_optimise,
        constraints,
        mother_ss_top=None,
        has_middle_ply=False):
    '''
    performs the pruning for disorientation, damage tolerance and contiguity
    design guidelines during ply orientation optimisation

    INPUTS:

    - child_ss: possible fibre orientations for the new ply
    - level: level in the beam search tree
    - constraints: set of design guidelines
    - n_plies_to_optimise: number of plies to optimise during BELLA step 2
    - mother_ss_bot: beginning of the partial lay-up of the ply group being
    optimised
    - mother_ss_bot: end of the partial lay-up of the ply group being optimised
    design
    - has_middle_ply: True if one panel at least has a middle ply
    '''
    # =========================================================================
    #       pruning for middle ply symmetry
    # =========================================================================
    if constraints.sym and level == n_plies_to_optimise - 1 and has_middle_ply:
        child_ss = np.array([0, 90], int)

    # =========================================================================
    #      pruning for damage tolerance
    # =========================================================================
    my_set = set([45, -45])

    if constraints.dam_tol:

        if level == 0:
            for ind in range(child_ss.size)[:: -1]:
                if child_ss[ind] not in my_set:
                    child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            if child_ss.size == 0:
                return None
            return child_ss

        elif not constraints.sym and level == n_plies_to_optimise - 1:
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

            if level == n_plies_to_optimise - 2 and not constraints.sym:

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

    # =========================================================================
    #       pruning for disorientation
    # =========================================================================
    if constraints.diso:
        if constraints.sym or level % 2 == 0: # plies at bottom part
            #   externally with mother_ss_bot
            if mother_ss_bot.size > 0:
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_bot[-1],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        else: # asymetric laminate top part
            #   externally with mother_ss_top
            if mother_ss_top.size > 0:
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_top[0],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        if not constraints.sym and level == n_plies_to_optimise - 1:
            # last ply asymmetric laminates
            if level % 2 == 1: # top part
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_bot[-1],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)
            else: # bottom part
                for ind in range(child_ss.size)[:: -1]:
                    if not is_diso(child_ss[ind], mother_ss_top[0],
                                   constraints.delta_angle):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        if child_ss.size == 0:
            return None
    # =========================================================================
    #     pruning for the contiguity constraint
    # =========================================================================
    if constraints.contig:

        # not last ply
        if not level == n_plies_to_optimise - 1:

            # general case
            if constraints.sym or level % 2 == 0: # bottom ply

                for ind in range(child_ss.size)[:: -1]:
                    #   externally with mother_ss_bot
                    test, _ = external_contig(
                        angle=np.array((child_ss[ind],)),
                        n_plies_group=1,
                        constraints=constraints,
                        ss_before=mother_ss_bot)
                    if test.size == 0:
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)
                        continue


        # last ply
        else:

            # symmetric laminate with no middle ply
            if constraints.sym and not has_middle_ply:
                ss_before = mother_ss_bot[
                    mother_ss_bot.size - constraints.n_contig:]
                for ind in range(child_ss.size)[:: -1]:
                    new_stack = np.hstack((
                        ss_before,
                        child_ss[ind],
                        child_ss[ind],
                        np.flip(ss_before, axis=0)))
                    if not internal_contig2(new_stack, constraints):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

            # symmetric laminate with middle ply
            elif constraints.sym and has_middle_ply:
                ss_before = mother_ss_bot[
                    mother_ss_bot.size - constraints.n_contig:]
                for ind in range(child_ss.size)[:: -1]:
                    new_stack = np.hstack((
                        ss_before,
                        child_ss[ind],
                        np.flip(ss_before, axis=0)))
                    if not internal_contig2(new_stack, constraints):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

            else: # not symmetric
                ss_before = mother_ss_bot[
                    mother_ss_bot.size - constraints.n_contig:]
                ss_after = mother_ss_top[:constraints.n_contig]
                for ind in range(child_ss.size)[:: -1]:
                    if not internal_contig2(
                            new_stack=np.hstack((
                                ss_before, child_ss[ind], ss_after)),
                            constraints=constraints):
                        child_ss = np.delete(child_ss, np.s_[ind], axis=0)

        if child_ss.size == 0:
            return None

    return child_ss