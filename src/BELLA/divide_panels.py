# -*- coding: utf-8 -*-
"""
Functions to divide the panels into groups of plies

- divide_panels
    divides the plies of each panel into groups of plies and calculates the
    lamination parameter coefficients for the objective functions that accounts
    for the normalised positive areas and first and second moments of areas

- forbidden_ply_counts
    determines the ply counts that lead to an  impossible partitioning of the
    plies into groups of allowed number of plies
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import math as ma
import numpy as np

sys.path.append(r'C:\BELLA')
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.panels import Panel
from src.BELLA.multipanels import MultiPanel

class PartitioningError(Exception):
    " Errors occuring during the partition of the panels into groups of plies"

def divide_panels(multipanel, parameters, constraints):
    """
    divides the plies of each panel into groups of plies and calculates the
    lamination parameter coefficients for the objective functions that accounts
    for the normalised positive areas and first and second moments of areas

    Guidelines:
    1:  The first two outer plies should not be stopped
    2:  The number of ply drops should be minimal (not butt joints)
    3:  The ply drops should be distributed as evenly as possible along the
        thickness of the laminates
    4:  If this is not exactly possible the ply drops should rather be
        concentrated in the larger groups (because smaller groups have a
        smaller design space)
    5:  Then ply drops away from the middle plane are prefered to limit fibre
        waviness

    INPUTS

    - multipanel: multipanel structure
    - parameters: optimiser parameters
    - constraints: set of design and manufacturing constraints
    """
    #===================
    # Division of the thickest panel into groups of plies
    #===================
    if constraints.sym:

        # evaluate the number of ply orientations for covering plies
        n_lam = multipanel.n_plies_max - 2*constraints.n_covering

        mini = ma.ceil(ma.floor(n_lam/2)/parameters.group_size_max)
        maxi = ma.floor(ma.floor(n_lam/2)/parameters.group_size_min)
        if mini > maxi:
            raise PartitioningError("""
Partitioning of the laminate not possible with the current group size
limitations !
    Try increasing the maximum number of plies per group or reducing the
    minimum number of plies per group.""")

        # The loop ensures that the division into groups is conformed to the
        # constraints on the group sizes, and groups of maximum size are
        # prefered
        for n_groups in np.arange(mini, maxi + 1):

            # ?
            missing = n_groups*parameters.group_size_max - ma.floor(n_lam/2)
            if missing > (parameters.group_size_max \
                          - parameters.group_size_min)*n_groups:
                continue

            #
            if n_groups == 0:
                continue

                # distribution of parameters.group_size_min plies in each group
            n_plies_in_groups = parameters.group_size_min * \
            np.ones((n_groups,), int)

            # n_extra: number of remaining plies to be distributed in the groups
            n_extra = ma.floor(n_lam/2) \
            - n_groups*parameters.group_size_min

            # n_full_groups: number of groups that can be totally filled by the
            # distribution of the remianing plies
            if n_extra >= parameters.group_size_max \
            - parameters.group_size_min and n_extra != 0:
                n_full_groups = n_extra \
                // (parameters.group_size_max-parameters.group_size_min)
                n_extra = n_extra \
                % (parameters.group_size_max-parameters.group_size_min)
            else:
                n_full_groups = 0

            # filling of the n_full_groups
            n_plies_in_groups[n_groups - n_full_groups:] \
            = parameters.group_size_max

            # distribution of the last remaining plies
            if n_extra != 0:
                n_plies_in_groups[n_groups - n_full_groups - 1] += n_extra

            # pos_first_ply_groups: position of the first ply of each group in
            # the order in which they appear in the stacking sequence
            pos_first_ply_groups = np.zeros((n_groups,), int)
    #        pos_first_ply_groups[0] = 1
            for ind in np.arange(1, n_groups):
                pos_first_ply_groups[ind] \
                = pos_first_ply_groups[ind - 1] + n_plies_in_groups[ind - 1]
            break

        # checking group sizes are correct (should not return an error!!!)
        if multipanel.thick_panel_has_middle_ply:
            if sum(n_plies_in_groups)*2 + 1 != n_lam:
                raise PartitioningError('Wrong partitioning!')
        elif sum(n_plies_in_groups)*2 != n_lam:
            raise PartitioningError('Wrong partitioning!')

#        if middle_ply != 0:
#            n_plies_in_groups[-1] += 1

        if n_groups > maxi:
            raise PartitioningError('''
No partition possible of the plies into groups of smaller size
parameters.group_size_min and bigger size parameters.group_size_max.
Increase parameters.group_size_max or decrease parameters.group_size_min.
''')

    else: # for non symmetric laminates

        # Evaluate the number of ply orientations for covering plies
        n_lam = multipanel.n_plies_max - 2*constraints.n_covering

        mini = ma.ceil(n_lam/parameters.group_size_max)
        maxi = ma.floor(n_lam/parameters.group_size_min)

        if mini > maxi:
            raise PartitioningError("""
Partitioning of the laminate not possible with the current group size
limitations !
    Try increasing the maximum number of plies per group or reducing the
    minimum number of plies per group.""")

        # iteration with increasing number of groups
        for n_groups in np.arange(mini, maxi + 1):
            # ?
            missing = n_groups * parameters.group_size_max \
            - n_lam
            if missing > (parameters.group_size_max \
                          - parameters.group_size_min)*n_groups:
                continue

            #
            if n_groups == 0:
                continue

            # distribution of parameters.group_size_min plies in each group
            n_plies_in_groups = parameters.group_size_min \
            * np.ones((n_groups,), int)

            # n_extra: number of remaining plies to be distributed in groups
            n_extra = n_lam - n_groups*parameters.group_size_min

            # n_full_groups: number of groups that can be totally filled by the
            # distribution of the remaining plies
            if n_extra >= parameters.group_size_max \
            - parameters.group_size_min and n_extra != 0:
                n_full_groups = n_extra // (
                    parameters.group_size_max - parameters.group_size_min)
                n_extra %= (
                    parameters.group_size_max - parameters.group_size_min)
            else:
                n_full_groups = 0

            if n_full_groups > 0:
                n_plies_in_groups[-n_full_groups:] \
                = parameters.group_size_max
            # Addition of the last other plies
            if n_extra != 0:
                n_plies_in_groups[-(n_full_groups + 1)] += n_extra

            # order_of_groups: group sizes in the order in which they
            # appear in the stacking sequence
            middle_point = ma.ceil(n_groups/2)
            order_of_groups = np.zeros((n_groups,), int)
            order_of_groups[:middle_point] \
            = n_plies_in_groups[0:2*middle_point:2]
            order_of_groups[middle_point:] = np.flip(
                n_plies_in_groups[1:n_groups:2], axis=0)

            # pos_of_groups: position of the first ply of each
            # group in the order they appear in the final stacking sequence
            pos_of_groups = np.zeros((n_groups,), int)
    #        pos_of_groups[0] = 1
            for ind in np.arange(1, n_groups):
                pos_of_groups[ind] = pos_of_groups[ind - 1] \
                + order_of_groups[ind - 1]

            pos_first_ply_groups = np.ones((n_groups,), int)
            pos_first_ply_groups[0:2*middle_point:2] \
            = pos_of_groups[:middle_point]
            pos_first_ply_groups[1:n_groups:2] = np.flip(
                pos_of_groups[middle_point:], axis=0)
            break

        # checking group sizes are correct (should not return an error!!!)
        if sum(n_plies_in_groups) != n_lam:
            raise PartitioningError('Wrong partitioning!')

        if n_groups > maxi:
            raise PartitioningError('''
No partition possible of the plies into groups of smaller size
parameters.group_size_minand bigger size parameters.group_size_max.
Increase parameters.group_size_max or decrease parameters.group_size_min.
''')

#        # correction for when middle ply not in largest laminate
#        if multipanel.has_middle_ply \
#        and multipanel.reduced.panels[-1].middle_ply == 0:
#            n_plies_in_groups += 1

    # number of plies per group for thickest laminates
    multipanel.reduced.n_plies_per_group = n_plies_in_groups
    # position of the group first plies for thickest laminates
    multipanel.reduced.n_first_plies = pos_first_ply_groups
    # number of groups
    multipanel.reduced.n_groups = n_groups
    # percent_ of the laminate thickness associated to each group
    # for thickest laminatess
    percent_thickness = multipanel.reduced.n_plies_per_group \
    / sum(multipanel.reduced.n_plies_per_group)
#    print('percent_thickness\n', percent_thickness, '\n')

#    print('multipanel.reduced.n_plies_per_group',
#          multipanel.reduced.n_plies_per_group)
#    print('multipanel.reduced.n_first_plies',
#          multipanel.reduced.n_first_plies)
#    print('multipanel.n_plies_max', multipanel.n_plies_max)

    #===================
    # Division of other panels into groups of plies
    #===================
    for ind_panel, panel in enumerate(multipanel.reduced.panels):
        # ----- if panel is one of the thickest panels ----- #
        if panel.n_plies == multipanel.n_plies_max:
            #--------------------------
            # number of plies per group
            #--------------------------
            panel.n_plies_per_group \
            = multipanel.reduced.n_plies_per_group.astype(int)
            #--------------------------
            # position of the group first plies
            #--------------------------
            panel.n_first_plies = multipanel.reduced.n_first_plies
            #--------------------------
            # Number of ply drops number of ply drops for each group compared
            # to the groups thickest of the thickest laminate
            #--------------------------
            panel.n_ply_drops = np.zeros((multipanel.reduced.n_groups,))
        else: # ----- if panel is not one of the thickest panels ----- #
            #--------------------------
            # Distribution of the ply drops in the groups of the panel
            # from comparison with the thickest panel
            #--------------------------
            n_plies_panel = panel.n_plies
            if constraints.sym: # middle plies not included!
                if n_plies_panel % 2:
                    n_plies_panel -= 1
                n_drops = (multipanel.n_plies_max - n_plies_panel)/2
            else:
                n_drops = multipanel.n_plies_max - n_plies_panel

            #print('n_plies', n_plies_panel, '\n')

            #print('n_drops', n_drops, '\n')
            #print('percent_thickness', percent_thickness)

            # attempt of a uniform distribution of ply drops
            n_ply_drops = np.floor(n_drops*percent_thickness).astype(float)

#            # correction for potential middle ply
#            if panel.middle_ply != 0:
#                if n_ply_drops[-1] > 0:
#                    n_ply_drops[-1] -= 0.5
#                elif np.allclose(np.zeros((multipanel.reduced.n_groups,)),
#                                 n_ply_drops):
#                    n_ply_drops[-1] = 0.5
#                else:
#                    for index_group in range(
#                            multipanel.reduced.n_groups -1)[::-1]:
#                        if n_ply_drops[index_group] != 0:
#                            n_ply_drops[index_group] -= 1
#                            n_ply_drops[-1] = 0.5
#                            break

            missing = int(n_drops - np.sum(n_ply_drops))
            #print('Ply drops per group:', n_ply_drops, '\n')
            #print('missing ply drops:', missing, '\n')
            # Addition of the missing ply drops
            iteration = 0
            while missing > 0:
                for index_group in range(multipanel.reduced.n_groups):
                    if multipanel.reduced.n_plies_per_group[index_group] == \
                    parameters.group_size_max - iteration:
                        n_ply_drops[index_group] += 1
                        missing -= 1
                        if missing == 0:
                            break
                iteration += 1
            # Here, panel.n_ply_drops = chosen combination of ply drops
            panel.n_ply_drops = n_ply_drops
            #--------------------------
            # number of plies per group
            #--------------------------
            panel.n_plies_per_group = (multipanel.reduced.n_plies_per_group \
                                   - panel.n_ply_drops).astype(int)
            #--------------------------
            # position of each group first ply
            #--------------------------
            if constraints.sym:
                n_first_plies = np.zeros(
                    (multipanel.reduced.n_groups,), dtype=int)
                n_first_plies[0] = constraints.n_covering

                for jjj in range(1, multipanel.reduced.n_groups):
                    n_first_plies[jjj] = n_first_plies[jjj - 1] \
                                        + panel.n_plies_per_group[jjj - 1]
                panel.n_first_plies = n_first_plies
                #print(panel.n_plies_per_group)
                #print(n_first_plies)
            else:
                n_first_plies_b = np.zeros(
                    (multipanel.reduced.n_groups,), dtype=int)
                n_first_plies = np.zeros(
                    (multipanel.reduced.n_groups,), dtype=int)
                n_first_plies[0] = constraints.n_covering
                n_first_plies_b[0] = constraints.n_covering
                # order_of_groups: group sizes in the order
                # they appear in the final stacking sequence - in a column
                middle_point = ma.ceil(multipanel.reduced.n_groups/2)
                order_of_groups = np.zeros((multipanel.reduced.n_groups,), int)
                order_of_groups[:middle_point] = \
                panel.n_plies_per_group[0:2*middle_point:2]
                order_of_groups[middle_point: ] = np.flip(
                    panel.n_plies_per_group[1:multipanel.reduced.n_groups:2],
                    axis=0)
                for index_group in np.arange(1, multipanel.reduced.n_groups):
                    n_first_plies_b[index_group] \
                   = n_first_plies_b[index_group - 1] \
                    + order_of_groups[index_group-1]
                # Filling the start positions for each group in the table
                n_first_plies[0:2*middle_point:2] = \
                n_first_plies_b[:middle_point]
                n_first_plies[1:multipanel.reduced.n_groups:2] \
                = np.flip(n_first_plies_b[middle_point: ], axis=0)
                panel.n_first_plies = n_first_plies
                #print(panel.n_first_plies)


def forbidden_ply_counts(constraints, parameters):
    """
    determines the ply counts that lead to an  impossible partitioning of the
    plies into groups of allowed number of plies

    INPUTS

    - constraints.n_plies_min: minimum number of plies for a laminate
    - constraints.n_plies_max: maximum number of plies for a laminate
    - group_size_min: minimum number of plies for a group
    - group_size_max: maximum number of plies for a group
    - constraints.sym = True for symmetric laminates
    """
    result = []
    for n_plies_test in range(constraints.n_plies_min,
                              constraints.n_plies_max + 1):
        try:
            #===================
            # Division of the thickest panel into groups of plies
            #===================
            if constraints.sym:
                n_lam = n_plies_test - 2*constraints.n_covering
                mini = ma.ceil(ma.floor(n_lam/2)/parameters.group_size_max)
                maxi = ma.floor(ma.floor(n_lam/2)/parameters.group_size_min)
                if mini > maxi:
                    raise PartitioningError("""
Partitioning of the laminate not possible with the current group size
limitations !
    Try increasing the maximum number of plies per group or reducing the
    minimum number of plies per group.""")
            else: # for non symmetric laminates
                n_lam = n_plies_test - 2*constraints.n_covering
                mini = ma.ceil(n_lam/parameters.group_size_max)
                maxi = ma.floor(n_lam/parameters.group_size_min)
                if mini > maxi:
                    raise PartitioningError("""
Partitioning of the laminate not possible with the current group size
limitations !
    Try increasing the maximum number of plies per group or reducing the
    minimum number of plies per group.""")
        except PartitioningError:
            result.append(n_plies_test)
        else:
            pass
    return np.array(result)


if __name__ == "__main__":
    print('*** Test for the function divide_panels ***\n')
    constraints = Constraints(n_covering=0, 
                              covering=False, 
                              sym=True)
    parameters = Parameters(constraints=constraints,
                            group_size_max=10, 
                            group_size_min=6)
    n_plies_target1 = 81
    n_plies_target2 = 70
    n_plies_target3 = n_plies_target2
    n_plies_target4 = n_plies_target2
    panel_1 = Panel(ID=1,
                    n_plies=n_plies_target1,
                    constraints=constraints,
                    neighbour_panels=[2])
    panel_2 = Panel(ID=2,
                    n_plies=n_plies_target2,
                    constraints=constraints,
                    neighbour_panels=[1])
    panel_3 = Panel(ID=3,
                    n_plies=n_plies_target3,
                    constraints=constraints,
                    neighbour_panels=[2])
    panel_4 = Panel(ID=4,
                    n_plies=n_plies_target4,
                    constraints=constraints,
                    neighbour_panels=[2])
    multipanel = MultiPanel(panels=[panel_1, panel_2, panel_3, panel_4])
    multipanel.from_mp_to_blending_strip(constraints)
    
    print(f'Panel 1: {multipanel.reduced.panels[0].n_plies} plies')
    print(f'Panel 2: {multipanel.reduced.panels[1].n_plies} plies')
#    print(f'Panel 3: {multipanel.reduced.panels[2].n_plies} plies')
    divide_panels(multipanel, parameters, constraints)
    print('\n')
    print(f'Panel 1: number of plies per groups = {multipanel.reduced.panels[0].n_plies_per_group}',
          sum(multipanel.reduced.panels[0].n_plies_per_group))
    print(f'Panel 2: number of plies per groups = {multipanel.reduced.panels[1].n_plies_per_group}',
          sum(multipanel.reduced.panels[1].n_plies_per_group))
#    print(f'Panel 3: number of plies per groups = {multipanel.reduced.panels[2].n_plies_per_group}',
#          sum(multipanel.reduced.panels[2].n_plies_per_group))
    print('\n')
    print(f'Panel 1: number of ply drops per group = {multipanel.reduced.panels[0].n_ply_drops}')
    print(f'Panel 2: number of ply drops per group = {multipanel.reduced.panels[1].n_ply_drops}')
#    print(f'Panel 3: number of ply drops per group = {multipanel.reduced.panels[2].n_ply_drops}')
    print('\n')
    print(f'Panel 1: position group first plies = {multipanel.reduced.panels[0].n_first_plies}')
    print(f'Panel 2: position group first plies = {multipanel.reduced.panels[1].n_first_plies}')
#    print(f'Panel 3: position group first plies = {multipanel.reduced.panels[2].n_first_plies}')

    print('\n*** Test for the function forbidden_ply_counts ***\n')
    constraints = Constraints(
        n_covering=1, 
        covering=True, 
        sym=True, 
        n_plies_min=10, 
        n_plies_max=100)
    parameters = Parameters(constraints, 
                            group_size_min=5, 
                            group_size_max=7)
    result = forbidden_ply_counts(constraints, parameters)
    print(f'Forbidden ply counts: {result}\n')
