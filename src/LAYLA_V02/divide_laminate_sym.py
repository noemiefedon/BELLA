# -*- coding: utf-8 -*-
"""
Function partitioning laminates into groups of plies for symmetric laminates
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import math as ma
import numpy as np

class PartitioningError(Exception):
    "Exceptions during partitioning laminates into groups of plies"

def divide_laminate_sym(parameters, targets, step=0):
    '''
    performs the partitioning of a symmetric laminate into groups of plies.
    The number of groups is minised with a preference for large groups of
    plies, especially for the last ply groups corresponding to the inside of
    the laminate.

    OUTPUTS

    - n_plies_in_groups: number of plies in each group of plies
    - pos_first_ply_groups: position of the first ply of each group
    with a numbering starting from the bottom to the top of the laminate
    - n_groups: the number of steps to be performed by the algorithm

    INPUTS

    - parameters: parameters of the optimiser
    - targets.n_plies: number of plies
    - step: number of the outer steps
    '''

    # middle_ply = 0 if there is no ply overlapping the mid-surface, otherwise
    # middle_ply is equal to the position number of this ply
    if targets.n_plies % 2 == 1:
        middle_ply = int((targets.n_plies + 1) / 2)
    else:
        middle_ply = 0

    mini = ma.ceil(ma.floor(targets.n_plies/2)/parameters.group_size_max[step])
    maxi = ma.floor(ma.floor(targets.n_plies/2)/parameters.group_size_min)
    if mini > maxi:
        raise PartitioningError('''
Partitioning of the laminate not possible with for symmetric laminates with
the current ply-group size limitations.
Try increasing the maximum number of plies per group or reducing the minimum
number of plies per group.
''')

    # iteration with increasing number of groups
    for n_groups in np.arange(mini, maxi + 1):

        # ?
        missing = n_groups * parameters.group_size_max[step] \
        - ma.floor(targets.n_plies/2)
        if missing > (parameters.group_size_max[step] \
                      - parameters.group_size_min)*n_groups:
            continue

        #
        if n_groups == 0:
            continue

        # distribution of parameters.group_size_min plies in each group
        n_plies_in_groups = parameters.group_size_min * \
        np.ones((n_groups,), int)

        # n_extra: number of remaining plies to be distributed in the groups
        n_extra = ma.floor(targets.n_plies/2) \
        - n_groups*parameters.group_size_min

        # n_full_groups: number of groups that can be totally filled by the
        # distribution of the remianing plies
        if n_extra >= parameters.group_size_max[step] \
        - parameters.group_size_min and n_extra != 0:
            n_full_groups = n_extra \
            // (parameters.group_size_max[step]-parameters.group_size_min)
            n_extra = n_extra \
            % (parameters.group_size_max[step]-parameters.group_size_min)
        else:
            n_full_groups = 0

        # filling of the n_full_groups
        n_plies_in_groups[n_groups - n_full_groups:] \
        = parameters.group_size_max[step]

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
    if middle_ply == 0 and sum(n_plies_in_groups)*2 != targets.n_plies:
        raise PartitioningError('Wrong partitioning!')
    if middle_ply != 0 and sum(n_plies_in_groups)*2 + 1 != targets.n_plies:
        raise PartitioningError('Wrong partitioning!')

    if middle_ply != 0:
        n_plies_in_groups[-1] += 1

#    if n_groups > 1:
#        if parameters.group_size_max[step] < 5:
#            print('''
#The number of plies of the last group (parameters.group_size_max) is
#recommended to be equal to or greater than 5.
#''')
#        if parameters.group_size_min < 4:
#            print('''
#The number of plies of the smaller groups (parameters.group_size_min) is
#recommended to be equal to or greater than 4.
#''')
#    elif n_groups == 1:
#        if parameters.group_size_max[step] < 4:
#            print('''
#The number of plies of the last group (parameters.group_size_max) is
#recommended to be equal to or greater than 4.
#''')
#        if parameters.group_size_min < 4:
#            print('''
#The number of plies of the smaller groups (parameters.group_size_min) is
#recommended to be equal to or greater than 4.
#''')
    if n_groups > maxi:
        raise PartitioningError('''
No partition possible of the plies into groups of smaller size
parameters.group_size_minand bigger size parameters.group_size_max.
Increase parameters.group_size_max or decrease parameters.group_size_min.
''')
    return n_plies_in_groups, pos_first_ply_groups, n_groups

if __name__ == "__main__":
    import sys
    sys.path.append(r'C:\BELLA_and_LAYLA')
    from src.LAYLA_V02.parameters import Parameters
    from src.LAYLA_V02.constraints import Constraints
    from src.LAYLA_V02.targets import Targets
    constraints = Constraints(sym=True)
    targets = Targets(n_plies=21)
    parameters = Parameters(constraints=constraints,
                            n_outer_step=2,
                            group_size_min=3,
                            group_size_max=5)
    n_plies_in_groups, pos_first_ply_groups, n_groups \
    = divide_laminate_sym(parameters, targets)
    print(n_plies_in_groups)
    print(pos_first_ply_groups)
    print(n_groups)
