# -*- coding: utf-8 -*-
"""
Function partitioning laminates into groups of plies for  asymmetric laminates
"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import math as ma
import numpy as np

class PartitioningError(Exception):
    "Exceptions during partitioning laminates into groups of plies"

def divide_laminate_asym(parameters, targets, step=0):
    '''
    performs the partitioning of the plies for U-laminates.

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
    mini = ma.ceil(targets.n_plies/parameters.group_size_max[step])
    maxi = ma.floor(targets.n_plies/parameters.group_size_min)
    if mini > maxi:
        raise PartitioningError('''
Partitioning of the laminate not possible with for asymmetric laminates with
the current ply-group size limitations.
Try increasing the maximum number of plies per group or reducing the minimum
number of plies per group.
''')

    # iteration with increasing number of groups
    for n_groups in np.arange(mini, maxi + 1):
        # ?
        missing = n_groups * parameters.group_size_max[step] \
        - targets.n_plies
        if missing > (parameters.group_size_max[step] \
                      - parameters.group_size_min)*n_groups:
            continue

        #
        if n_groups == 0:
            continue

        # distribution of parameters.group_size_min plies in each group
        n_plies_in_groups = parameters.group_size_min \
        * np.ones((n_groups,), int)

        # n_extra: number of remaining plies to be distributed in the groups
        n_extra = targets.n_plies - n_groups*parameters.group_size_min

        # n_n_full_groups: number of groups that can be totally filled by the
        # distribution of the remianing plies
        if n_extra >= parameters.group_size_max[step] \
        - parameters.group_size_min and n_extra != 0:
            n_full_groups = n_extra // (
                parameters.group_size_max[step] - parameters.group_size_min)
            n_extra %= (
                parameters.group_size_max[step] - parameters.group_size_min)
        else:
            n_full_groups = 0

        if n_full_groups > 0:
            n_plies_in_groups[-n_full_groups:] \
            = parameters.group_size_max[step]
        # Addition of the last other plies
        if n_extra != 0:
            n_plies_in_groups[-(n_full_groups + 1)] += n_extra

        # order_of_groups: group sizes in the order in which they
        # appear in the stacking sequence
        middle_point = ma.ceil(n_groups/2)
        order_of_groups = np.zeros((n_groups,), int)
        order_of_groups[:middle_point] = n_plies_in_groups[0:2*middle_point:2]
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
        pos_first_ply_groups[0:2*middle_point:2] = pos_of_groups[:middle_point]
        pos_first_ply_groups[1:n_groups:2] = np.flip(
            pos_of_groups[middle_point:], axis=0)
        break

    # checking group sizes are correct (should not return an error!!!)
    if sum(n_plies_in_groups) != targets.n_plies:
        raise PartitioningError('Wrong partitioning!')

#    if n_groups == 1:
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
#    elif n_groups > 1:
#        if parameters.group_size_min < 4:
#            print('''
#The number of plies of the smaller groups (parameters.group_size_min) is
#recommended to be equal to or greater than 4.
#''')
#        if parameters.group_size_max[step] < 5:
#            print('''
#The number of plies of the last group (parameters.group_size_max) is
#recommended to be equal to or greater than 5.
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
    from src.LAYLA_V02.constraints import Constraints
    from src.LAYLA_V02.parameters import Parameters
    from src.LAYLA_V02.targets import Targets
    constraints = Constraints(sym=False, bal=True)
    targets = Targets(n_plies=200)
    parameters = Parameters(
        constraints=constraints,
        n_outer_step=5,
        group_size_min=3,
        group_size_max=6)
    n_plies_in_groups, pos_first_ply_groups, n_groups \
    = divide_laminate_asym(parameters, targets)
    print(n_plies_in_groups)
    print(pos_first_ply_groups)
    print(n_groups)
