# -*- coding: utf-8 -*-
"""
Functions to calculate the order in which plies are optimised

"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def calc_ply_order(constraints, targets):
    """
    calulates the order in which plies are optimised

    OUTPUTS

    - ply_order: array of the ply indices sorted in the order in which plies
    are optimised (middle ply of symmetric laminates included)

    INPUTS

    - constraints: lay-up design guidelines
    - targets: target lamination parameters and ply counts
    """
    if constraints.sym:
        ply_order = np.arange(targets.n_plies // 2 + targets.n_plies % 2)
        return ply_order

    order_before_sorting = np.arange(targets.n_plies)
    ply_order = np.zeros((targets.n_plies,), int)
    ply_order[0::2] = order_before_sorting[
        :targets.n_plies // 2 + targets.n_plies % 2]
    ply_order[1::2] = order_before_sorting[
        targets.n_plies // 2 + targets.n_plies % 2:][::-1]
    return ply_order

def  calc_levels(ply_order, n_plies_in_groups, n_groups):
    """
    calulates the indices of the plies for each ply group optimisation

    INPUTS

    - ply_order: array of the ply indices sorted in the order in which plies
    are optimised (middle ply of symmetric laminates included)
    - n_plies_in_groups: number of plies in each ply group
    - n_groups: number of ply groups
    """
    levels_in_groups = [None]*n_groups
    for ind_group in range(n_groups):
        levels_in_groups[ind_group] = []

    ind_all_plies = 0
    for ind_group in range(n_groups):
        for ind_plies in range(n_plies_in_groups[ind_group]):
            levels_in_groups[ind_group].append(ply_order[ind_all_plies])
            ind_all_plies += 1

    return levels_in_groups


if __name__ == "__main__":
    import sys
    sys.path.append(r'C:\BELLA_and_LAYLA')
    from src.LAYLA_V02.parameters import Parameters
    from src.LAYLA_V02.constraints import Constraints
    from src.LAYLA_V02.targets import Targets
    constraints = Constraints(sym=False)
    targets = Targets(n_plies=6)
    ply_order = calc_ply_order(constraints, targets)
    print(ply_order)
