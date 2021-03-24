# -*- coding: utf-8 -*-
"""
Functions to calculate moments of areas

"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np


def calc_mom_of_areas(constraints, targets, ply_order, n_plies_in_groups):
    """
    calulates ply moments of areas

    OUTPUS

    - mom_areas[ply_index, 0]: signed area of ply of index 'ply_index'
    - mom_areas[ply_index, 1]: signed first moment of area of ply of index
    'ply_index'
    - mom_areas[ply_index, 2]: signed second moment of area of ply of index
    'ply_index'

    - cummul_mom_areas[:, 0/1/2]: cummulated areas/first/second moments of
    areas of the plies in the order in which plies are optimised

    - group_mom_areas[:, 0/1/2]: cummulated areas/first/second moments of
    areas of ply groups in the order in which plies are optimised

    INPUTS

    - constraints: lay-up design guidelines
    - targets: target lamination parameters and ply counts
    - ply_order: ply indices sorted in the order in which plies are optimised
    - n_plies_in_groups: number of plies in each group of plies
    """
    group_mom_areas = np.zeros((n_plies_in_groups.size, 3), float)

    if constraints.sym:

        ply_indices = np.arange(targets.n_plies // 2 + targets.n_plies % 2)
        mom_areas = np.zeros((
            targets.n_plies // 2 + targets.n_plies % 2, 3), float)

        pos_bot = (2 / targets.n_plies) * ply_indices - 1
        pos_top = (2 / targets.n_plies) * (ply_indices + 1) - 1

        if targets.n_plies % 2:
            pos_top[-1] = 0

        mom_areas[:, 0] = pos_top - pos_bot
        mom_areas[:, 1] = pos_top**2 - pos_bot**2
        mom_areas[:, 2] = pos_top**3 - pos_bot**3

        n_plies_in_group = 0
        ind_ply_group = 0
        mom_areas_ply_group = np.zeros((3,), float)

        cummul_mom_areas = np.zeros((
            targets.n_plies // 2 + targets.n_plies % 2, 3), float)

        for ply_index in range(targets.n_plies // 2 + targets.n_plies % 2):

            cummul_mom_areas[ply_index:, :] += abs(mom_areas[ply_index, :])

            n_plies_in_group += 1
            mom_areas_ply_group += abs(mom_areas[ply_index, :])

            if n_plies_in_group == n_plies_in_groups[ind_ply_group]:
                group_mom_areas[ind_ply_group, :] = mom_areas_ply_group
                ind_ply_group += 1
                n_plies_in_group = 0
                mom_areas_ply_group = np.zeros((3,), float)

    else:
        mom_areas = np.zeros((targets.n_plies, 3), float)
        cummul_mom_areas = np.zeros((targets.n_plies, 3), float)

        ply_indices = np.arange(targets.n_plies)
        pos_bot = ((2 / targets.n_plies) * ply_indices - 1)[ply_order]
        pos_top = ((2 / targets.n_plies) * (ply_indices + 1) - 1)[ply_order]

        mom_areas[:, 0] = pos_top - pos_bot
        mom_areas[:, 1] = pos_top**2 - pos_bot**2
        mom_areas[:, 2] = pos_top**3 - pos_bot**3
        mom_areas /= 2

        n_plies_in_group = 0
        ind_ply_group = 0
        mom_areas_ply_group = np.zeros((3,), float)

        for ply_index in range(targets.n_plies - 1):

            cummul_mom_areas[ply_index:, :] += abs(mom_areas[ply_index, :])

            n_plies_in_group += 1
            mom_areas_ply_group[:] += abs(mom_areas[ply_index, :])

            if n_plies_in_group == n_plies_in_groups[ind_ply_group]:
                group_mom_areas[ind_ply_group, :] = mom_areas_ply_group
                ind_ply_group += 1
                n_plies_in_group = 0
                mom_areas_ply_group = np.zeros((3,), float)

        pos_mom_areas = np.array([
            (abs(pos_top[-1]) + abs(pos_bot[-1])) / 2,
            (abs(pos_top[-1]**2) + abs(pos_bot[-1]**2)) / 2,
            (abs(pos_top[-1]**3) + abs(pos_bot[-1]**3)) / 2])

        cummul_mom_areas[-1, :] += pos_mom_areas
        mom_areas_ply_group += pos_mom_areas
        group_mom_areas[-1, :] += mom_areas_ply_group

    return mom_areas, cummul_mom_areas, group_mom_areas


if __name__ == "__main__":
    print('*** Test for the functions calc_moment_of_areas ***\n')
    import sys
    sys.path.append(r'C:\BELLA_and_LAYLA')

    from src.LAYLA_V02.constraints import Constraints
    from src.LAYLA_V02.targets import Targets
    from src.LAYLA_V02.ply_order import calc_ply_order
    constraints = Constraints(sym=True)
    targets = Targets(n_plies=21)
    ply_order = calc_ply_order(constraints, targets)
    n_plies_in_groups = np.array([5, 6])
    mom_areas, cummul_mom_areas, group_mom_areas = calc_mom_of_areas(
        constraints, targets, ply_order, n_plies_in_groups)
    print(mom_areas)
    print(cummul_mom_areas)
    print(group_mom_areas, sum(group_mom_areas))
