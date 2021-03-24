# -*- coding: utf-8 -*-
"""
Functions to calculate the order in which plies are optimised

"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np

def calc_ply_order(multipanel, constraints):
    """
    calulates the order in which plies are optimised

    OUTPUTS

    - ply_order[ind_panel]: array of the ply indices sorted in the order in
    which plies are optimised (middle ply of symmetric laminates included)

    INPUTS

    - constraints: lay-up design guidelines
    - multipanel: multi-panel structure
    """
    ply_order = []

    for panel in multipanel.reduced.panels:
        if constraints.sym:
            ply_order.append(
                np.arange(panel.n_plies // 2 + panel.n_plies % 2))
        else:
            order_before_sorting = np.arange(panel.n_plies)
            ply_order_new = np.zeros((panel.n_plies,), int)
            ply_order_new[0::2] = order_before_sorting[
                :panel.n_plies // 2 + panel.n_plies % 2]
            ply_order_new[1::2] = order_before_sorting[
                panel.n_plies // 2 + panel.n_plies % 2:][::-1]
            ply_order.append(ply_order_new)

    return ply_order
