# -*- coding: utf-8 -*-
"""
All possible ply lamination parameters
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
sys.path.append(r'C:\BELLA')
from src.CLA.lampam_functions import calc_delta_lampam_mp_3
from src.CLA.lampam_functions import calc_delta_lampam


def calc_delta_lampams(multipanel, constraints, mom_areas, ply_order):
    """
    calulates all ply partial lamination parameters in a multipanel structure

    INPUTS

    - multipanel: multi-panel structure
    - constraints: set of constraints
    - mom_areas[panel_index, ply_index, 0]:
        area of ply of index 'ply_index' in panel of index 'panel_index'
    - mom_areas[panel_index, ply_index, 1]:
        first moment of area of ply of index 'ply_index' in panel of index
        'panel_index'
    - mom_areas[panel_index, ply_index, 2]:
        second moment of area of ply of index 'ply_index' in panel of index
        'panel_index'
    - ply_order: ply indices sorted in the order in which plies are optimised
    """
    delta_lampams = []


    for ind_panel, panel in enumerate(multipanel.reduced.panels):

        if constraints.sym:
            n_plies_panel = panel.n_plies // 2 + panel.n_plies % 2
            delta_lampams_panel = np.empty((
                n_plies_panel,
                constraints.n_set_of_angles,
                12), float)
            for ind_ply in range(n_plies_panel):


                delta_lampams_panel[ind_ply, :, 0:4] \
                = mom_areas[ind_panel][ind_ply, 0] * constraints.cos_sin
                delta_lampams_panel[ind_ply, :, 4:8] = 0
                delta_lampams_panel[ind_ply, :, 8:12] \
                = mom_areas[ind_panel][ind_ply, 2] * constraints.cos_sin
            if panel.n_plies % 2 and ind_ply == panel.middle_ply_index:
                delta_lampams_panel[ind_ply, :, :] /= 2

        else:
            delta_lampams_panel = np.empty((
                panel.n_plies, constraints.n_set_of_angles,
                12), float)
            for ind_ply in range(delta_lampams_panel.shape[0]):
                delta_lampams_panel[ind_ply, :, 0:4] \
                = mom_areas[ind_panel][ind_ply, 0] * constraints.cos_sin
                delta_lampams_panel[ind_ply, :, 4:8] \
                = mom_areas[ind_panel][ind_ply, 1] * constraints.cos_sin
                delta_lampams_panel[ind_ply, :, 8:12] \
                = mom_areas[ind_panel][ind_ply, 2] * constraints.cos_sin

        delta_lampams.append(delta_lampams_panel)

    return delta_lampams


def calc_delta_lampams2(
        multipanel, constraints, delta_lampams, pdl, n_plies_to_optimise):
    """
    calulates all ply partial lamination parameters in a multipanel structure
    that correspond to a specific ply drop layout

    INPUTS

    - multipanel: multipanel structure
    - delta_lampams: all ply partial lamination parameters in the multipanel
    structure
    - pdl: ply drop layout
    - constraints: set of constraints
    - n_plies_to_optimise: number of plies to optimise during BELLA step 2
    """
    lampam_matrix = np.zeros((
            multipanel.reduced.n_panels,
            constraints.n_set_of_angles,
            n_plies_to_optimise,
            12), float)

    for ind_panel, panel in enumerate(multipanel.reduced.panels):

        for ind_angle in range(constraints.n_set_of_angles):
            counter_plies = -1
            for index_ply in range(n_plies_to_optimise):
                if pdl[ind_panel, index_ply] != -1:
                    counter_plies += 1
#                    print('ind_panel, ind_angle, index_ply',
#                          ind_panel, ind_angle, index_ply)
#                    print('counter_plies', counter_plies)
                    lampam_matrix[ind_panel, ind_angle, index_ply, :] \
                    = delta_lampams[ind_panel][counter_plies, ind_angle]
    return lampam_matrix


if __name__ == "__main__":
    print('*** Test for the functions calc_delta_lampams ***\n')
    import sys
    sys.path.append(r'C:\BELLA')

    from src.BELLA.constraints import Constraints
    from src.BELLA.panels import Panel
    from src.BELLA.multipanels import MultiPanel
    from src.BELLA.parameters import Parameters
    from src.BELLA.obj_function import ObjFunction
    from src.BELLA.ply_order import calc_ply_order
    from src.BELLA.moments_of_areas import calc_mom_of_areas
    from src.BELLA.pdl_ini import create_initial_pdls
    from src.BELLA.divide_panels import  divide_panels

    constraints = Constraints(sym=False)
    constraints = Constraints(sym=True)
    obj_func_param = ObjFunction(constraints)

    parameters = Parameters(constraints)
    panel1 = Panel(1, constraints, neighbour_panels=[], n_plies=6)
    multipanel = MultiPanel([panel1])

    parameters = Parameters(constraints)
    panel1 = Panel(1, constraints, neighbour_panels=[1], n_plies=16)
    panel2 = Panel(2, constraints, neighbour_panels=[1], n_plies=18)
    multipanel = MultiPanel([panel1, panel2])

    ply_order = calc_ply_order(multipanel, constraints)
    indices = ply_order[-1]
    n_plies_to_optimise = indices.size
    mom_areas_plus, mom_areas = calc_mom_of_areas(
        multipanel, constraints, ply_order)

    delta_lampams = calc_delta_lampams(
        multipanel, constraints, mom_areas, ply_order)

    print(delta_lampams[0][0][0])
    print(delta_lampams[0][-1][0])


    print('*** Test for the functions calc_delta_lampams2 ***\n')
    divide_panels(multipanel, parameters, constraints)
    pdl = create_initial_pdls(
        multipanel, constraints, parameters, obj_func_param)[0]
    lampam_matrix = calc_delta_lampams2(
        multipanel, constraints, delta_lampams, pdl, n_plies_to_optimise)