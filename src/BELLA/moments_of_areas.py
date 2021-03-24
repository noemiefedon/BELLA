# -*- coding: utf-8 -*-
"""
Functions to calculate moments of areas

"""
__version__ = '2.0'
__author__ = 'Noemie Fedon'

import numpy as np


def calc_mom_of_areas(multipanel, constraints, ply_order):
    """
    calulates ply moments of areas

    OUTPUS

    - mom_areas_plus[panel_index, ply_index, 0]:
        signed area of ply of index 'ply_index' in panel of index 'panel_index'
    - mom_areas_plus[panel_index, ply_index, 1]:
        signed first moment of area of ply of index 'ply_index' in panel of
        index 'panel_index'
    - mom_areas_plus[panel_index, ply_index, 2]:
        signed second moment of area of ply of index 'ply_index' in panel of
        index 'panel_index'

    - mom_areas[panel_index, ply_index, 0]:
        area of ply of index 'ply_index' in panel of index 'panel_index'
    - mom_areas[panel_index, ply_index, 1]:
        first moment of area of ply of index 'ply_index' in panel of index
        'panel_index'
    - mom_areas[panel_index, ply_index, 2]:
        second moment of area of ply of index 'ply_index' in panel of index
        'panel_index'

    INPUTS

    - constraints: lay-up design guidelines
    - multipanel: multi-panel structure
    - ply_order: ply indices sorted in the order in which plies are optimised
    """
    mom_areas_plus = []
    mom_areas = []

    for ind_panel, panel in enumerate(multipanel.reduced.panels):

        if constraints.sym:

            ply_indices = np.arange(panel.n_plies // 2 + panel.n_plies % 2)
            mom_areas_panel = np.zeros((
                panel.n_plies // 2 + panel.n_plies % 2, 3), float)
            mom_areas_plus_panel = np.zeros((
                panel.n_plies // 2 + panel.n_plies % 2, 3), float)

            pos_bot = (2 / panel.n_plies) * ply_indices - 1
            pos_top = (2 / panel.n_plies) * (ply_indices + 1) - 1

            if panel.n_plies % 2:
                pos_top[-1] = 0

#            print(pos_bot)
#            print(pos_top)

            mom_areas_panel[:, 0] = pos_top - pos_bot
            mom_areas_panel[:, 1] = 0
            mom_areas_panel[:, 2] = pos_top**3 - pos_bot**3

            mom_areas_plus_panel[:, 0] = pos_top - pos_bot
            mom_areas_plus_panel[:, 1] = abs(pos_top**2 - pos_bot**2)
            mom_areas_plus_panel[:, 2] = pos_top**3 - pos_bot**3

        else:
            mom_areas_panel = np.zeros((panel.n_plies, 3), float)
            mom_areas_plus_panel = np.zeros((panel.n_plies, 3), float)

            ply_indices = np.arange(panel.n_plies)

            print(ply_indices, type(ply_indices))
            print(ply_order, type(ply_order))
            pos_bot = ((2 / panel.n_plies) \
                       * ply_indices - 1)[ply_order[ind_panel]]
            pos_top = ((2 / panel.n_plies) \
                       * (ply_indices + 1) - 1)[ply_order[ind_panel]]

            mom_areas_panel[:, 0] = pos_top - pos_bot
            mom_areas_panel[:, 1] = pos_top**2 - pos_bot**2
            mom_areas_panel[:, 2] = pos_top**3 - pos_bot**3
            mom_areas_panel /= 2


            for ind in range(panel.n_plies):

                if pos_top[ind] * pos_bot[ind] >= 0:
                    mom_areas_plus_panel[
                        ind, 0] = abs(pos_top[ind] - pos_bot[ind])
                    mom_areas_plus_panel[
                        ind, 1] = abs(pos_top[ind]**2 - pos_bot[ind]**2)
                    mom_areas_plus_panel[
                        ind, 2] = abs(pos_top[ind]**3 - pos_bot[ind]**3)
                else:
                    mom_areas_plus_panel[
                        ind, 0] = abs(pos_top[ind]) + abs(pos_bot[ind])
                    mom_areas_plus_panel[
                        ind, 1] = abs(pos_top[ind]**2) + abs(pos_bot[ind]**2)
                    mom_areas_plus_panel[
                        ind, 2] = abs(pos_top[ind]**3) + abs(pos_bot[ind]**3)

            mom_areas_plus_panel /= 2

        mom_areas_plus.append(mom_areas_plus_panel)
        mom_areas.append(mom_areas_panel)

    return mom_areas_plus, mom_areas


def calc_mom_of_areas2(multipanel, constraints, mom_areas_plus, pdl,
                       n_plies_to_optimise):
    """
    calulates ply moments of areas

    OUTPUS

    - cummul_areas[panel_index][ply_index, 0]:


    INPUTS

    - constraints: lay-up design guidelines
    - multipanel: multi-panel structure
    - pdl: ply drop layout
    - mom_areas_plus[panel_index, ply_index, 0]:
        signed area of ply of index 'ply_index' in panel of index 'panel_index'
    - mom_areas_plus[panel_index, ply_index, 1]:
        signed first moment of area of ply of index 'ply_index' in panel of
        index 'panel_index'
    - mom_areas_plus[panel_index, ply_index, 2]:
        signed second moment of area of ply of index 'ply_index' in panel of
        index 'panel_index'
    - n_plies_to_optimise: number of plies to optimise during BELLA step 2
    """

    cummul_areas = np.zeros(
        (multipanel.reduced.n_panels, n_plies_to_optimise), float)
    cummul_first_mom_areas = np.zeros(
        (multipanel.reduced.n_panels, n_plies_to_optimise), float)
    cummul_sec_mom_areas = np.zeros(
        (multipanel.reduced.n_panels, n_plies_to_optimise), float)

    for ind_panel, panel in enumerate(multipanel.reduced.panels):
        counter_plies = -1
        for index_ply in range(n_plies_to_optimise):
            if pdl[ind_panel, index_ply] != -1:
                counter_plies += 1

                cummul_areas[ind_panel, index_ply:] \
                += mom_areas_plus[ind_panel][counter_plies][0]

                cummul_first_mom_areas[ind_panel, index_ply:] \
                += mom_areas_plus[ind_panel][counter_plies][1]

                cummul_sec_mom_areas[ind_panel, index_ply:] \
                += mom_areas_plus[ind_panel][counter_plies][2]

    return cummul_areas, cummul_first_mom_areas, cummul_sec_mom_areas


if __name__ == "__main__":
    print('*** Test for the functions calc_moment_of_areas ***\n')
    import sys
    sys.path.append(r'C:\BELLA')

    from src.BELLA.constraints import Constraints
    from src.BELLA.panels import Panel
    from src.BELLA.multipanels import MultiPanel
    from src.BELLA.parameters import Parameters
    from src.BELLA.obj_function import ObjFunction
    from src.BELLA.ply_order import calc_ply_order
    from src.BELLA.pdl_ini import create_initial_pdls
    from src.BELLA.divide_panels import  divide_panels

    constraints = Constraints(sym=False)
#    constraints = Constraints(sym=True)
    obj_func_param = ObjFunction(constraints)

    parameters = Parameters(constraints)
    panel1 = Panel(1, constraints, neighbour_panels=[], n_plies=6)
    multipanel = MultiPanel([panel1])

    parameters = Parameters(constraints)
    panel1 = Panel(1, constraints, neighbour_panels=[1], n_plies=10)
    panel2 = Panel(2, constraints, neighbour_panels=[1], n_plies=8)
    multipanel = MultiPanel([panel1, panel2])
    multipanel.from_mp_to_blending_strip(
        constraints, parameters.n_plies_ref_panel)
    ply_order = calc_ply_order(multipanel, constraints)

    indices = ply_order[-1]
    n_plies_to_optimise = indices.size
    mom_areas_plus, mom_areas = calc_mom_of_areas(
        multipanel, constraints, ply_order)

    print('mom_areas_plus')
    print(mom_areas_plus[0])
    print(mom_areas_plus[1])
    print(sum(mom_areas_plus[0]))
    print(sum(mom_areas_plus[1]))
    print('mom_areas')
    print(mom_areas[0])
    print(mom_areas[1])
    print(sum(mom_areas[0]))
    print(sum(mom_areas[1]))

    print('*** Test for the functions calc_mom_of_areas2 ***\n')
    divide_panels(multipanel, parameters, constraints)
    pdl = create_initial_pdls(
        multipanel, constraints, parameters, obj_func_param)[0]

    cummul_areas, cummul_first_mom_areas, cummul_sec_mom_areas = \
    calc_mom_of_areas2(
        multipanel, constraints, mom_areas_plus, pdl, n_plies_to_optimise)

    print('cummul_areas')
    print(cummul_areas)

    print('cummul_first_mom_areas')
    print(cummul_first_mom_areas)

    print('cummul_sec_mom_areas')
    print(cummul_sec_mom_areas)

