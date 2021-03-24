# -*- coding: utf-8 -*-
"""
- calc_penalty_spacing_1ss
    calculates the ply drop spacing penalty of a panel regarding one of its
    boundary

- calc_penalty_spacing
    calculates the penalties of a ply drop layout for the ply-drop spacing rule
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
from src.BELLA.format_pdl import reduce_for_guide_based_blending


def calc_penalty_spacing(
        pdl,
        multipanel,
        constraints,
        obj_func_param=None,
        on_blending_strip=False,
        pdl_before=None,
        pdl_after=None):
    """
    calculates the penalties of a ply drop layout for the ply-drop spacing rule

    INPUTS

    - multipanel: multi-panel structure
    - pdl: matrix of ply drop layouts of the current group
    - pdl_before: matrix of ply drop layout of the previous group
    - pdl_after: matrix of ply drop layout of the group placed afterwards
    - constraints: lay-up design guidelines
    - obj_func_param: objective function parameters
    - on_blending_strip: to be calculated on blending strip
     """
    if not constraints.pdl_spacing:
        return 0

    ### blending strip --------------------------------------------------------
    if on_blending_strip:

        if pdl_before is None:
            pdl_before = np.zeros((multipanel.reduced.n_panels, 0))
        if pdl_after is None:
            pdl_after = np.zeros((multipanel.reduced.n_panels, 0))

        penalty_spacing = 0

        for ind_panel1, ind_panel2 in multipanel.reduced.boundaries:

            if pdl[ind_panel1] is None or pdl[ind_panel2] is None \
            or pdl[ind_panel1].size == 1 or pdl[ind_panel2].size == 1:
                continue

            # stack the ply drop layouts before/current/after
            layout1 = np.hstack((
                pdl_before[ind_panel1],
                pdl[ind_panel1],
                pdl_after[ind_panel1]))
            layout2 = np.hstack((
                pdl_before[ind_panel2],
                pdl[ind_panel2],
                pdl_after[ind_panel2]))

            # delete plies that does not cover the panels
            to_keep = [layout1[ind] != -1 or layout1[ind] != layout2[ind]
                       for ind in range(layout1.size)]
            layout1 = layout1[to_keep]
            if -1 not in layout1:
                layout1 = layout2[to_keep]

    #        print('ind_panel1, ind_panel2', ind_panel1, ind_panel2)
    #        print(layout1)

            penalty_spacing += (multipanel.reduced.boundary_weights[
                (ind_panel1, ind_panel2)] * calc_penalty_spacing_1ss(
                layout1, constraints.min_drop))

        return penalty_spacing

    ### multi-panel -----------------------------------------------------------
    if not hasattr(multipanel, 'reduced'):
        multipanel.from_mp_to_blending_strip(constraints, n_plies_ref_panel=1)

    reduced_pdl = reduce_for_guide_based_blending(multipanel, pdl)

#    print('reduced_pdl.shape', reduced_pdl.shape)
    return calc_penalty_spacing(
        pdl=reduced_pdl,
        multipanel=multipanel,
        constraints=constraints,
        obj_func_param=obj_func_param,
        on_blending_strip=True,
        pdl_before=pdl_before,
        pdl_after=pdl_after)


def is_same_pdl(pdl1, pdl2, ind_ref, thick_to_thin=True):
    """
    returns True if pdl1 == pdl2
    """
    if thick_to_thin:
        for index in range(ind_ref):
            if pdl1[index] is None:
                pass
            elif (pdl1[index] != pdl2[index]).any():
                return False
        return True

    for index in range(len(pdl1)):
        if pdl1[index] is None:
            pass
        elif (pdl1[index] != pdl2[index]).any():
            return False
    return True


def calc_penalty_spacing_1ss(ss, min_drop):
    """
    returns the penalty for the ply drop spacing rule by considering the
    stacking sequence in one panel in reference to one of its neighboor panel

    Sum of the missing continuous plies between ply drops divided by the length
    of the ply drop layout

    - min_drop: minimum number of continuous plies required between two block
    of dropped plies
    """
    penal = 0
    # indentify index of the successive -1 elements in the array
    index1 = 0
    while index1 < ss.size:
        if ss[index1] == -1:
            break
        else:
            index1 += 1
    index2 = index1 + 1
    while index2 < ss.size:
        while index2 < ss.size:
            if ss[index2] == -1:
                break
            else:
                index2 += 1
        # test for the ply drop spacing condition
        if index2 < ss.size and index2 - index1 - 1 < min_drop:
            penal += min_drop - (index2 - index1 - 1)
        index1 = index2
        index2 += 1
    return penal / ss.size


if __name__ == "__main__":


    print('\n*** Test for the function calc_penalty_spacing_1ss ***')
    ss = np.array([-1, 0, 0, -1, 0., 1., -1., 2., -1., 4., -1, 1., -1])
    min_drop = 2
    p = calc_penalty_spacing_1ss(ss, min_drop)
    print(ss)
    print(f'Penalty for the ply drop spacing rule spacing {p}')

    print(calc_penalty_spacing_1ss(np.array([
            1., 2., 3., 4., 5., 6.,-1., 7., 8., 9., 10., -1.,
            -1.,11.,11.,-1.,-1.,10.,
            9., 8., 7.,-1., 6., 5., 4., 3., 2., 1.]), 2))

