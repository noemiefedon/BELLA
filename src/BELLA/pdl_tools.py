# -*- coding: utf-8 -*-
"""
Functions used to generate manufacturable ply drop layouts with guide-based
blending

- format_ply_drops and format_ply_drops2
    format the ply drop layouts

- ply_drops_rules
    deletes the ply drop layouts that does not satisfy the ply drop guidelines

- global_pdl_from_local_pdl
    combines group ply drop layouts to form the ply drop layouts of an entire
    laminate structure

- input_pdl_from_sst
    recovers the ply drop layout based on a stacking sequence table

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
"""
import sys
import numpy as np

sys.path.append(r'C:\BELLA')
from src.guidelines.ply_drop_spacing import calc_penalty_spacing_1ss


def ply_drops_at_each_boundaries(
        new_pdl, n_ply_drops_unique, indices_unique, n_ply_drops):
    """
    formats the ply drop layout:
        - Initially, the ply drop layout is formatted for each number of ply
        drops
        - Then, the ply drop layout is formatted for each panel boundary
    """
    pdl = np.zeros((n_ply_drops.size, new_pdl.shape[1]), int)
    for line in range(pdl.shape[0]):
        pdl[line] = new_pdl[indices_unique[n_ply_drops[line]]]
    return pdl


def ply_drops_rules(
    pdl,
    min_drop,
    boundaries,
    pdl_before=None,
    pdl_after=None,
    pdl_spacing=False):
    """
    deletes the ply drop layouts that does not satisfy the ply drop spacing and
    stacking rules.

    INPUTS

    - pdl: matrix of ply drop layouts of the current group
    - pdl_spacing to activate the ply drop spacing rule
    - pdl_before: matrix of ply drop layout of the previous group
    - pdl_after: matrix of ply drop layout of the group placed afterwards
    - min_drop: minimum number of continuous plies required between two block
    of dropped plies
    - boundaries: panels connectivity matrix, each row of the array stores the
    indices of two adjacent panels
     """
    if pdl_before is None:
        pdl_before = np.ones((pdl.shape[1],1))
    if pdl_after is None:
        pdl_after = np.ones((pdl.shape[1], 1))

    pdl_to_keep = np.ones((pdl.shape[0],), dtype=bool)
#    print('pdl_before', pdl_before)
#    print('pdl_after', pdl_after)
    if pdl_spacing:
        for ind_pdl in range(pdl.shape[0]):
            for ii1, ii2 in boundaries:
                # stack the ply drop layouts before/current/after
                layout1 = np.hstack((pdl_before[ii1],
                                     pdl[ind_pdl, ii1],
                                     pdl_after[ii1]))
                layout2 = np.hstack((pdl_before[ii2],
                                     pdl[ind_pdl, ii2],
                                     pdl_after[ii2]))
                # delete plies that does not cover none of the two adjacent
                # panels
                to_keep = [layout1[ii] >= 0 \
                         or layout1[ii] != layout2[ii] \
                             for ii in range(layout1.size)]
                layout1 = layout1[to_keep]
                layout2 = layout2[to_keep]
                # check if the ply drop spacing rule is verified
                if pdl_spacing:
                    if calc_penalty_spacing_1ss(layout1, min_drop) \
                    + calc_penalty_spacing_1ss(layout2, min_drop) > 0:
                        pdl_to_keep[ind_pdl] = False
                        break
    return pdl[pdl_to_keep]


def format_ply_drops(my_list, n_max):
    """
    formats the matrix of ply drop layout 'my_list' so that each panel is
    described with a list of 'n_max' numbers such as for the thickest panel:
        pdl_final = [0, 1, 2, 3, ..., n_max - 1]
    and for another thinner panel:
       - if a ply of index in the thicker panel belongs to the thinner panel:
           pdl_final[index] = index
       - otherwise:
           pdl_final[index] = -1
     """
    result = np.matlib.repmat(np.arange(n_max), len(my_list), 1)
    for index1, el1 in enumerate(my_list):
        for el2 in el1:
            result[index1, el2] = -1
    return result


def format_ply_drops2(ss):
    """
    formats the matrix of ply drop layout 'my_list' so that each panel is
    described with a list of 'n_max' numbers such as for the thickest panel:
        pdl_final = [0, 1, 2, 3, ..., n_max - 1]
    and for another thinner panel:
       - if a ply of index in the thicker panel belongs to the thinner panel:
           pdl_final[index] = number of the plies stacked on the panel
       - otherwise:
           pdl_final[index] = -1
     """
    for ind_panel in range(ss.shape[0]):
        index_plyLocal = 0
        for index_plyGlobal in range(ss.shape[1]):
            if ss[ind_panel, index_plyGlobal] != -1:
                ss[ind_panel, index_plyGlobal] = index_plyLocal
                index_plyLocal += 1
    return ss


def global_pdl_from_local_pdl(multipanel, sym, pdl_before_cummul,
                              pdl_after_cummul=None):
    """
    combines group ply drop layouts to form the ply drop layouts of an entire
    laminate structure

    INPUTS
    multipanel: multi-panel structure
    sym for a symmetric panel
    pdl_before_cummul and pdl_after_cummul: array of the ply drop layouts for
    the successive groups
    """
#    print('pdl_before_cummul')
#    print(pdl_before_cummul[0].shape, pdl_before_cummul[1].shape)
#    print('pdl_after_cummul')
#    print(pdl_after_cummul)
    # assemble the pdl with -1 for ply drops and 1 for non-dropped plies
    pdl = [None]*(multipanel.reduced.n_panels)
    for ind_panel in range(multipanel.reduced.n_panels):
        pdl[ind_panel] = np.array([], dtype=int)
    if sym: # for symmetric laminates
        for ind_panel in range(multipanel.reduced.n_panels):
            index_in_ss = 0
            for local_pdl in pdl_before_cummul:
                if local_pdl is not None:
                    for index_ply \
                    in range(local_pdl[ind_panel].size):
                        if local_pdl[ind_panel][index_ply] == -1:
                            pdl[ind_panel] = np.hstack((
                                pdl[ind_panel],
                                np.array([-1]).astype(int)))
                        else:
                            pdl[ind_panel] = np.hstack((
                                pdl[ind_panel],
                                np.array([1]).astype(int)))
                            index_in_ss += 1
#            if multipanel.has_middle_ply:
#                if multipanel.middle_ply[ind_panel]:
#                    pdl[ind_panel] = np.hstack((
#                        pdl[ind_panel],
#                        np.array([1]).astype(int)))
#                else:
#                    pdl[ind_panel] = np.hstack((
#                        pdl[ind_panel],
#                        np.array([-1]).astype(int)))
            pdl[ind_panel] = np.hstack((
                pdl[ind_panel], np.flip(pdl[ind_panel], axis=0)))
        pdl = np.array(pdl)
        if multipanel.has_middle_ply:
            pdl = np.delete(pdl, np.s_[pdl.shape[1] // 2], axis=1)
    else: # for asymmetric laminates
        pdl_end = [None]*(multipanel.reduced.n_panels)
        for ind_panel in range(multipanel.reduced.n_panels):
            pdl_end[ind_panel] = np.array([], dtype=int)
        for ind_panel in range(multipanel.reduced.n_panels):
            index_in_ss = 0
            index_in_ss_end = 1
            for ind_local_pdl, local_pdl in enumerate(pdl_before_cummul):
                if ind_local_pdl % 2 == 0 and pdl_before_cummul \
                and local_pdl is not None:
                    for index_ply in range(local_pdl[ind_panel].size):
                        if pdl_before_cummul[
                                ind_local_pdl][ind_panel][index_ply] == -1:
                            pdl[ind_panel] = np.hstack((
                                pdl[ind_panel],
                                np.array([-1]).astype(int)))
                        else:
                            pdl[ind_panel] = np.hstack((
                                pdl[ind_panel],
                                np.array([1]).astype(int)))
                            index_in_ss += 1
                elif ind_local_pdl % 2 == 1 and pdl_after_cummul \
                and pdl_after_cummul[ind_local_pdl] is not None:
                    for index_ply in range(pdl_after_cummul[
                            ind_local_pdl][ind_panel].size)[::-1]:
                        if pdl_after_cummul[
                                ind_local_pdl][ind_panel][index_ply] == -1:
                            pdl_end[ind_panel] = np.hstack((
                                np.array([-1]).astype(int),
                                pdl_end[ind_panel]))
                        else:
                            pdl_end[ind_panel] = np.hstack((
                                np.array([1]).astype(int),
                                pdl_end[ind_panel]))
                            index_in_ss_end += 1
            pdl[ind_panel] = np.hstack((
                pdl[ind_panel], pdl_end[ind_panel]))
    pdl = np.array(pdl)
#    print('pdl', pdl)
    # value for the non-dropped plies changed to the position of the plies
    if sym: # for symmetric laminates
        for ind_panel in range(multipanel.reduced.n_panels):
            to_add = 0
            for ind_ply in range(pdl.shape[1]):
                if pdl[ind_panel, ind_ply] != -1:
                    pdl[ind_panel, ind_ply] += to_add
                    to_add += 1
        pdl[:, (pdl.shape[1] + 1)//2:] \
        = np.flip(pdl[:, :pdl.shape[1] // 2], axis=1)
    else: # for asymmetric laminates
        for ind_panel in range(multipanel.reduced.n_panels):
            to_add = 0
            for ind_ply in range(pdl.shape[1]):
                if pdl[ind_panel, ind_ply] != -1:
                    pdl[ind_panel, ind_ply] += to_add
                    to_add += 1
    if pdl.shape[1] != multipanel.n_plies_max \
    and pdl.shape[1] != multipanel.n_plies_max + 1:
        print('pdl.shape', pdl.shape)
        raise Exception('This should not happen')
    return pdl


def input_pdl_from_sst(sst, multipanel, constraints):
    """
    recovers the ply drop layout baed on a stacking sequence table
    """
    sst_to_mod = np.copy(sst)

    if constraints.sym:
        pdl_before_cummul = [None] * 2
        for index_pdl in range(len(pdl_before_cummul)):
            pdl_before_cummul[index_pdl] = None

        if constraints.n_covering == 1:
            pdl_before_cummul[0] = np.matlib.repmat(
                np.array([0], dtype=int), multipanel.n_panels, 1)
            sst_to_mod = sst_to_mod[:, 1:sst.shape[1] // 2]
        elif constraints.n_covering == 2:
            pdl_before_cummul[0] = np.matlib.repmat(
                np.array([0, 1], dtype=int), multipanel.n_panels, 1)
            sst_to_mod = sst_to_mod[:, 2:sst.shape[1] // 2]
        else:
            sst_to_mod = sst_to_mod[:, :sst.shape[1] // 2]

        for ind_panel in range(sst_to_mod.shape[0]):
            counter = 0
            for ind_angle in range(sst_to_mod.shape[1]):
                if sst_to_mod[ind_panel, ind_angle] != - 1:
                    sst_to_mod[ind_panel, ind_angle] = counter
                    counter += 1

        pdl_before_cummul[1] = sst_to_mod

        my_pdl = global_pdl_from_local_pdl(
            multipanel, constraints.sym, pdl_before_cummul)
        return(my_pdl, pdl_before_cummul)

    pdl_before_cummul = [None] * 3
    pdl_after_cummul = [None] * 3
    for index_pdl in range(len(pdl_before_cummul)):
        pdl_before_cummul[index_pdl] = None
        pdl_after_cummul[index_pdl] = None

    if constraints.n_covering == 1:
        pdl_before_cummul[0] = np.matlib.repmat(
            np.array([0], dtype=int), multipanel.n_panels, 1)
        pdl_after_cummul[1] = np.matlib.repmat(
            np.array([0], dtype=int), multipanel.n_panels, 1)
        sst_to_mod = sst_to_mod[:, 1:-1]
    elif constraints.n_covering == 2:
        pdl_before_cummul[0] = np.matlib.repmat(
            np.array([0, 1], dtype=int), multipanel.n_panels, 1)
        pdl_after_cummul[1] = np.matlib.repmat(
            np.array([0, 1], dtype=int), multipanel.n_panels, 1)
        sst_to_mod = sst_to_mod[:, 2:-2]

    for ind_panel in range(sst_to_mod.shape[0]):
        counter = 0
        for ind_angle in range(sst_to_mod.shape[1]):
            if sst_to_mod[ind_panel, ind_angle] != - 1:
                sst_to_mod[ind_panel, ind_angle] = counter
                counter += 1

    pdl_before_cummul[2] = sst_to_mod

    my_pdl = global_pdl_from_local_pdl(
        multipanel, constraints.sym, pdl_before_cummul, pdl_after_cummul)
    return(my_pdl, pdl_before_cummul, pdl_after_cummul)



if __name__ == "__main__":
    print('\n*** Test for the function format_ply_drops ***')
#    print('Input list:\n')
#    my_list = ((), (0, 1), (0, 1, 2))
#    print(my_list, '\n')
#    print('Input for the maximum number of plies for the group:\n')
#    n_max = 5
#    print(n_max, '\n')
#    print('output:\n')
#    print(format_ply_drops(my_list, n_max))

    print('\n*** Test for the function format_ply_drops2 ***')
#    print('Input array:\n')
#    ss = np.array([[0, 1, 2, 3, 4],
#                   [-1, -1, 2, 3, 4],
#                   [-1, 1, -1, 3, -1]])
#    print(ss, '\n')
#    print('output:\n')
#    print(format_ply_drops2(ss))

    print('\n*** Test for the function ply_drops_rules ***')
#    print('Input ply drop layouts:\n')
#    pdl = np.array([[[0, 1, 2, 3, 4, 5],
#                     [-1, -1, 2, -1, 4, 5]],
#                    [[0, 1, 2, 3, 4, 5],
#                     [0, -1, 2, 3, -1, -1]],
#                     [[0, 1, 2, 3, 4, 5],
#                      [-1, 2, -1, -1, 4, 5]]], dtype=int)
#    pdl_before = np.array([[0, -1, 2, 3, 4],
#                          [0, -1, 2, 3, 1]])
#    pdl_after = None
#    min_drop = 2
#    boundaries = np.array([[0, 1]])
#    print(pdl, '\n')
#    print('Input ply drop layout of the previous group:\n')
#    print(pdl_before, '\n')
#    print('Input ply drop layout of the next group:\n')
#    print(pdl_after, '\n')
#    print(f'MinDrop  = {min_drop}\n')
#    print('output:\n')
#    print(ply_drops_rules(
#        pdl, pdl_before, pdl_after, min_drop,boundaries,
#        pdl_spacing=False))

    print('\n*** Test for the function ply_drops_at_each_boundaries ***')
#    new_pdl = np.array([[ 0, 1, 2, 3, 4, 5],
#                        [-1, 1, 2, 3, 4, 5],
#                        [ 0, 1, -1, -1, 4, 5]])
#    n_ply_drops_unique = np.array([0, 1, 2])
#    indices_unique = {0: 0, 1: 1, 2: 2}
#    n_ply_drops = np.array([0, 1, 2, 1])
#    print(ply_drops_at_each_boundaries
#         (new_pdl, n_ply_drops_unique, indices_unique, n_ply_drops))
