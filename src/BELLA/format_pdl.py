# -*- coding: utf-8 -*-
"""
- convert_ss_to_sst
    produces a stacking sequence table from stacking sequences and a ply drop
    layout

- convert_sst_to_ss
    produces a list of stacking sequences from a stacking sequence table

- convert_ss_ref_to_reduced_sst
    retrieves a reduced stacking sequence table from a modified panel stacking
    sequence, a reduced ply drop layout scheme and the previous reduced
    stacking sequence list

- convert_ss_guide_to_sst
    retrieves a stacking sequence table from a guide stacking sequence and the
    ply drop layout scheme

- reduce_for_guide_based_blending
    returns smaller data structures for a structure with guide-based blending,
    either a reduced stacking sequence table, a reduced list of stacking
    sequences or a reduced ply drop layout

- extend_after_guide_based_blending
    returns the complete stacking sequence table or the list of stacking
    sequences for a structure with guide-based blending

- pos_in_ss_ref_to_pos_in_sst
    converts the positions of plies in the reference panel stacking sequence
    into ply positions related to the stacking sequence table

- pos_in_sst_to_pos_in_panel, pos_in_sst_to_pos_in_panel_2
    converts the positions of plies related to a stacking sequence table
    into ply positions in a specific panel stacking sequence
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

def pos_in_ss_ref_to_pos_in_sst(pos_ref, pdl_ref):
    """
    converts the positions of plies in the reference panel stacking sequence
    into ply positions related to the stacking sequence table

    INPUTS

    pos_ref: positions of the plies in the reference panel
    pdl_ref: line of the ply drop layout table corresponding to the reference
    panel
    """
    if not len(pos_ref):
        return []
    # to order the positions
    ind_sort = np.argsort(pos_ref)
    pos_ref = np.copy(pos_ref)[ind_sort]
    pos_sst = np.array((), dtype='int16')
    counter_pos = 0
    counter = 0
    for ind_pdl_ply, pdl_ply in enumerate(pdl_ref):
        if pdl_ply != -1:
            counter += 1
            if counter == pos_ref[counter_pos]:
                pos_sst = np.hstack((pos_sst, ind_pdl_ply + 1))
                counter_pos += 1
                if counter_pos == pos_ref.size:
                    break

    # to retrieve original order of the positions
    reorder = [np.where(ind_sort == ind)[0][0] for ind in range(pos_sst.size)]
    return pos_sst[reorder]


def pos_in_sst_to_pos_in_panel(pos_sst, pdl_panel):
    """
    converts the positions of plies related to a stacking sequence table
    into ply positions in a specific panel stacking sequence

    INPUTS

    - pos_sst: ordered positions of the plies related to the stacking sequence
    table
    - pdl_panel: line of the ply drop layout table corresponding to the panel
    """
#    print('pos_sst', pos_sst)
#    print('pdl_panel', pdl_panel)
    pos_panel = np.array((), dtype='int16')
    if not len(pos_sst):
        return pos_panel
    counter_pos_sst = 0
    counter_ply_panel = 0
    for ind_ply in range(0, pos_sst[-1]):
#        print('ind_ply', ind_ply)
#        print('counter_pos_sst', counter_pos_sst)
#        print('pos_sst[counter_pos_sst]', pos_sst[counter_pos_sst],
#              'ind_ply + 1', ind_ply + 1)
        if pos_sst[counter_pos_sst] == ind_ply + 1:
            if pdl_panel[ind_ply] != -1:
                counter_ply_panel += 1
                pos_panel = np.hstack((pos_panel, counter_ply_panel))
            counter_pos_sst += 1
        else:
            if pdl_panel[ind_ply] != -1:
                counter_ply_panel += 1
#        print('pos_panel', pos_panel)
    return pos_panel


def pos_in_sst_to_pos_in_panel_2(pos_sst, pdl_panel):
    """
    converts the positions of plies related to a stacking sequence table
    into ply positions in a specific panel stacking sequence

    This function is used for input ply positions not necessarily oredered,
    the output are either the ply positions in the panel, or 1e10 if the plies
    are not in the panel.

    INPUTS

    - pos_sst: positions of the plies related to the stacking sequence
    table
    - pdl_panel: line of the ply drop layout table corrsponding to the panel
    """
#    print('pos_sst', pos_sst)
#    print('pdl_panel', pdl_panel)
    pos_panel = np.array((), dtype='int16')
    if not len(pos_sst):
        return pos_panel
    for input_ply_pos in pos_sst:
        if pdl_panel[input_ply_pos - 1] == -1:
            pos_panel = np.hstack((pos_panel, 1e10))
        else:
            pos_panel = np.hstack((pos_panel, len(
                list(filter(lambda x: (x != -1), pdl_panel[:input_ply_pos])))))
    return pos_panel

def convert_sst_to_ss(ss_tab):
    """
    converts a stacking sequence table to a list of stacking sequences
    """
    liste = []
    for ind_panel in range(len(ss_tab)):
        ss_new = np.array((), dtype=int)
        for ind_ply in range(ss_tab[ind_panel].size):
            if ss_tab[ind_panel][ind_ply] != -1:
                ss_new = np.hstack((ss_new, ss_tab[ind_panel][ind_ply]))
        liste.append(ss_new)
    return liste


def convert_ss_guide_to_sst(ss_guide, pdl):
    """
    retrieves a stacking sequence table from a guide stacking sequence and the
    ply drop layout scheme
    """
    sst = -np.ones((pdl.shape[0], pdl.shape[1]), dtype='int16')
    for ind in range(ss_guide.size):
        to_change = np.where(pdl[:, ind] != -1)[0]
        for elem in to_change:
            sst[elem][ind] = ss_guide[ind]
    return sst

def convert_ss_ref_to_reduced_sst(
        ss_ref, ind_ref, reduced_pdl, reduced_ss_before):
    """
    retrieves a reduced stacking sequence table from a modified panel stacking
    sequence, a reduced ply drop layout scheme and the previous reduced
    stacking sequence list
    """
    reduced_sst = convert_ss_to_sst(reduced_ss_before, reduced_pdl)
    ind_in_ss_ref = 0
    for ind_in_sst in range(reduced_sst.shape[1]):
        if reduced_sst[ind_ref, ind_in_sst] != -1:
            to_change = np.where(reduced_sst[:, ind_in_sst] != - 1)[0]
            for index in to_change:
                reduced_sst[index, ind_in_sst] = ss_ref[ind_in_ss_ref]
            ind_in_ss_ref += 1
    return reduced_sst


def convert_ss_to_sst(ss, pdl):
    """
    retrieves a stacking sequence table from a list of stacking sequence and
    the ply drop layout scheme
    """
    return convert_ss_guide_to_sst(ss[-1], pdl)


def reduce_for_guide_based_blending(multipanel, data):
    """
    returns smaller data structures for a structure with guide-based blending,
    either a reduced stacking sequence table, a reduced list of stacking
    sequences or a reduced ply drop layout

    INPUTS
    - multipanel: multi-panel structure
    - data: data to be reduced
    """
    if isinstance(data, list):
        return [data[multipanel.reduced.ind_for_reduc[ind]] \
                for ind in range(multipanel.reduced.n_panels)]

    return np.array([data[multipanel.reduced.ind_for_reduc[ind]] \
            for ind in range(multipanel.reduced.n_panels)])


def extend_after_guide_based_blending(multipanel, reduced_ss):
    """
    returns the complete stacking sequence table or the list of stacking
    sequences for a structure with guide-based blending

    INPUTS
    - multipanel: multi-panel structure
    - reduced_ss: reduced laminate stacking sequences or stacking sequence
    table
    """
    if isinstance(reduced_ss, list):
        return [reduced_ss[multipanel.reduced.ind_panels_guide[ind]] \
                for ind in range(multipanel.n_panels)]
    return np.array([reduced_ss[multipanel.reduced.ind_panels_guide[ind]] \
            for ind in range(multipanel.n_panels)])




if __name__ == "__main__":

    print('\n*** Test for the function convert_ss_guide_to_sst ***')
    ss_guide = np.array([0, 45, 90, -45])
    pdl = np.array([[-1, 1, 2, -1], [0, 1, 2, 3]])
    print('SS_guide', ss_guide)
    print('pdl', pdl)
    print(convert_ss_guide_to_sst(ss_guide, pdl))

    print('\n*** Test for the function convert_sst_to_ss ***')
    ss_tab = np.array([
        [45, 45, -45, -45, 0, -45, -45, 90, -45, 0, -45, 0, 0, 90],
        [45, 45, -45, -45, 0, -1, -45, -45, 90, -45, -45, 0, 0, 90],
        [45, 45, -45, -1, 0, -45, -1, -45, 90, -45, -45, 0, 0, 90]])
    print_list_ss(convert_sst_to_ss(ss_tab))

    print('\n*** Test for the function pos_in_ss_ref_to_pos_in_sst ***')
#    pos_ref = [3, 5, 7]
#    reduced_pdl = np.array([
#        [0, -1, -1, 3, 4, 5, -1, -1, 8, -1, 10, 11, 12, 13],
#        [0, 1, -1, 3, -1, 5, -1, -1, 8, 9, 10, 11, 12, 13],
#        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]])
#    ind_ref = 1
#    expected_answer = [4, 9, 11]
#    print(pos_in_ss_ref_to_pos_in_sst(pos_ref, reduced_pdl[ind_ref]))
#    print('expected_answer', expected_answer)

#    pos_ref = [7, 3, 5]
#    reduced_pdl = np.array([
#        [0, -1, -1, 3, 4, 5, -1, -1, 8, -1, 10, 11, 12, 13],
#        [0, 1, -1, 3, -1, 5, -1, -1, 8, 9, 10, 11, 12, 13],
#        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]])
#    ind_ref = 1
#    expected_answer = [11, 4, 9]
#    print(pos_in_ss_ref_to_pos_in_sst(pos_ref, reduced_pdl[ind_ref]))
#    print('expected_answer', expected_answer)

    print('\n*** Test for the function pos_in_sst_to_pos_in_panel ***')
#    pos_sst = [4, 9, 11]
#    pdl_panel = [0, -1, -1, 3, 4, 5, -1, -1, 8, -1, 10, 11, 12, 13]
#    expected_answer = [2, 5, 6]
#    print(pos_in_sst_to_pos_in_panel(pos_sst, pdl_panel))
#    print('expected_answer', expected_answer)

#    pos_sst = [4, 9, 11]
#    pdl_panel = [0, -1, -1, -1, 4, 5, -1, -1, 8, -1, 10, 11, 12, 13]
#    expected_answer = [4, 5]
#    print(pos_in_sst_to_pos_in_panel(pos_sst, pdl_panel))
#    print('expected_answer', expected_answer)

#    pos_sst = [2, 3, 4]
#    pdl_panel= [0, 1, 2, -1, -1, 5, 6, 7, 8, 9, -1, -1, 12, 13]
#    expected_answer = [2]
#    print(pos_in_sst_to_pos_in_panel(pos_sst, pdl_panel))
#    print('expected_answer', expected_answer)

    print('\n*** Test for the function pos_in_sst_to_pos_in_panel_2 ***')
    pos_sst = [6, 5]
    pdl_panel= [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    expected_answer = [6, 5]
    print(pos_in_sst_to_pos_in_panel_2(pos_sst, pdl_panel))
    print('expected_answer', expected_answer)