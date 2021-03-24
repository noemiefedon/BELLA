# -*- coding: utf-8 -*-
"""
Functions used to generate manufacturable ply drop layouts with guide-based
blending

- format_ply_drops and format_ply_drops2
    format the ply drop layouts

- ply_drops_rules
    deletes the ply drop layouts that does not satisfy the ply drop guidelines

- randomly_pdl_guide
    randomly generates manufacturable ply drop layouts

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
import time
import random
import numpy as np
import scipy.special

sys.path.append(r'C:\BELLA')
from src.BELLA.parameters import Parameters
from src.BELLA.constraints import Constraints
from src.BELLA.obj_function import ObjFunction
from src.guidelines.ply_drop_spacing import calc_penalty_spacing
from src.BELLA.pdl_tools import format_ply_drops
from src.BELLA.pdl_tools import ply_drops_at_each_boundaries
from src.BELLA.pdl_tools import format_ply_drops2

from src.divers.pretty_print import print_lampam, print_ss, print_list_ss

def randomly_pdl_guide(
        boundaries,
        n_ply_drops,
        n_max,
        parameters,
        obj_func_param,
        constraints,
        multipanel,
        n_pdl_max=1,
        pdl_before=None,
        pdl_after=None,
        last_group=False,
        covering_top=False,
        covering_bottom=False,
        has_middle_ply=False,
        middle_ply_indices=np.array((), dtype='int16')):
    """
    randomly generates ply drop layouts that best satisfy the spacing and
    stacking rules within a limited time.

    INPUTS
    - n_ply_drops: list of the number of ply drops for each group compared to
    the groups thickest of the thickest laminate
    - n_max: maximum number of plies for the group
    - n_pdl_max: number of ply drop layouts asked
    - pdl_before: matrix of ply drop layouts for the group placed above
    - pdl_after: matrix of ply drop layouts for the group placed below
    - last_group: true for the last groups
    - constraints: design guidelines
    - parameters: optimiser parameters
    - multpanel: multi-panel structure
    - obj_func_param: objective function parameters
    - if covering_top = True, the top ply cannot be dropped
    - if covering_bottom = True, the bottom ply cannot be dropped
    - has_middle_ply: True if a panel has a middle ply
    """
#    print('boundaries', boundaries)
#    print('n_ply_drops', n_ply_drops)
#    print('n_max', n_max)
#    print('n_pdl_max', n_pdl_max)
#    print('last_group', last_group)
#    print('pdl_before', pdl_before)
#    print('has_middle_ply', has_middle_ply)
#    print('middle_ply_indices', middle_ply_indices)

#    # boundaries one panel to another by increasing order of thickness
#    boundaries = np.zeros((0, 2), dtype='int16')
#    for ind_panel in range(n_ply_drops.size - 1):
#        boundaries = np.vstack((
#            boundaries, np.array([ind_panel, ind_panel + 1], dtype='int16')))

    n_ply_drops_unique = np.unique(n_ply_drops)
    n_unique = n_ply_drops_unique.size
    # dictionary to retrieve indices related to the number of ply drops
    indices_unique = dict()
    for index, unique_index in enumerate(n_ply_drops_unique):
        indices_unique[unique_index] = index

    combi = []
    for drops in n_ply_drops_unique:
        # combi = list of the position that can take the ply drops per panel
        combi.append(scipy.special.comb(n_max, drops))

    n_pdl = int(min(np.product(combi), n_pdl_max))
    #print('n_pdl', n_pdl)
    # length of the group ply drop layout
    if last_group and has_middle_ply:
        n_maxx = n_max + 1
    else:
        n_maxx = n_max

    pdl_perfect = np.zeros((n_pdl, n_ply_drops.size, n_maxx), dtype=int)
    pdl_imperfect = np.zeros((n_pdl, n_ply_drops.size, n_maxx), dtype=int)
    p_spacing_imperfect = np.zeros((n_pdl,), dtype=float)

    ind_imperfect = 0
    ind_perfect = 0
    t_ini = time.time()
    elapsed_time = 0


    #print('n_pdl', n_pdl)
    while ind_perfect < n_pdl \
    and elapsed_time < parameters.time_limit_group_pdl:
        #print('ind_perfect', ind_perfect)
        #print('ind_imperfect', ind_imperfect)

        # randomly chose a pdl
        new_pdl = [[]]*n_unique
        if covering_top and covering_bottom:
            new_pdl[n_unique - 1] = random.sample(
                range(1, n_max - 1), n_ply_drops_unique[n_unique - 1])
        elif covering_top:
            new_pdl[n_unique - 1] = random.sample(
                range(1, n_max), n_ply_drops_unique[n_unique - 1])
        elif covering_bottom:
            new_pdl[n_unique - 1] = random.sample(
                range(n_max - 1), n_ply_drops_unique[n_unique - 1])
        else:
            new_pdl[n_unique - 1] = random.sample(
                range(n_max), n_ply_drops_unique[n_unique - 1])

        # for guide-based blending, not generalised blending
        for ind_panel in range(n_unique - 1)[::-1]:
            new_pdl[ind_panel] = new_pdl[ind_panel + 1][:]
            n_to_del = n_ply_drops_unique[ind_panel + 1] \
            - n_ply_drops_unique[ind_panel]
            to_del = random.sample(
                list(range(len(new_pdl[ind_panel]))), n_to_del)
            new_pdl[ind_panel] = [
                elem for ind_elem, elem in enumerate(new_pdl[ind_panel]) \
                if ind_elem not in to_del]

        # Formatting the ply drop layout in the form as in the example:
        #   [[0  1  2  3]
        #    [-1 -1  2  3]
        #    [-1 -1 -1  3]]
        # for a pdl with three panels
        #       the first panel having the four plies of index 0, 1, 2, 3
        #       the second panel having 2 plies of index 2 and 3
        #       the last panel having only the ply of index 3
        #print('new_pdl', new_pdl)
        new_pdl = format_ply_drops(new_pdl, n_max)
        # <class 'numpy.ndarray'>
#        print('new_pdl')
#        print(new_pdl)

        if last_group and has_middle_ply:
            middle = -(middle_ply_indices[:-1] != 0).astype(int)
#            print('middle', middle)
            middle = middle.reshape((new_pdl.shape[0], 1))
#            print('middle', middle)
            new_pdl = np.hstack((new_pdl, middle))


        # Formatting the ply drop layout so that a ply drop scheme is
        # associated to each panel boundary
        new_pdl = ply_drops_at_each_boundaries(
            new_pdl, n_ply_drops_unique, indices_unique, n_ply_drops)
        # <class 'numpy.ndarray'>
#        print('new_pdl')
#        print(new_pdl)

        # Application ply drop spacing and stacking rules:
        #   - Ply drops should be separated by at least min_drop plies

        # for the last groups of symmetric laminates
        if last_group and constraints.sym:
            if has_middle_ply:
                pdl_after = np.flip(np.copy(new_pdl[:, :-1]), axis=1)
                # print(new_pdl)
            else:
                pdl_after = np.flip(np.copy(new_pdl), axis=1)

        p_spacing = calc_penalty_spacing(
            pdl=new_pdl,
            pdl_before=pdl_before,
            pdl_after=pdl_after,
            multipanel=multipanel,
            obj_func_param=obj_func_param,
            constraints=constraints,
            on_blending_strip=True)

#        print('p_spacing', p_spacing)
        # <class 'numpy.ndarray'>
#        print('new_pdl1')
#        print(new_pdl)

        # Formatting the ply drop layout in the form as in the example:
        #                      [[0  1  2  3]
        #                       [-1 -1  1  2]
        #                       [-1 -1 -1  1]]
        # for a pdl:  with three panels
        #     the first panel having the four plies of index 0, 1, 2, 3
        #     the second panel having 2 plies of index 2 and 3
        #     the last panel having only the ply of index 3
        new_pdl = format_ply_drops2(new_pdl).astype(int)
        # <class 'numpy.ndarray'>
#        print('new_pdl', new_pdl)

        elapsed_time = time.time() - t_ini

        # Store the new pdl if it is perfect (no violation of manufacturing
        # constraint) or if it is among the n_pdl best unmanufacturable
        # solutions found so far
        if p_spacing == 0:
            # To remove duplicates
            is_double = False
            for ind in range(ind_perfect):
                if np.allclose(new_pdl, pdl_perfect[ind]):
                    is_double = True
                    break
            if is_double:
                continue
            #print('is_double', is_double)
            pdl_perfect[ind_perfect] = new_pdl
            ind_perfect += 1
        else:
            # To only keep the imperfect pdl with the smallest penalties
            if ind_imperfect >= n_pdl:
                if p_spacing < max(p_spacing_imperfect):
                    # To remove duplicates
                    is_double = False
                    for ind in range(ind_imperfect):
                        if np.allclose(new_pdl, pdl_imperfect[ind]):
                            is_double = True
                            break
                    if is_double:
                        continue
                    #print('is_double', is_double)
                    indexx = np.argmin(p_spacing_imperfect)
                    pdl_imperfect[indexx] = new_pdl
                    p_spacing_imperfect[indexx] = p_spacing
            else:
                # To remove duplicates
                is_double = False
                for ind in range(ind_imperfect):
                    if np.allclose(new_pdl, pdl_imperfect[ind]):
                        is_double = True
                        break
                if is_double:
                    continue
                #print('is_double', is_double)
#                print(new_pdl)
#                print(ind_imperfect)
#                print(pdl_imperfect.shape)
                pdl_imperfect[ind_imperfect] = new_pdl
                p_spacing_imperfect[ind_imperfect] = p_spacing
                ind_imperfect += 1


    # if the time limit is reached
    if elapsed_time >= parameters.time_limit_group_pdl:

        pdl_imperfect = pdl_imperfect[:ind_imperfect]
        p_spacing_imperfect = p_spacing_imperfect[:ind_imperfect]
        #print('pdl_perfect', pdl_perfect)
        #print('pdl_imperfect', pdl_imperfect)

        if not ind_imperfect + ind_perfect:
            print('n_ply_drops', n_ply_drops)
            print('n_max', n_max)
            print('min_drop', constraints.min_drop)
            print('pdl_before', pdl_before)
            print('pdl_after', pdl_after)
            raise Exception("""
No conform ply drop layout can be generated.
Too many ply drops between two adjacent panels.""")


        # add the non-manufacturable ply drop layouts
        for ind in range(ind_perfect, n_pdl_max):
            indexx = np.argmin(p_spacing_imperfect)
            pdl_perfect[ind] = pdl_imperfect[indexx]
            p_spacing_imperfect[indexx] = 10e6

        return pdl_perfect
    # if enough manufacturable ply drop layouts have been found
    return pdl_perfect
