# -*- coding: utf-8 -*-
"""
- repair_membrane_1_no_ipo:
    repair for membrane properties only accounting for one panel when the
    laminate does not have to remain balanced
    """
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import math as ma
import numpy as np

sys.path.append(r'C:\BELLA')
from src.divers.sorting import sortAccording
from src.RELAY.repair_10_bal import calc_ind_plies
from src.RELAY.repair_10_bal import calc_lampamA_ply_queue
from src.RELAY.repair_membrane_1_no_ipo import calc_objA_options_3
from src.RELAY.repair_10_bal import calc_lampamA_options_3
from src.guidelines.ten_percent_rule_Abdalla import calc_distance_Abdalla

def repair_membrane_1_no_ipo_Abdalla(
        ss_ini, ply_queue_ini, in_plane_coeffs,
        p_A, lampam_target, constraints):
    """
    modifies a stacking sequence to better converge towards the in-plane target
    lamination parameters. The modifications preserves the satisfaction to the
    10% rule, to the balance requirements and to the damage tolerance
    constraints.

    The fibre orientations are modified one by one.

    INPUTS

    - ss_ini: partially retrieved stacking sequence
    - ply_queue_ini: queue of plies for innermost plies
    - in_plane_coeffs: coefficients in the in-plane objective function
    - p_A: coefficient for the proportion of the laminate thickness that can be
        modified during the repair for membrane properties
    - lampam_target: lamination parameter targets
    - constraints: design and manufacturing constraints
    """
    n_plies = ss_ini.size

    ss = np.copy(ss_ini)
    ply_queue = ply_queue_ini[:]

    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)
    objA = sum(in_plane_coeffs * ((lampamA - lampam_target[0:4]) ** 2))
#    print('objA', objA)

    ss_list = [np.copy(ss)]
    ply_queue_list = [ply_queue[:]]
    lampamA_list = [lampamA]
    objA_list = [objA]

    indices_1, indices_per_angle = calc_ind_plies(
        ss, n_plies, ply_queue, constraints, p_A)
    indices_to_sort = list(indices_1)
    indices_to_sort.insert(0, -1)
#    print('indices_1', list(indices_1))
#    print('indices_per_angle', list(indices_per_angle))
#    print('indices_to_sort', indices_to_sort)

    lampamA_options = calc_lampamA_options_3(n_plies, constraints)
    objA_options = calc_objA_options_3(
        lampamA, lampamA_options, lampam_target, constraints, in_plane_coeffs)
#    print('objA_options', objA_options)

    while np.min(objA_options) + 1e-20 < objA and objA > 1e-10:
        # attempts at modifying a couple of angled plies
        ind_angle1, ind_angle2 = np.unravel_index(
            np.argmin(objA_options, axis=None), objA_options.shape)
        angle1 = constraints.set_of_angles[ind_angle1]
        angle2 = constraints.set_of_angles[ind_angle2]
#        print('test  angle1', angle1, 'to angle2', angle2)
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle', indices_per_angle)

        # if no ply to be deleted
        if len(indices_per_angle[ind_angle1]) < 1:
            objA_options[ind_angle1, ind_angle2] = 1e10
            continue

        # attention to not break the 10% rule
        LPs = lampamA + lampamA_options[ind_angle2] \
        - lampamA_options[ind_angle1]
        if calc_distance_Abdalla(LPs, constraints) > 1e-10:
            objA_options[ind_angle1, ind_angle2] = 1e10
            continue

#        print(angle1, ' plies changed into ', angle2, 'plies')
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle[ind_angle1]', indices_per_angle[ind_angle1])
#        print('indices_per_angle[ind_angle2]', indices_per_angle[ind_angle2])

        lampamA = LPs
        objA = objA_options[ind_angle1, ind_angle2]

        # modification of the stacking sequence
        ind_ply_1 = indices_per_angle[ind_angle1].pop(0)
#        print('ind_ply_1', ind_ply_1)

        if ind_ply_1 == 6666: # ply from the queue
            ply_queue.remove(angle1)
            ply_queue.append(angle2)
        else:
            ss[ind_ply_1] = angle2
            if constraints.sym:
                ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

        ss_list.insert(0, np.copy(ss))
        ply_queue_list.insert(0, ply_queue[:])
        lampamA_list.insert(0, np.copy(lampamA))
        objA_list.insert(0, objA)

        indices_per_angle[ind_angle2].append(ind_ply_1)
        if constraints.sym:
            indices_per_angle[ind_angle2].sort(reverse=True)
        else:
            sortAccording(indices_per_angle[ind_angle2], indices_to_sort)
            indices_per_angle[ind_angle2].reverse()

#        print('indices_per_angle', indices_per_angle)
#        print('objA', objA)
        if objA < 1e-10:
            break

        objA_options = calc_objA_options_3(
            lampamA, lampamA_options, lampam_target, constraints,
            in_plane_coeffs)
#        print('objA_options', objA_options)

    return ss_list, ply_queue_list, lampamA_list, objA_list


