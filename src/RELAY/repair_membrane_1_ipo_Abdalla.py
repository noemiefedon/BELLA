# -*- coding: utf-8 -*-
"""
- calc_objA_options_1
    calculates the possible in-plane objective function values achievable by
    modifying the fibre orientations of couples of angled plies

- repair_membrane_1_ipo:
    repair for membrane properties only accounting for one panel when the
    laminate must remain balanced
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np

sys.path.append(r'C:\BELLA')
from src.divers.sorting import sortAccording
from src.RELAY.repair_10_bal import calc_ind_plies
from src.RELAY.repair_10_bal import calc_lampamA_ply_queue
from src.RELAY.repair_10_bal import calc_lampamA_options_1
from src.RELAY.repair_membrane_1_ipo import calc_objA_options_1
from src.guidelines.ten_percent_rule_Abdalla import calc_distance_Abdalla

def repair_membrane_1_ipo_Abdalla(
        ss_ini, ply_queue_ini, in_plane_coeffs,
        p_A, lampam_target, constraints):
    """
    modifies a stacking sequence to better converge towards the in-plane target
    lamination parameters. The modifications preserves the satisfaction to the
    10% rule formulated by Abdalla, to the balance requirements and to the
    damage tolerance constraints.

    First the fibre orientations of couples of angled plies situated in the
    middle of the laminate are modified to ensure convergence for the
    in-plane lamination parameters.
    Then the fibre orientations of 0 and 90 deg plies may be modified to the
    other orientation.

    INPUTS

    - ss_ini: partially retrieved stacking sequence
    - ply_queue_ini: queue of plies for innermost plies
    - in_plane_coeffs: coefficients in the in-plane objective function
    - p_A: coefficient for the proportion of the laminate thickness that can be
        modified during the repair for membrane properties
    - lampam_target: lamination parameter targets
    - constraints: design and manufacturing constraints
    - in_plane_coeffs: coefficients in the in-plane objective function
    """

    n_plies = ss_ini.size

    ss = np.copy(ss_ini)
    ply_queue = ply_queue_ini[:]


    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)
    objA = sum(in_plane_coeffs * ((lampamA - lampam_target[0:4]) ** 2))
#    print()
#    print('lampamA', lampamA)
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

    lampamA_options = calc_lampamA_options_1(n_plies, constraints)
    objA_options = calc_objA_options_1(
        lampamA, lampamA_options, lampam_target, constraints, in_plane_coeffs)

    while np.min(objA_options) + 1e-20 < objA  and objA > 1e-10:
#        print('objA', objA)
#        print('objA_options', objA_options)

        # attempts at modifying a couple of angled plies
        ind_pos_angle1, ind_pos_angle2 = np.unravel_index(
            np.argmin(objA_options, axis=None), objA_options.shape)
        angle1 = constraints.pos_angles[ind_pos_angle1]
        angle2 = constraints.pos_angles[ind_pos_angle2]
#        print('angle1', angle1, 'angle2', angle2)
        ind_angle1 = constraints.ind_angles_dict[angle1]
        ind_angle1_minus = constraints.ind_angles_dict[-angle1]
        ind_angle2 = constraints.ind_angles_dict[angle2]
        ind_angle2_minus = constraints.ind_angles_dict[-angle2]
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)

#        print('indices_per_angle', indices_per_angle)

        # if plies +-theta exist in the middle of the laminate
        if angle1 in [0, 90]:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 2:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
        else:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 1 \
            or len(indices_per_angle[ind_angle1_minus]) < 1:
                objA_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue

        # attention to not break the 10% rule
        LPs = lampamA + lampamA_options[ind_pos_angle2] - lampamA_options[
            ind_pos_angle1]
        if calc_distance_Abdalla(LPs, constraints) > 1e-10:
            objA_options[ind_angle1, ind_angle2] = 1e10
            continue

#        print('objA_options after clean', objA_options)
#        print('+-', angle1, ' plies changed into +-', angle2, 'plies')
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle[ind_angle1]', indices_per_angle[ind_angle1])
#        print('indices_per_angle[ind_angle2]', indices_per_angle[ind_angle2])

        lampamA = LPs
        objA = objA_options[ind_pos_angle1, ind_pos_angle2]
#        print()
#        print('lampamA', lampamA)
#        print('objA', objA)

        # modification of the stacking sequence
        ind_ply_1 = indices_per_angle[ind_angle1].pop(0)
        ind_ply_2 = indices_per_angle[ind_angle1_minus].pop(0)
#        print('ind_ply_1', ind_ply_1)
#        print('ind_ply_2', ind_ply_2)
#        print('ply_queue', ply_queue)

        if ind_ply_1 == 6666: # ply from the queue
            ply_queue.remove(angle1)
            ply_queue.append(angle2)
        else:
            ss[ind_ply_1] = angle2
            if constraints.sym:
                ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

        if ind_ply_2 == 6666: # ply from the queue
            if angle1 == 90:
                ply_queue.remove(90)
            else:
                ply_queue.remove(-angle1)

            if angle2 == 90:
                ply_queue.append(90)
            else:
                ply_queue.append(-angle2)
        else:
            if angle2 != 90:
                ss[ind_ply_2] = -angle2
            else:
                ss[ind_ply_2] = 90
            if constraints.sym:
                ss[ss.size - ind_ply_2 - 1] = ss[ind_ply_2]

#        lampamA_check = calc_lampamA_ply_queue(
#            ss, ss.size, ply_queue, constraints)

        ss_list.insert(0, np.copy(ss))
        ply_queue_list.insert(0, ply_queue[:])
        lampamA_list.insert(0, np.copy(lampamA))
        objA_list.insert(0, objA)

        indices_per_angle[ind_angle2].append(ind_ply_1)
        indices_per_angle[ind_angle2_minus].append(ind_ply_2)
        if constraints.sym:
            indices_per_angle[ind_angle2].sort(reverse=True)
            indices_per_angle[ind_angle2_minus].sort(reverse=True)
        else:
            sortAccording(indices_per_angle[ind_angle2], indices_to_sort)
            sortAccording(indices_per_angle[ind_angle2_minus], indices_to_sort)
            indices_per_angle[ind_angle2].reverse()
            indices_per_angle[ind_angle2_minus].reverse()

#        print('indices_per_angle', indices_per_angle)
#        print('objA', objA)
        if objA < 1e-10:
            break

        objA_options = calc_objA_options_1(
            lampamA, lampamA_options, lampam_target, constraints,
            in_plane_coeffs)

#    print('objA', objA)

    # attempt at changing a 0 deg ply into a 90 deg ply
    ind_0 = np.where(constraints.pos_angles == 0)[0][0]
    ind_90 = np.where(constraints.pos_angles == 90)[0][0]

    if indices_per_angle[constraints.ind_angles_dict[0]]:

        LPs = lampamA + (lampamA_options[ind_90] - lampamA_options[ind_0])/2
        obj_0_to_90 = sum(in_plane_coeffs*((LPs - lampam_target[0:4])**2))

#        print('obj_0_to_90', obj_0_to_90)

        if obj_0_to_90 + 1e-20 < objA \
        and calc_distance_Abdalla(LPs, constraints) == 0:
#            print('excess_10[0]', excess_10[0])
#            print('0 deg ply changed to 90 deg ply')
            objA = obj_0_to_90
            lampamA += (lampamA_options[ind_90] - lampamA_options[ind_0])/2
            ind_ply_1 = indices_per_angle[constraints.index0].pop(0)
            if ind_ply_1 == 6666: # ply from the queue
                ply_queue.remove(0)
                ply_queue.append(90)
            else:
                ss[ind_ply_1] = 90
                if constraints.sym:
                    ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

            ss_list.insert(0, np.copy(ss))
            ply_queue_list.insert(0, ply_queue[:])
            lampamA_list.insert(0, np.copy(lampamA))
            objA_list.insert(0, objA)
#            print()
#            print('lampamA', lampamA)
#            print('objA', objA)

            return ss_list, ply_queue_list, lampamA_list, objA_list

    # attempt at changing a 90 deg ply into a 0 deg ply
    if indices_per_angle[constraints.ind_angles_dict[90]]:
        LPs = lampamA + (lampamA_options[ind_0] - lampamA_options[ind_90])/2
        obj_90_to_0 = sum(in_plane_coeffs * ((LPs - lampam_target[0:4])**2))
#        print('obj_90_to_0', obj_90_to_0)

        if obj_90_to_0 + 1e-20 < objA\
        and calc_distance_Abdalla(LPs, constraints) == 0:
#            print('90 deg ply changed to 0 deg ply')
            objA = obj_90_to_0
            lampamA += (lampamA_options[ind_0] - lampamA_options[ind_90])/2
            ind_ply_1 = indices_per_angle[constraints.index90].pop(0)
            if ind_ply_1 == 6666: # ply from the queue
                ply_queue.remove(90)
                ply_queue.append(0)
            else:
                ss[ind_ply_1] = 0
                if constraints.sym:
                    ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]

            ss_list.insert(0, np.copy(ss))
            ply_queue_list.insert(0, ply_queue[:])
            lampamA_list.insert(0, np.copy(lampamA))
            objA_list.insert(0, objA)
#            print()
#            print('lampamA', lampamA)
#            print('objA', objA)

    return ss_list, ply_queue_list, lampamA_list, objA_list
