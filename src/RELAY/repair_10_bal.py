# -*- coding: utf-8 -*-
"""
repair for 10% rule and balance

- repair_10_bal
    repairs a laminate regarding the 10% rule and balance

- calc_mini_10:
    returns the minimum number of plies in the 0/90/+45/-45 directions to
    satisfy the 10% rule

- calc_current_10_2:
    returns the current number of plies in the 0/90/+45/-45 directions

- is_equal
    returns True if the set of partial stacking sequence + ply queue matches
    the initial stacking sequence
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import math as ma
import numpy as np

sys.path.append(r'C:\BELLA')
from src.divers.sorting import sortAccording
from src.BELLA.constraints import Constraints
from src.divers.pretty_print import print_ss
from src.guidelines.ten_percent_rule_Abdalla import calc_distance_Abdalla


def repair_10_bal(ss_ini, mini_10, constraints):
    """
    repairs a laminate regarding the 10% rule and balance
    """
    if not (constraints.rule_10_percent and constraints.rule_10_Abdalla):
        ss, ply_queue = repair_10_bal_2(ss_ini, mini_10, constraints)
        return ss, ply_queue

#    print('initial')
#    print_ss(ss_ini)

    ## repair for balance
    mini_10 = calc_mini_10(constraints, ss_ini.size)
    ss, ply_queue = repair_10_bal_2(ss_ini, mini_10, constraints)

    ## repair for 10% rule
    if constraints.bal:
        ss, ply_queue = repair_10_Abdalla_ipo(ss, ply_queue, constraints)
    else:
        ss, ply_queue = repair_10_Abdalla_no_ipo(ss, ply_queue, constraints)

    return ss, ply_queue


def repair_10_Abdalla_ipo(ss, ply_queue, constraints):
    """
    repairs a balanced laminate regarding the 10% rule of Abdalla

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    """
    n_plies = ss.size
    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)
    dist_10 = calc_distance_Abdalla(lampamA, constraints)

    indices_1, indices_per_angle = calc_ind_plies(
        ss, n_plies, ply_queue, constraints)

    indices_to_sort = list(indices_1)
    indices_to_sort.insert(0, -1)
#    print('indices_1', list(indices_1))
#    print('indices_per_angle', list(indices_per_angle))
#    print('indices_to_sort', indices_to_sort)

    lampamA_options = calc_lampamA_options_1(n_plies, constraints)

    dist_10_options = calc_dist_10_options_1(
        lampamA, lampamA_options, constraints)


    while dist_10 > 1e-10:
#        print('dist_10', dist_10)
#        print('dist_10_options', dist_10_options)

        # attempts at modifying a couple of angled plies
        ind_pos_angle1, ind_pos_angle2 = np.unravel_index(
            np.argmin(dist_10_options, axis=None), dist_10_options.shape)
        angle1 = constraints.pos_angles[ind_pos_angle1]
        angle2 = constraints.pos_angles[ind_pos_angle2]
#        print('angle1', angle1, 'angle2', angle2)
        ind_angle1 = constraints.ind_angles_dict[angle1]
        ind_angle1_minus = constraints.ind_angles_dict[-angle1]
        ind_angle2 = constraints.ind_angles_dict[angle2]
        ind_angle2_minus = constraints.ind_angles_dict[-angle2]
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)

#        print('indices_per_angle', indices_per_angle)

        # if no couple of plies to be deleted
        if angle1 in [0, 90]:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 2:
                dist_10_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue
        else:
            # if no couple of plies to be deleted
            if len(indices_per_angle[ind_angle1]) < 1 \
            or len(indices_per_angle[ind_angle1_minus]) < 1:
                dist_10_options[ind_pos_angle1, ind_pos_angle2] = 1e10
                continue

#        print('dist_10_options after clean', dist_10_options)
#        print('+-', angle1, ' plies changed into +-', angle2, 'plies')
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle[ind_angle1]', indices_per_angle[ind_angle1])
#        print('indices_per_angle[ind_angle2]', indices_per_angle[ind_angle2])

        lampamA += lampamA_options[ind_pos_angle2] - lampamA_options[
            ind_pos_angle1]
        dist_10 = dist_10_options[ind_pos_angle1, ind_pos_angle2]
#        print()
#        print('lampamA', lampamA)
#        print('dist_10', dist_10)

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
#        print('dist_10', dist_10)
        if dist_10 < 1e-10:
            break

        dist_10_options = calc_dist_10_options_1(
            lampamA, lampamA_options, constraints)

#    print('dist_10', dist_10)

    # attempt at changing a 0 deg ply into a 90 deg ply
    ind_0 = np.where(constraints.pos_angles == 0)[0][0]
    ind_90 = np.where(constraints.pos_angles == 90)[0][0]

    if indices_per_angle[constraints.ind_angles_dict[0]]:

        dist_10_0_to_90 = calc_distance_Abdalla(
            lampamA + (lampamA_options[ind_90] - lampamA_options[ind_0])/2,
            constraints)
#        print('dist_10_0_to_90', dist_10_0_to_90)

        if dist_10_0_to_90 + 1e-20 < dist_10:
#            print('excess_10[0]', excess_10[0])
#            print('0 deg ply changed to 90 deg ply')
            dist_10 = dist_10_0_to_90
            lampamA += (lampamA_options[ind_90] - lampamA_options[ind_0])/2
            ind_ply_1 = indices_per_angle[constraints.index0].pop(0)
            if ind_ply_1 == 6666: # ply from the queue
                ply_queue.remove(0)
                ply_queue.append(90)
            else:
                ss[ind_ply_1] = 90
                if constraints.sym:
                    ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]
#            print('lampamA', lampamA)
#            print('dist_10', dist_10)

            return ss, ply_queue

    # attempt at changing a 90 deg ply into a 0 deg ply
    if indices_per_angle[constraints.ind_angles_dict[90]]:

        dist_10_90_to_0 = calc_distance_Abdalla(
            lampamA + (lampamA_options[ind_0] - lampamA_options[ind_90])/2,
            constraints)
#        print('dist_10_90_to_0', dist_10_90_to_0)

        if dist_10_90_to_0 + 1e-20 < dist_10:
#            print('90 deg ply changed to 0 deg ply')
            dist_10 = dist_10_90_to_0
            lampamA += (lampamA_options[ind_0] - lampamA_options[ind_90])/2
            ind_ply_1 = indices_per_angle[constraints.index90].pop(0)
            if ind_ply_1 == 6666: # ply from the queue
                ply_queue.remove(90)
                ply_queue.append(0)
            else:
                ss[ind_ply_1] = 0
                if constraints.sym:
                    ss[ss.size - ind_ply_1 - 1] = ss[ind_ply_1]
#            print('lampamA', lampamA)
#            print('dist_10', dist_10)

    return ss, ply_queue


def repair_10_Abdalla_no_ipo(ss, ply_queue, constraints):
    """
    repairs a non-balanced laminate regarding the 10% rule of Abdalla

    INPUTS

    - ss: partially retrieved stacking sequence
    - ply_queue: queue of plies for innermost plies
    - constraints: design and manufacturing constraints
    """
    n_plies = ss.size
    lampamA = calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints)
    dist_10 = calc_distance_Abdalla(lampamA, constraints)
#    print('dist_10', dist_10)

    indices_1, indices_per_angle = calc_ind_plies(
        ss, n_plies, ply_queue, constraints)
    indices_to_sort = list(indices_1)
    indices_to_sort.insert(0, -1)
#    print('indices_1', list(indices_1))
#    print('indices_per_angle', list(indices_per_angle))
#    print('indices_to_sort', indices_to_sort)

    lampamA_options = calc_lampamA_options_3(n_plies, constraints)
    dist_10_options = calc_dist_10_options_3(
        lampamA, lampamA_options, constraints)
#    print('dist_10_options', dist_10_options)

    while dist_10 > 1e-10:
        # attempts at modifying a couple of angled plies
        ind_angle1, ind_angle2 = np.unravel_index(
            np.argmin(dist_10_options, axis=None), dist_10_options.shape)
        angle1 = constraints.set_of_angles[ind_angle1]
        angle2 = constraints.set_of_angles[ind_angle2]
#        print('test  angle1', angle1, 'to angle2', angle2)
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle', indices_per_angle)

        # if no ply to be deleted
        if len(indices_per_angle[ind_angle1]) < 1:
            dist_10_options[ind_angle1, ind_angle2] = 1e10
            continue

#        print(angle1, ' plies changed into ', angle2, 'plies')
#        print('ind_angle1', ind_angle1, 'ind_angle2', ind_angle2)
#        print('indices_per_angle[ind_angle1]', indices_per_angle[ind_angle1])
#        print('indices_per_angle[ind_angle2]', indices_per_angle[ind_angle2])

        lampamA += lampamA_options[ind_angle2] - lampamA_options[ind_angle1]
        dist_10 = dist_10_options[ind_angle1, ind_angle2]

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

        indices_per_angle[ind_angle2].append(ind_ply_1)
        if constraints.sym:
            indices_per_angle[ind_angle2].sort(reverse=True)
        else:
            sortAccording(indices_per_angle[ind_angle2], indices_to_sort)
            indices_per_angle[ind_angle2].reverse()

#        print('indices_per_angle', indices_per_angle)
#        print('dist_10', dist_10)
        if dist_10 < 1e-10:
            break

        dist_10_options = calc_dist_10_options_3(
            lampamA, lampamA_options, constraints)
#        print('dist_10_options', dist_10_options)

    return ss, ply_queue


def repair_10_bal_2(ss_ini, mini_10, constraints):
    """
    repairs a laminate regarding the balance guideline and the 10% rule applied
    on ply counts
    """

#    if not constraints.ipo and constraints.percent_45_135 != 0:
#        raise Exception("""
#Repair for 10% rule not implemented for laminates with no balance requirements
#and a limit percentage for the ply orientated in the combined +-45 direction!
#""")

#    print('initial')
#    print_ss(ss_ini)

    ss = np.copy(ss_ini)
    if not constraints.ipo and not constraints.rule_10_percent:
        return ss, []

    if constraints.sym and ss.size % 2 and ss[ss.size // 2] not in {0, 90}:
        ss[ss.size // 2] = 0

    if constraints.sym:
        ind_plies = np.array(range(0, ss.size // 2))
    else:
        ind_plies = np.arange(ss.size)
        beginning = np.copy(ind_plies[0:ind_plies.size // 2])
        ending = np.copy(ind_plies[ind_plies.size // 2:ss.size][::-1])
        ind_plies = np.zeros((ss.size,), int)
        ind_plies[::2] = ending
        ind_plies[1::2] = beginning
#    print('ind_plies', list(ind_plies))

    ply_queue = []
    if constraints.rule_10_percent:
        for elem in range(ma.ceil(mini_10[0])):
            ply_queue.append(0)
        for elem in range(ma.ceil(mini_10[1])):
            ply_queue.append(90)
        for elem in range(ma.ceil(mini_10[2])):
            ply_queue.append(45)
        for elem in range(ma.ceil(mini_10[3])):
            ply_queue.append(-45)
#    print('initial ply queue', ply_queue)

    if constraints.rule_10_percent and constraints.percent_45_135:
        missing_extra_45_135 = ma.ceil(mini_10[4]) \
        - ma.ceil(mini_10[2]) - ma.ceil(mini_10[3])
    else:
        missing_extra_45_135 = 0

    counter_remaining_plies = ind_plies.size

    change = False
    for counter, ind_ply in enumerate(ind_plies):
#        print()
#        print('ind_ply', ind_ply)
#        print('new_angle', ss[ind_ply])
#        print('ply_queue', ply_queue)
#        print_ss(ss)

        ply_queue_before = ply_queue[:]
        new_angle = ss[ind_ply]
#        print('ind_ply', ind_ply, 'new_angle', new_angle)
        counter_remaining_plies -= 1

        if new_angle in ply_queue:
            ply_queue.remove(new_angle)
        else:
            if constraints.ipo and not new_angle in (0, 90):
                ply_queue.append(-new_angle)
            if not constraints.ipo and new_angle in (45, -45):
                missing_extra_45_135 = max(0, missing_extra_45_135 - 1)

#        print('ply_queue', ply_queue, len(ply_queue))
        if counter_remaining_plies < len(ply_queue) + missing_extra_45_135:
            change = True
            last_ply = counter
            ply_queue = ply_queue_before
            for ind_ply in ind_plies[counter:]:
                ss[ind_ply] = 666
                if constraints.sym:
                    ss[ss.size - ind_ply - 1] = 666
            break
#        print_ss(ss)
#        print('ply_queue', ply_queue)

#    print('last_ply', last_ply)
#    print('ply_queue', ply_queue, len(ply_queue))

#    print_ss(ss)

    for ind in range(missing_extra_45_135):
        if ind % 2:
            ply_queue.append(45)
        else:
            ply_queue.append(-45)

    if change and last_ply + len(ply_queue) != len(ind_plies):
#        ply_queue_2 = ply_queue[:]
#        ply_queue_2.append(90)
#        if (constraints.sym \
#            and np.isclose(np.sort(np.array(2*ply_queue_2)),
#                           np.sort(ss_ini[ss == 666])).all()) \
#        or (not constraints.sym \
#            and np.isclose(np.sort(np.array(ply_queue_2)),
#                           np.sort(ss_ini[ss == 666])).all()):
#            ply_queue.append(90)
#        else:
        ply_queue.append(0)

    return ss, ply_queue


def calc_dist_10_options_1(lampamA, lampamA_options, constraints):
    """
    calculates the possible distances away from the LP feasible region for the
    10% of Abdalla achievable by modifying the fibre orientations of couples of
    angled plies in a balanced laminate

    dist_10_options[ind_pos_angle1, ind_pos_angle2] for +-angle1 plies changed
    to +-angle2 plies

    The plies to be modified are chosen among the innermost layers to
    reduce the out-of-plane modifications
    """
    dist_10_options = 1e10*np.ones((constraints.pos_angles.size,
                                    constraints.pos_angles.size), float)
    for ind_pos_angle1 in range(constraints.pos_angles.size):
        for ind_pos_angle2 in range(constraints.pos_angles.size):
            if ind_pos_angle1 == ind_pos_angle2:
                continue
            LPs = lampamA - lampamA_options[ind_pos_angle1] \
            + lampamA_options[ind_pos_angle2]
            dist_10_options[ind_pos_angle1, ind_pos_angle2] \
            = calc_distance_Abdalla(LPs, constraints)

    return dist_10_options

def calc_dist_10_options_3(lampamA, lampamA_options, constraints):
    """
    calculates the possible in-plane objective function values achievable by
    modifying one fibre orientation in a non-balanced laminate

    dist_10_options[ind_pos_angle1, ind_pos_angle2] for angle1 ply changed
    to angle2 plies
    """
    dist_10_options = 1e10*np.ones((constraints.n_set_of_angles,
                                    constraints.n_set_of_angles), float)
    for ind_pos_angle1 in range(constraints.n_set_of_angles):
        for ind_pos_angle2 in range(constraints.n_set_of_angles):
            if ind_pos_angle1 == ind_pos_angle2:
                continue
            LPs = lampamA - lampamA_options[ind_pos_angle1] \
            + lampamA_options[ind_pos_angle2]
            dist_10_options[ind_pos_angle1, ind_pos_angle2] \
            = calc_distance_Abdalla(LPs, constraints)
    return dist_10_options


def calc_mini_10(constraints, n_plies):
    """
    returns the minimum number of plies in the 0/90/+45/-45/+-45 directions to
    satisfy the 10% rule (array)

    INPUTS

        ss: stacking sequence (array)
        constraints: constraints (instance of the class Constraints)
    """
    mini_10 = np.zeros((5,), float)
    mini_10[0] = ma.ceil(constraints.percent_0 * n_plies)
    mini_10[1] = ma.ceil(constraints.percent_90 * n_plies)
    mini_10[2] = ma.ceil(constraints.percent_45 * n_plies)
    mini_10[3] = ma.ceil(constraints.percent_135 * n_plies)
    mini_10[4] = ma.ceil(constraints.percent_45_135 * n_plies)
    if constraints.ipo:
        mini_10[2] = max(mini_10[2], mini_10[3])
        if mini_10[4] % 2:
            mini_10[4] += 1
            mini_10[4] = max(mini_10[4], 2 * mini_10[2])
        mini_10[2] = max(mini_10[2], mini_10[4] // 2)
        mini_10[3] = mini_10[2]

    if constraints.sym:
        mini_10 /= 2
        # middle ply can only be oriented at 0 or 90 degrees
        if n_plies % 2:
            mini_10[2:] = np.ceil(mini_10[2:])
        else:
            mini_10 = np.ceil(mini_10)

        if constraints.ipo:
            if mini_10[4] % 2:
                mini_10[4] += 1
                mini_10[4] = max(mini_10[4], 2 * mini_10[2])
            mini_10[2] = max(mini_10[2], mini_10[4] // 2)
            mini_10[3] = mini_10[2]

    return mini_10

def is_equal(ss, ply_queue, ss_ini, sym):
    """
    returns True if the set of partial stacking sequence + ply queue matches
    the initial stacking sequence
    """
    if not np.isclose(ss[ss != 666], ss_ini[ss != 666] ).all():
        return False

    if sym:
        if not np.isclose(np.sort(np.array(2*ply_queue)),
                          np.sort(ss_ini[ss == 666])).all():
            return False
    else:
        if not np.isclose(np.sort(np.array(ply_queue)),
                          np.sort(ss_ini[ss == 666])).all():
            return False
    return True

def calc_ind_plies(ss, n_plies, ply_queue, constraints, p_A=100):
    """
    makes:
        - a list of all ply indices which can be modified during the refinement
        for membrane properties, sorted by starting with innermost plies
        - a list of the ply indices in each fibre direction which can be
        modified during the refinement for membrane properties, sorted by
        starting with innermost plies

    Notes:
        - al lplies from the queue of plies are included.
        - middle plies of symmetric laminates are not included.
        - the rest of the plies are included only if they are part of the
        inner part of laminate representing p_A %
        of the overall laminate thickness.
    """
    ind_min = ma.floor(
        (1 - p_A/100)*(n_plies/2))
#    print('ind_min', ind_min)
    if constraints.sym:
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                indices_1 = range(max(2, ind_min), n_plies // 2)[::-1]
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                indices_1 = range(max(2, ind_min), n_plies // 2)[::-1]
            else:
                indices_1 = range(max(1, ind_min), n_plies // 2)[::-1]
        else:
            indices_1 = range(ind_min, n_plies // 2)[::-1]
    else:
        if constraints.dam_tol:
            if hasattr(constraints, 'dam_tol_rule') \
            and constraints.dam_tol_rule in {2, 3}:
                ind_1 = list(range(max(2, ind_min), n_plies // 2)[::-1])
                ind_2 = list(range(
                    n_plies // 2, min(n_plies - ind_min, n_plies - 2)))
            elif not hasattr(constraints, 'dam_tol_rule') and \
            constraints.n_plies_dam_tol == 2:
                ind_1 = list(range(max(2, ind_min), n_plies // 2)[::-1])
                ind_2 = list(range(
                    n_plies // 2, min(n_plies - ind_min, n_plies - 2)))
            else:
                ind_1 = list(range(max(1, ind_min), n_plies // 2)[::-1])
                ind_2 = list(range(
                    n_plies // 2, min(n_plies - ind_min, n_plies - 1)))
        else:
            ind_1 = list(range(ind_min, n_plies // 2)[::-1])
            ind_2 = list(range(n_plies // 2, n_plies - ind_min))
        #print('ind_1', ind_1, 'ind_2', ind_2)
        indices_1 = np.zeros((len(ind_1) + len(ind_2),), 'int16')
        indices_1[::2] = ind_2
        indices_1[1::2] = ind_1
#    print('indices_1', list(indices_1))

    indices_per_angle = []
    for ind_angle in range(constraints.n_set_of_angles):
        indices_per_angle.append([])

    for ind_ply_1 in indices_1:
        if ss[ind_ply_1] != 666:
            indices_per_angle[
                constraints.ind_angles_dict[ss[ind_ply_1]]].append(ind_ply_1)

    for angle in ply_queue:
        indices_per_angle[constraints.ind_angles_dict[angle]].insert(0, 6666)

    return indices_1, indices_per_angle

def calc_lampamA_options_1(n_plies, constraints):
    """
    calculates the elementary changes of in-plane lamination parameters
    when modifying the fibre orientations of couples of angled plies
    """
    lampamA_options = np.empty((constraints.pos_angles.size, 4), float)
    for ind_angle, angle in enumerate(constraints.pos_angles):
        lampamA_options[ind_angle] = np.copy(constraints.cos_sin[
            constraints.ind_angles_dict[angle]]).reshape(4)
        lampamA_options[ind_angle] += constraints.cos_sin[
            constraints.ind_angles_dict[-angle]].reshape(4)
    if not constraints.sym:
        lampamA_options *= (1 / n_plies)
    else:
        lampamA_options *= (2 / n_plies)
    return lampamA_options

def calc_lampamA_options_3(n_plies, constraints):
    """
    calculates the elementary changes of in-plane lamination parameters
    when modifying a fibre orientation
    """
    lampamA_options = np.empty((constraints.n_set_of_angles, 4), float)
    for ind_angle, angle in enumerate(constraints.set_of_angles):
        lampamA_options[ind_angle] = np.copy(
            constraints.cos_sin[ind_angle]).reshape(4)
    if not constraints.sym:
        lampamA_options *= 1 / n_plies
    else:
        lampamA_options *= 2 / n_plies
    return lampamA_options

def calc_lampamA_ply_queue(ss, n_plies, ply_queue, constraints):
    """
    calculates in-plane lamination parameters based on a partially retrieved
    stacking sequence and a list of plies for the innermost plies whose
    positions are left to be determined
    """
    cos_sin = np.zeros((4,), float)

    if not constraints.sym:
        for angle in ss:
            if angle != 666:
                cos_sin += constraints.cos_sin[
                    constraints.ind_angles_dict[angle]].reshape((4, ))

        for angle in ply_queue:
            cos_sin += constraints.cos_sin[
                constraints.ind_angles_dict[angle]].reshape((4, ))

        return (1 / n_plies) * cos_sin

    for angle in ss[:np.size(ss) // 2]:
        if angle != 666:
            cos_sin += constraints.cos_sin[
                constraints.ind_angles_dict[angle]].reshape((4, ))

    for angle in ply_queue:
        cos_sin += constraints.cos_sin[
            constraints.ind_angles_dict[angle]].reshape((4, ))

    if np.size(ss) % 2:
        cos_sin += 0.5 * constraints.cos_sin[
            constraints.ind_angles_dict[ss[n_plies // 2]]].reshape((4,))

    return (2 / n_plies) * cos_sin

if __name__ == "__main__":
    constraints = Constraints(
        sym=True,
        bal=True,
        dam_tol=False,
        rule_10_percent=True,
        n_contig=4,
        percent_0=10,
        percent_45=0,
        percent_90=10,
        percent_135=0,
        percent_45_135=10,
        set_of_angles=[0, 45, -45, 90])


    print('\n\n*** Test for the function is_equal ***')
    ss_ini = np.array([
        -45, 45, 60, 15, -15, 60, 30, 45, 0])
    ss = np.array([
        -45, 45, 60, 15, -15, 60, 30, 45, 0])
    ply_queue = []
    print(is_equal(ss, ply_queue, ss_ini, constraints.sym))


    print('\n\n*** Test for the function calc_mini_10 ***')
    mini_10 = calc_mini_10(constraints, 40)
    print('\nmini_10', mini_10)
    n_45 = ma.ceil(mini_10[2])
    n_135 = ma.ceil(mini_10[3])
    n_45_135 = ma.ceil(mini_10[4])
    if constraints.rule_10_percent and constraints.percent_45_135:
        missing_extra_45_135 = ma.ceil(mini_10[4]) \
        - ma.ceil(mini_10[2]) - ma.ceil(mini_10[3])
    else:
        missing_extra_45_135 = 0
    print('n_45', n_45)
    print('n_135', n_135)
    print('n_45_135', n_45_135)
    print('missing_extra_45_135', missing_extra_45_135)


    print('\n*** Test for the function repair_10_bal***')
    ss_ini = np.array([60, 45, 60, 0], int)
    ss_ini = np.array([60, 45, 60, 15, -15, 30, 45], int)
    ss_ini = np.array([-45, 45, -45, 0, 0, 0, -45, 90, 45, 45], int)
    if constraints.sym:
        ss_ini = np.hstack((ss_ini, np.flip(ss_ini)))
    print('\nInitial stacking sequence')
    print_ss(ss_ini, 2000)
    print('ss_ini.zize', ss_ini.size)
    mini_10 = calc_mini_10(constraints, ss_ini.size)
    print('mini_10', mini_10)
    ss, ply_queue = repair_10_bal(ss_ini, mini_10, constraints)
    print('\nSolution stacking sequence')
    print_ss(ss, 2000)
    print('ply_queue', ply_queue)

    print('\n*** Test for the function calc_ind_plies ***')
    constraints = Constraints(
        sym=True,
        bal=True,
        dam_tol=False,
        rule_10_percent=True,
        percent_0=10,
        percent_45=5,
        percent_90=10,
        percent_135=5,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    p_A = 50
    ss = np.array([
      45,  90,  45, 0, -45, 0, 666, 666], int)
    ss = np.hstack((ss, np.flip(ss)))
    ply_queue = [90, -45]
    n_plies = ss.size
    indices_1, indices_per_angle = calc_ind_plies(
        ss, n_plies, ply_queue, constraints, p_A)
    print('indices_1', list(indices_1))
    print('indices_per_angle', indices_per_angle)

    print('\n*** Test for the function calc_lampamA_ply_queue ***')
    constraints = Constraints(
        sym=False,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    ply_queue = [90, 90, 90, -45, -45, -45, 90, 90, 90, -45, -45, -45]
    n_plies = ss.size
    lampamA = calc_lampamA_ply_queue(
        ss, n_plies, ply_queue, constraints)

    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   90, 90, 90, -45, -45, -45, 90, 90, 90, -45, -45, -45,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    lampamA_check = calc_lampam(ss, constraints)[0:4]
    print('lampamA', lampamA)
    print('lampamA_check', lampamA_check)

    print()
    constraints = Constraints(
        sym=True,
        set_of_angles=[0, 45, 30, -30, -45, 60, -60, 90])
    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666, 666,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    ply_queue = [90, 90, 90, -45, -45, -45]
    n_plies = ss.size
    lampamA = calc_lampamA_ply_queue(
        ss, n_plies, ply_queue, constraints)

    ss = np.array([45,   0,  45,   0, -45,  45,  45,   0,  45,   0, -45,  45,  45, -45, -45, -45,  45, -45,  45,
                   90, 90, 90, -45, -45, -45, 90, 90, 90, -45, -45, -45,
                   45, -45,  45, -45, -45, -45,  45,  45, -45,   0,  45,   0,  45,  45, -45,   0,  45,   0,  45], int)
    lampamA_check = calc_lampam(ss, constraints)[0:4]
    print('lampamA', lampamA)
    print('lampamA_check', lampamA_check)
