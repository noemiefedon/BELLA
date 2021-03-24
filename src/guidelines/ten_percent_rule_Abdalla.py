# -*- coding: utf-8 -*-
"""
Application of the 10% rule for the design of a laminate

- display_ply_counts
    displays the ply counts in each fibre direction for a laminate lay-up

- is_ten_percent_rule
    returns True for a panel stacking sequence satisfying the 10% rule,
    False otherwise

- calc_penalty_10_ss and calc_penalty_10_pc
    returns the stacking sequence penalty for 10% rule

- ten_percent_rule
    returns only the stacking sequences that satisfy the 10% rule when added to
    plies for which the ply orientations have been previously determined

- calc_n_plies_per_angle
    returns the ply counts in each fibre direction
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np
import math

def calc_distance_2_points(point1, point2):
    """
    calculates the distance between two points in 2D

    Args:
        point1: coordinates of point 1
        point2: coordinates of point 2

    Returns:
        float: distance between point 1 and point 2

    Examples:
        >>> calc_distance_2_points([0, 0], [0, 2])
        2
    """
    return np.linalg.norm(point1 - point2, 2)


def calc_distance_Abdalla(LPs, constraints, num=10000):
    """
    calculates the distance between a lamination arameter point and the
    feasible lamination-parameter region for the 10% rule of Abdalla

    Args:
        LPs: lamination parameters of a laminate
        constraints (instance of the class Constraints): set of design
            guidelines
        num: number of points taken for respresenting the boundaries of the
            lamination-parameter feasible region

    Returns:
        float: distance between point (LP1, LP2) and the feasible
        lamination-parameter region for the 10% rule of Abdalla

    Examples:
        >>> constraints=Constraints(rule_10_percent=True, rule_10_Abdalla=True,
            percent_Abdalla=10)
        >>> calc_distance_Abdalla(LPs=np.array([0, 0]), constraints)
        0
    """

    if constraints.rule_10_percent and constraints.rule_10_Abdalla:
        if math.pow((1 - 4 * constraints.percent_Abdalla), 2) \
        + (1 - 4 * constraints.percent_Abdalla) * LPs[1] \
        - 2 * math.pow(LPs[0], 2) + 1e-15 > 0 \
        and 1 - 4 * constraints.percent_Abdalla - LPs[1] + 1e-15 > 0:
            return 0

        LP1_max = math.sqrt(1 - 4 * constraints.percent_Abdalla)
        LP1_min = - LP1_max

        def point_parabola_Abdalla(LP1, constraints):
            LP2 = 2 * (LP1 / (1 - 4 * constraints.percent_Abdalla)) **2 \
            - (1 - 4 * constraints.percent_Abdalla)
            return np.array([LP1, LP2])

        def point_straight_curve_Abdalla(LP1, constraints):
            LP2 = 1 - 4 * constraints.percent_Abdalla
            return np.array([LP1, LP2])


        min1 = min((calc_distance_2_points(
            LPs[0:2], point_parabola_Abdalla(LP1, constraints)) \
            for LP1 in np.linspace(LP1_min, LP1_max, num)))

        min2 = min((calc_distance_2_points(
            LPs[0:2], point_straight_curve_Abdalla(LP1, constraints)) \
            for LP1 in np.linspace(LP1_min, LP1_max, num)))

        return min(min1, min2)

#import sys
#sys.path.append(r'C:\BELLA')
#from src.BELLA.constraints import Constraints
#constraints =Constraints(rule_10_percent=True,
#                         rule_10_Abdalla=True,
#                         percent_Abdalla=10)
#print(calc_distance_Abdalla(
#    np.array([(1 - 4 * constraints.percent_Abdalla)/math.sqrt(2), 0]),
#    constraints))
