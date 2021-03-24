#  - * -  coding: utf - 8  - * -
"""
This module test the functions for the 10% rule.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import numpy as np
import math as ma
import pytest

sys.path.append(r'C:\BELLA')
from src.BELLA.constraints import Constraints
from src.guidelines.ten_percent_rule import is_ten_percent_rule
from src.guidelines.ten_percent_rule_Abdalla import calc_distance_2_points
from src.guidelines.ten_percent_rule_Abdalla import calc_distance_Abdalla

@pytest.mark.parametrize(
    """constraints, stack, ply_queue, n_plies_per_angle, equality_45_135,
equality_0_90, expect""", [
        (Constraints(rule_10_percent=True, percent_0=50),
         np.array([0, 45, 90], int), [], None, False, False, False),
        (Constraints(rule_10_percent=True, percent_0=50),
         np.array([0, 45, 90], int), None, None, False, False, False),
        (Constraints(rule_10_percent=True, percent_0=50),
         np.array([0, 666, 666], int), [0, 45], None, False, False, True),
        (Constraints(rule_10_percent=True, percent_0=50),
         None, None, np.array([3, 0, 3, 0]), False, False, False),
        (Constraints(rule_10_percent=True, percent_0=50),
         None, None, np.array([0, 3, 0, 3]), False, False, True)
        ])

def test_is_ten_percent_rule(
        constraints, stack, ply_queue, n_plies_per_angle, equality_45_135,
        equality_0_90, expect):
    output = is_ten_percent_rule(
        constraints, stack, ply_queue, n_plies_per_angle, equality_45_135,
        equality_0_90)
    assert output == expect

@pytest.mark.parametrize(
    """point1, point2, expect""", [
        (np.array([0, 0]), np.array([0, 0]), 0),
        (np.array([0, 0]), np.array([0, 2]), 2),
        (np.array([0, 0]), np.array([1, 1]), ma.sqrt(2)),
        ])

def test_calc_distance_2_points(point1, point2, expect):
    output = calc_distance_2_points(point1, point2)
    assert output == expect


@pytest.mark.parametrize(
    """LPs, constraints, expect""", [
        (np.array([0, 0]), Constraints(rule_10_percent=True,
         rule_10_Abdalla=True, percent_Abdalla=10), 0),
         (np.array([0, 0.7]), Constraints(rule_10_percent=True,
         rule_10_Abdalla=True, percent_Abdalla=10), 0.1),
        ])

def test_calc_distance_Abdalla(LPs, constraints, expect):
    output = calc_distance_Abdalla(LPs, constraints)
    assert abs(output - expect) < 1e-5


