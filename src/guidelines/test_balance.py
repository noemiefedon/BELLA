#  - * -  coding: utf - 8  - * -
"""
This module tests the functions in balance.py.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import pytest
import numpy as np

sys.path.append(r'C:\BELLA')
from src.LAYLA_V02.constraints import Constraints
from src.guidelines.balance import is_balanced
from src.guidelines.balance import calc_penalty_bal

@pytest.mark.parametrize(
    "stack, constraints, expect", [
        (np.array([0, 45, 90]), Constraints(), False),
        (np.array([0, 45, 90, -45]), Constraints(), True)
        ])

def test_is_balanced(stack, constraints, expect):
    output = is_balanced(stack, constraints)
    assert output == expect

@pytest.mark.parametrize(
    "n_plies_per_angle, constraints, cummul_areas, expect", [
        (np.array([0, 0, 1, 0]), Constraints(), 1, 1.),
        (np.array([0, 0, 2, 0]), Constraints(), 0.5, 0.5),
        (np.array([[0, 0, 2, 0],
                   [0, 2, 2, 0]]), Constraints(), 1, np.array([1., 0.5]))
        ])

def test_calc_penalty_bal(
        n_plies_per_angle, constraints, cummul_areas, expect):
    output = calc_penalty_bal(n_plies_per_angle, constraints, cummul_areas)
    assert (output == expect).all()
