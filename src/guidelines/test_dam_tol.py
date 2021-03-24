#  - * -  coding: utf - 8  - * -
"""
This module tests the functions in dam_tol.py.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import sys
import pytest
import numpy as np

sys.path.append(r'C:\BELLA')
from src.LAYLA_V02.constraints import Constraints
from src.guidelines.dam_tol import is_dam_tol

@pytest.mark.parametrize(
    "stack, constraints, expect", [
        (np.array([45, 0, -45]), Constraints(dam_tol=True, dam_tol_rule=1), True),
        (np.array([45, 0, 0]), Constraints(dam_tol=True, dam_tol_rule=1), False),
        (np.array([45, -45, 0, 45, -45]), Constraints(dam_tol=True, dam_tol_rule=2), True),
        (np.array([45, -45, 0, 90, -45]), Constraints(dam_tol=True, dam_tol_rule=2), False),
        ])

def test_is_dam_tol(stack, constraints, expect):
    output = is_dam_tol(stack, constraints)
    assert output == expect
