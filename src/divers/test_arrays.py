# -*- coding: utf-8 -*-
"""
This module test the functions for manipulating and combining Python arrays.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import pytest
import numpy as np

from arrays import max_arrays

@pytest.mark.parametrize(
    "array1, array2, expect", [
        (np.array([1]), np.array([4, 3, 5]), np.array([4, 3, 5])),
        (np.array([1, 2, 3]), np.array([4, 3, 5]), np.array([4, 3, 5]))
        ])

def test_max_arrays(array1, array2, expect):
    output = max_arrays(array1, array2)
    assert (output == expect).all()

def test_max_arrays_error():
    array1 = np.array([1, 2])
    array2 = np.array([1, 2, 3])
    with pytest.raises(ValueError):
        max_arrays(array1, array2)