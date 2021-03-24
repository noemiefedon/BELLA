# -*- coding: utf-8 -*-
"""
This module contains functions for manipulating and combining Python arrays.
"""
__version__ = '1.0'
__author__ = 'Noemie Fedon'

import numpy as np

def max_arrays(array1, array2):
    """
    returns an array collating the element-wise maximum values from two arrays

    Args:
        array1 (numpy array): The first array
        array2 (numpy array): The second array

    Returns:
        array: the element-wise maximum values of ``array1`` and ``array2``.

    Raises:
        ValueError: If the size of the arrays ```array1`` and ``array2`` are
        different.

    Examples:
        >>> max_arrays(np.array([1., 4., 5.]), np.array([4., 3., 5.]))
        array([4., 4., 5.])
    """
    if isinstance(array1, (int, float)):
        array1 = np.array([array1])
    if isinstance(array2, (int, float)):
        array2 = np.array([array2])

    if array1.size == 1:
        if array2.size == 1:
            return np.array([max(array1[0], array2[0])], float)

        array1 = array1[0] * np.ones(array2.shape)
        return np.array([max(array1[ind], array2[ind]) \
                         for ind in range(array2.size)], float)

    if array2.size == 1:
        array2 = array2[0] * np.ones(array1.shape)
        return np.array([max(array1[ind], array2[ind]) \
                         for ind in range(array2.size)], float)

    if array1.size != array2.size:
        raise ValueError("Both arrays must have the same length.")

    return np.array([max(array1[ind], array2[ind]) \
                     for ind in range(array2.size)], float)
