#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2020 Simon R Proud
#
# This file is part of PyCoxMunk.
#
# PyCoxMunk is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PyCoxMunk is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PyCoxMunk.  If not, see <http://www.gnu.org/licenses/>.
"""Utility functions used in multiple places in the code."""
import numpy as np


def check_and_reshape(arr, good_shape):
    """If array is single value then scale to match other arrays."""
    if arr.shape == (1,):
        arr = np.full(good_shape, arr[0])
    return arr


def check_type(in_val, var_typ):
    """Check that input variable is correct type.
    All inputs should be numpy array or a float."""
    if type(in_val) == float:
        return np.array([in_val])
    elif type(in_val) == np.ndarray:
        return in_val
    else:
        raise TypeError(f'{var_typ} must be a single float or numpy array! Got: {type(in_val)}')
