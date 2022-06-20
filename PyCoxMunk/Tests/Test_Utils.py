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
"""Test the scene geometry module."""

from PyCoxMunk.src import CM_Utils
import numpy as np
import unittest


class TestUtils(unittest.TestCase):
    def test_typechecker(self):
        """Checks for array/value types."""
        CM_Utils.check_type(1., 'test')
        CM_Utils.check_type(np.array([1.]), 'test')
        CM_Utils.check_type(np.array([1]), 'test')
        with self.assertRaises(TypeError):
            CM_Utils.check_type([1.], 'test')
        with self.assertRaises(TypeError):
            CM_Utils.check_type(1, 'test')

    def test_check_and_reshape(self):
        """Test the array resizer"""

        test_arr = np.random.uniform(0, 1, 10)
        out_arr = CM_Utils.check_and_reshape(test_arr, (10,))

        self.assertEqual(out_arr.shape, test_arr.shape)
        with self.assertRaises(ValueError):
            CM_Utils.check_and_reshape(test_arr, 100)

        test_arr = np.random.uniform(0, 1, 1)
        out_arr = CM_Utils.check_and_reshape(test_arr, 10)
        self.assertEqual(out_arr.shape, (10,))


    def test_gauss_leg(self):
        """Test the Gauss-Legendre calculations."""

        x, w = CM_Utils.gauss_leg_quadx(5, 0, 10)
        expected = np.array([1.18463443, 2.39314335, 2.84444444, 2.39314335, 1.18463443])
        np.testing.assert_almost_equal(x, expected)

        x, w = CM_Utils.gauss_leg_quadx(50, -10, 10)
        expected = np.array([0.02908623, 0.06759799, 0.10590548, 0.14380823, 0.18115561, 0.21780243,
                             0.25360674, 0.28842994, 0.32213728, 0.35459836, 0.38568757, 0.41528463,
                             0.44327504, 0.46955051, 0.49400938, 0.51655703, 0.53710622, 0.55557745,
                             0.57189926, 0.5860085, 0.59785059, 0.60737971, 0.614559, 0.61936067,
                             0.62176617, 0.62176617, 0.61936067, 0.614559, 0.60737971, 0.59785059,
                             0.5860085, 0.57189926, 0.55557745, 0.53710622, 0.51655703, 0.49400938,
                             0.46955051, 0.44327504, 0.41528463, 0.38568757, 0.35459836, 0.32213728,
                             0.28842994, 0.25360674, 0.21780243, 0.18115561, 0.14380823, 0.10590548,
                             0.06759799, 0.02908623])
        np.testing.assert_almost_equal(x, expected)

        x, w = CM_Utils.gauss_leg_quadx(12, 1, 100)
        expected = np.array([2.33517915, 5.29349664, 7.92387726, 10.05678762, 11.55788056, 12.33277877,
                             12.33277877, 11.55788056, 10.05678762, 7.92387726, 5.29349664, 2.33517915])
        np.testing.assert_almost_equal(x, expected)

        with self.assertRaises(ValueError):
            CM_Utils.gauss_leg_quadx(0, -10, 10)

