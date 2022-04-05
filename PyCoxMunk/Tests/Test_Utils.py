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

from PyCoxMunk.src.CM_Utils import check_type
import numpy as np
import unittest


class TestUtils(unittest.TestCase):
    def test_typechecker(self):
        """Checks for array/value types."""
        check_type(1., 'test')
        check_type(np.array([1.]), 'test')
        check_type(np.array([1]), 'test')
        with self.assertRaises(TypeError):
            check_type([1.], 'test')
        with self.assertRaises(TypeError):
            check_type(1, 'test')
