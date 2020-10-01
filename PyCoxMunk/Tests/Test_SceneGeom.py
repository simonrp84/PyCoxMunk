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

from PyCoxMunk.Utils import CM_SceneGeom as CM_SG
import numpy as np
import unittest


class TestWaterClass(unittest.TestCase):
    def SetUp(self):
        """Set up some common variables for the tests."""
        self.sza = np.array([12.0, 18.9, 85.0, 92.0, 170.0])

    def Test_RAA_Calc(self):
        """Test the calculation of relative azimuth angles."""


if __name__ == '__main__':
    unittest.main()
