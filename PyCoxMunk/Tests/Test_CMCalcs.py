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
"""Test the Cox-Munk calculations module."""

from PyCoxMunk.src.CM_SceneGeom import CMSceneGeom as Cm_sg
from PyCoxMunk.src.CM_SceneGeom import cm_calcangles
from unittest import mock
import numpy as np
import unittest


class TestCMCalcs(unittest.TestCase):
    def setUp(self):
        """Set up some common variables for the tests."""
        self.sza = np.array([58.655, 167.412, 3.173, 63.215])
        self.vza = np.array([57.733, 21.248, 57.777, 99.768])
        self.saa = np.array([285.129, 133.130, 57.125, 223.555])
        self.vaa = np.array([38.505, 185.793, 133.017, 341.330])
        self.raa = np.array([113.376, 52.663, 75.892, 117.775])
        self.lats = np.array([-5.247, -58.473, 17.666, 48.237])
        self.lons = np.array([-116.184, 73.047, 120.744, -121.522])


if __name__ == '__main__':
    unittest.main()
