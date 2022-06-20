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

from PyCoxMunk.CM_Shared_Wind import CMSharedWind
from PyCoxMunk.CM_SceneGeom import CMSceneGeom
import numpy as np
import unittest


class TestSharedWind(unittest.TestCase):
    def setUp(self):
        """Set up some common variables for the tests."""
        self.u10 = np.array([15., -4., -1, 10.7])
        self.v10 = np.array([8.1, 4.0, 11.9, -8.1])

        self.ws = np.array([17.04728717, 5.65685425, 11.94194289, 13.42013413])
        self.wd = np.array([2.38240853, 4.74502664, 5.37000261, 4.60018284])
        self.wcfrac = np.array([0.06388454, 0.00131533, 0.01825043, 0.02752158])
        self.ergodic = np.array([0.53386542, 0.93202052, -9.53116332, -9.9822338])
        self.cosomega = np.array([-0.00881482, 0.30019419,  0.10046262, -0.37514133])
        self.cosbeta = np.array([0.74863662, 0.63201908, 0.73177281, 0.48107697])
        self.w = np.array([0.78980563, 0.63295005, 0.73508197, 0.97767278])
        self.p = np.array([2.08958099e-03, 6.53917194e-23, 2.72598232e-05, 2.21290767e-19])

        self.sza = np.array([58.655, 167.412, 3.173, 63.215])
        self.vza = np.array([57.733, 21.248, 87.777, 99.768])
        self.saa = np.array([285.129, 133.130, 57.125, 223.555])
        self.vaa = np.array([38.505, 185.793, 133.017, 341.330])
        self.raa = np.array([113.376, 52.663, 75.892, 117.775])
        self.lats = np.array([-5.247, -58.473, 17.666, 48.237])
        self.lons = np.array([-116.184, 73.047, 120.744, -121.522])

        self.scene_geom = CMSceneGeom(self.sza, self.saa, self.vza, self.vaa, self.lats, self.lons, raa=self.raa)

    def test_sharedwind(self):
        """Test the calculation of relative azimuth angles."""
        sh_wind = CMSharedWind(self.scene_geom, self.u10, self.v10)
        np.testing.assert_allclose(sh_wind.ws, self.ws, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.wd, self.wd, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.wcfrac, self.wcfrac, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.ergodic, self.ergodic, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.cosomega, self.cosomega, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.cosbeta, self.cosbeta, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.w, self.w, rtol=1e-5)
        np.testing.assert_allclose(sh_wind.p, self.p, rtol=1e-5)


if __name__ == '__main__':
    unittest.main()
