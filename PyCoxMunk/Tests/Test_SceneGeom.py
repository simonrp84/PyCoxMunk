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

from PyCoxMunk.Utils.CM_SceneGeom import CMSceneGeom as Cm_sg
from mock import patch
import numpy as np
import unittest


class TestWaterClass(unittest.TestCase):
    def setUp(self):
        """Set up some common variables for the tests."""
        self.sza = np.array([58.655, 167.412, 3.173, 63.215])
        self.vza = np.array([57.733, 21.248, 57.777, 99.768])
        self.saa = np.array([285.129, 133.130,  57.125, 223.555])
        self.vaa = np.array([38.505, 185.793, 133.017, 341.330])
        self.raa = np.array([113.376, 52.663,  75.892, 117.775])
        self.lats = np.array([-5.247, -58.473,  17.666, 48.237])
        self.lons = np.array([-116.184, 73.047, 120.744, -121.522])

    def test_raa_calc(self):
        """Test the calculation of relative azimuth angles."""

        # Check the calculation function is called when no RAA is supplied
        with patch('PyCoxMunk.Utils.CM_SceneGeom.CMSceneGeom.calc_relazi') as calc_func:
            calc_func.side_effect = np.array([1.])
            Cm_sg(1., 1., 1., 1., 1., 1.)
            calc_func.assert_called_once()

        # Check the calculation function not called if RAA supplied
        with patch('PyCoxMunk.Utils.CM_SceneGeom.CMSceneGeom.calc_relazi') as calc_func:
            Cm_sg(1., 1., 1., 1., 1., 1., 1.)
            calc_func.assert_not_called()

        # Check returned values are correct in case of single values
        geom = Cm_sg(1., 50., 1., 30., 1., 1.)
        self.assertAlmostEqual(geom.raa, 20.)

        # And likewise for arrays
        geom = Cm_sg(self.sza, self.saa, self.vza,
                                 self.vaa, self.lats, self.lons)
        np.testing.assert_allclose(geom.raa, self.raa)

    def test_bounds_checker(self):
        """Ensure that bounds checks catch bad input data"""
        # SZA
        with self.assertRaises(ValueError):
            Cm_sg(Cm_sg.zenith_min-5, 1., 1., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(Cm_sg.zenith_max+100., 1., 1., 1., 1., 1.)
        Cm_sg(Cm_sg.zenith_max-10., 1., 1., 1., 1., 1.)

        # SAA
        with self.assertRaises(ValueError):
            Cm_sg(1., Cm_sg.azimuth_min-20., 1., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., Cm_sg.azimuth_max+20., 1., 1., 1., 1.)
        Cm_sg(1., Cm_sg.azimuth_max-20., 1., 1., 1., 1.)

        # VZA
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., Cm_sg.zenith_min-15., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., Cm_sg.zenith_max+15., 1., 1., 1.)
        Cm_sg(1., 1., Cm_sg.zenith_max-10., 1., 1., 1.)

        # VAA
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., Cm_sg.azimuth_min-25., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., Cm_sg.azimuth_max+50., 1., 1.)
        Cm_sg(1., 1., 1., Cm_sg.azimuth_max-1., 1., 1.)

        # Lats
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., Cm_sg.lat_min-25., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., Cm_sg.lat_max+25., 1.)
        Cm_sg(1., 1., 1., 1., Cm_sg.lat_max-15., 1.)

        # Lons
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., Cm_sg.lon_min-0.01)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., Cm_sg.lon_max+10.)
        Cm_sg(1., 1., 1., 1., 1., Cm_sg.lon_min+0.01)

        # RAA
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., 1., Cm_sg.relazi_min-10.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., 1., Cm_sg.relazi_max+11.)
        Cm_sg(1., 1., 1., 1., 1., 1., Cm_sg.relazi_min+10.)

        # Also do one test for the case where we pass arrays
        with self.assertRaises(ValueError):
            tmp_arr = np.array([1., 1.])
            Cm_sg(np.array([Cm_sg.zenith_min-5., 1.]), tmp_arr, tmp_arr,
                            tmp_arr, tmp_arr, tmp_arr)
        Cm_sg(tmp_arr, tmp_arr, tmp_arr,
              tmp_arr, tmp_arr, tmp_arr)




if __name__ == '__main__':
    unittest.main()
