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

from PyCoxMunk.CM_SceneGeom import CMSceneGeom as Cm_sg
from PyCoxMunk.CM_SceneGeom import cm_calcangles
from unittest import mock
import numpy as np
import unittest


class TestSceneGeom(unittest.TestCase):
    def setUp(self):
        """Set up some common variables for the tests."""
        self.sza = np.array([58.655, 167.412, 3.173, 63.215])
        self.vza = np.array([57.733, 21.248, 57.777, 99.768])
        self.saa = np.array([285.129, 133.130, 57.125, 223.555])
        self.vaa = np.array([38.505, 185.793, 133.017, 341.330])
        self.raa = np.array([113.376, 52.663, 75.892, 117.775])
        self.lats = np.array([-5.247, -58.473, 17.666, 48.237])
        self.lons = np.array([-116.184, 73.047, 120.744, -121.522])

    def test_raa_calc(self):
        """Test the calculation of relative azimuth angles."""

        # Check the calculation function is called when no RAA is supplied
        with mock.patch('PyCoxMunk.CM_SceneGeom.CMSceneGeom.calc_relazi') as calc_func:
            calc_func.side_effect = np.array([1.])
            Cm_sg(1., 1., 1., 1., 1., 1.)
            calc_func.assert_called_once()

        # Check the calculation function not called if RAA supplied
        with mock.patch('PyCoxMunk.CM_SceneGeom.CMSceneGeom.calc_relazi') as calc_func:
            Cm_sg(1., 1., 1., 1., 1., 1., 1.)
            calc_func.assert_not_called()

        # Check returned values are correct in case of single values
        geom = Cm_sg(1., 50., 1., 30., 1., 1.)
        self.assertAlmostEqual(geom.raa, 20.)

        # And likewise for arrays
        geom = Cm_sg(self.sza, self.saa, self.vza,
                     self.vaa, self.lats, self.lons)
        np.testing.assert_allclose(geom.raa, self.raa)

    @mock.patch('satpy.modifiers.angles._get_sensor_angles')
    @mock.patch('satpy.modifiers.angles._get_sun_angles')
    def test_calcangles(self, get_sun, get_sat):
        """Test user-requested angle calculation."""
        from satpy import Scene
        import xarray as xr
        get_sat.return_value = (self.vaa, self.vza)
        get_sun.return_value = (self.saa, self.sza)

        test_scn = Scene()
        test_scn['VIS006'] = xr.DataArray(np.array([5., 4., 4., 5.]))

        retr = cm_calcangles(test_scn, 'VIS006')
        np.testing.assert_allclose(retr['solar_zenith_angle'].data, self.sza)
        np.testing.assert_allclose(retr['solar_azimuth_angle'].data, self.saa)
        np.testing.assert_allclose(retr['satellite_zenith_angle'].data, self.vza)
        np.testing.assert_allclose(retr['satellite_azimuth_angle'].data, self.vaa)

        # Check what happens if bad refband is passed
        with self.assertRaises(KeyError):
            cm_calcangles(test_scn, 'VIS008')

        # Check if we pass bad data
        with self.assertRaises(KeyError):
            test_scn['VIS006'] = xr.DataArray(np.array([[5., 4.], [4., 5.]]))
            cm_calcangles(test_scn, 'VIS008')

    def test_shapechecker(self):
        """Tests for the shape checker that ensures arrays are same size."""
        # Default case
        Cm_sg(self.sza, self.saa, self.vza,
              self.vaa, self.lats, self.lons)

        # Mismatched arrays
        with self.assertRaises(ValueError):
            Cm_sg(self.sza, self.saa, self.vza,
                  self.vaa, np.array([1., 2.]), self.lons)

        # One array, rest scalar
        retr = Cm_sg(self.sza, 5., 3., -10., 15., -21.)
        self.assertEqual(len(retr.sza), 4)
        self.assertEqual(len(retr.saa), 4)

    def test_bounds_checker(self):
        """Ensure that bounds checks catch bad input data"""
        # SZA
        test_cmsg = Cm_sg(0., 0., 0., 0., 0., 0.)
        Cm_sg(test_cmsg.zenith_min - 5, 1., 1., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(test_cmsg.zenith_min - 5, 1., 1., 1., 1., 1., fix_angs=False)

        Cm_sg(test_cmsg.zenith_max + 100., 1., 1., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(test_cmsg.zenith_max + 100., 1., 1., 1., 1., 1., fix_angs=False)
        Cm_sg(test_cmsg.zenith_max - 10., 1., 1., 1., 1., 1.)

        # SAA
        Cm_sg(1., test_cmsg.azimuth_min - 20., 1., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., test_cmsg.azimuth_min - 20., 1., 1., 1., 1., fix_angs=False)
        Cm_sg(1., test_cmsg.azimuth_min + 20., 1., 1., 1., 1.)

        Cm_sg(1., test_cmsg.azimuth_max + 20., 1., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., test_cmsg.azimuth_max + 20., 1., 1., 1., 1., fix_angs=False)
        Cm_sg(1., test_cmsg.azimuth_max - 20., 1., 1., 1., 1.)

        # VZA
        Cm_sg(1., 1., test_cmsg.zenith_min - 15., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., test_cmsg.zenith_min - 15., 1., 1., 1., fix_angs=False)
        Cm_sg(1., 1., test_cmsg.zenith_min + 10., 1., 1., 1.)

        Cm_sg(1., 1., test_cmsg.zenith_max + 15., 1., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., test_cmsg.zenith_max + 15., 1., 1., 1., fix_angs=False)
        Cm_sg(1., 1., test_cmsg.zenith_max - 10., 1., 1., 1.)

        # VAA
        Cm_sg(1., 1., 1., test_cmsg.azimuth_min - 25., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., test_cmsg.azimuth_min - 25., 1., 1., fix_angs=False)
        Cm_sg(1., 1., 1., test_cmsg.azimuth_min + 1., 1., 1.)

        Cm_sg(1., 1., 1., test_cmsg.azimuth_max + 50., 1., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., test_cmsg.azimuth_max + 50., 1., 1., fix_angs=False)
        Cm_sg(1., 1., 1., test_cmsg.azimuth_max - 1., 1., 1.)

        # Lats
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., test_cmsg.lat_min - 25., 1.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., test_cmsg.lat_max + 25., 1.)
        Cm_sg(1., 1., 1., 1., test_cmsg.lat_max - 15., 1.)

        # Lons
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., test_cmsg.lon_min - 0.01)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., test_cmsg.lon_max + 10.)
        Cm_sg(1., 1., 1., 1., 1., test_cmsg.lon_min + 0.01)

        # RAA
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., 1., check_raa=True, raa=test_cmsg.relazi_min - 10.)
        with self.assertRaises(ValueError):
            Cm_sg(1., 1., 1., 1., 1., 1., raa=test_cmsg.relazi_max + 11., fix_angs=False, check_raa=True)
        Cm_sg(1., 1., 1., 1., 1., 1., raa=test_cmsg.relazi_min + 10., fix_angs=False, check_raa=True)

        # Also do one test for the case where we pass arrays
        with self.assertRaises(ValueError):
            tmp_arr = np.array([1., 1.])
            Cm_sg(np.array([test_cmsg.zenith_min - 5., 1.]), tmp_arr, tmp_arr,
                  tmp_arr, tmp_arr, tmp_arr, fix_angs=False)
        Cm_sg(tmp_arr, tmp_arr, tmp_arr,
              tmp_arr, tmp_arr, tmp_arr)


if __name__ == '__main__':
    unittest.main()
