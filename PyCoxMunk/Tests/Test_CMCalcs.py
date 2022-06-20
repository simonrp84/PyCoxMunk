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

import PyCoxMunk.src.CM_Calcs as CMCalcs
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

    def test_cmrefl(self):
        """Test the CM_Reflectance class."""
        with self.assertRaises(ValueError):
            CMCalcs.CM_Reflectance(cwvl='1.3')
        with self.assertRaises(ValueError):
            CMCalcs.CM_Reflectance(cwvl=1)
        with self.assertRaises(ValueError):
            CMCalcs.CM_Reflectance(cwvl=[1., 2.])

        cmref = CMCalcs.CM_Reflectance(cwvl=1.3)
        self.assertEqual(1.3, cmref.cwvl)
        self.assertEqual(None, cmref.rhowc)

        cmref = CMCalcs.CM_Reflectance(cwvl=1.3, rho=1.54, rhowc=12.3)
        self.assertEqual(1.54, cmref.rho)
        self.assertEqual(12.3, cmref.rhowc)

    def test_compute_bands(self):
        """Test code that selects Cox-Munk bands."""
        blist = [0.22, 0.47, 0.65, 0.83, 1.4, 1.92, 2.56, 3.4]
        ex_res = [(None, 0.47), (None, 0.47), (0.65, 0.87), (0.65, 0.87),
                  (1.375, 1.6), (1.6, 2.13), (2.13, 3.7), (2.13, 3.7)]
        for i in range(0, len(blist)):
            self.assertEqual(CMCalcs._compute_bands_to_use(blist[i]), ex_res[i])
        with self.assertRaises(ValueError):
            CMCalcs._compute_bands_to_use(0.18)
        with self.assertRaises(ValueError):
            CMCalcs._compute_bands_to_use(5.18)

    def test_interp_frac(self):
        """Test interpolation of wavelengths."""
        w1 = 0.45
        w2 = 0.65
        cw = [0.45, 0.47, 0.53, 0.59, 0.64]
        omins = [(1.0, 0.0), (0.9, 0.1), (0.6, 0.4), (0.3, 0.7), (0.05, 0.95)]
        for i in range(0, len(cw)):
            retr = CMCalcs._get_interp_frac(w1, w2, cw[i])
            self.assertAlmostEqual(retr[0], omins[i][0])
            self.assertAlmostEqual(retr[1], omins[i][1])

        with self.assertRaises(ValueError):
            CMCalcs._get_interp_frac(1, 2, 0.5)
        with self.assertRaises(ValueError):
            CMCalcs._get_interp_frac(1, 2, 2.5)


if __name__ == '__main__':
    unittest.main()
