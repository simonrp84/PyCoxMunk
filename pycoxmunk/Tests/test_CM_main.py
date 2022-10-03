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
from pyresample import create_area_def
from datetime import datetime
from pycoxmunk.CM_Main import PyCoxMunk
from unittest import mock
from satpy import Scene
import dask.array as da
import xarray as xr
import numpy as np
import unittest


class TestCMMain(unittest.TestCase):
    """Tests for the overarching CM_Main module."""

    def setup_method(self, cmca):
        # configure self.attribute
        pass

    @staticmethod
    def _get_lonlats():
        return np.array([10, 11]), np.array([5, 6])

    def test_init(self):
        """Test creation of the CM class."""

        # Create test Scene
        scn = Scene()
        good_bnd_names = ['IR_108']
        good_angle_names = {'sza': 'solar_zenith_angle',
                            'vza': 'satellite_zenith_angle',
                            'saa': 'solar_azimuth_angle',
                            'vaa': 'satellite_azimuth_angle'}

        targ_area = create_area_def("source_area",
                                    "EPSG:4326",
                                    area_extent=(-50, -50, 50, 50),
                                    width=10,
                                    height=10)

        for band in good_bnd_names:
            scn[band] = xr.DataArray(da.from_array(np.zeros((10, 10))),
                                     coords={'y': np.zeros(10), 'x': np.zeros(10)},
                                     attrs={'start_time': datetime.utcnow(),
                                            'area': targ_area})

        for band in good_angle_names.keys():
            scn[good_angle_names[band]] = xr.DataArray(da.from_array(np.zeros((10, 10))),
                                                       coords={'y': np.zeros(10), 'x': np.zeros(10)},
                                                       attrs={'start_time': datetime.utcnow()})

        # Check we get an error if we don't pass a bool to various options.

        # Bad BRDF
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, do_brdf='potato')

        # Bad masking
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, mask_bad='tomato')

        # Bad deletion
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, delete_when_done='parsnip')

        # Bad Scene
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None)

        # Bad band names
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, band_names='potato')
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, band_names=[])

        # Bad angle names
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, good_bnd_names, angle_names=1.0)
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, good_bnd_names, angle_names={})
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, good_bnd_names, angle_names={'sza': 1.})
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, good_bnd_names, angle_names={'sza': 1., 'vza': 1., })
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, good_bnd_names, angle_names={'sza': 1., 'vza': 1., 'saa': 1., })
        # Good angle names
        PyCoxMunk(scn, good_bnd_names, angle_names=good_angle_names)

        # cm_calcangles
        mocker = mock.MagicMock()
        mocker_lalo = mock.MagicMock()
        mocker_lalo.attrs['area'].get_lonlats = self._get_lonlats
        tmp_dict = {'IR_108': mocker_lalo,
                    'solar_zenith_angle': np.array([10, 10]),
                    'satellite_zenith_angle': np.array([10, 10]),
                    'solar_azimuth_angle': np.array([10, 10]),
                    'satellite_azimuth_angle': np.array([10, 10])}
        mocker.return_value = tmp_dict
        with mock.patch('pycoxmunk.CM_Main.cm_calcangles', mocker):
            PyCoxMunk(scn, good_bnd_names, angle_names='calc')
            mocker.assert_called()

        # Bad ocean color dir
        with self.assertRaises(ValueError):
            PyCoxMunk(scn, good_bnd_names, angle_names=good_angle_names, oc_dir=1.0)

        # No ocean color dir
        pcm = PyCoxMunk(scn, good_bnd_names, angle_names=good_angle_names, oc_dir=None)
        self.assertFalse(pcm.use_occci)
