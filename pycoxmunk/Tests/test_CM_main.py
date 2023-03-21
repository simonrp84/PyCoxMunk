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
import pytest
from pyresample import create_area_def
from datetime import datetime
from pycoxmunk.CM_Main import PyCoxMunk
from unittest import mock
from satpy import Scene
import dask.array as da
import xarray as xr
import numpy as np


class TestCMMain:
    """Tests for the overarching CM_Main module."""

    @staticmethod
    def create_test_scene(band_names, angle_names):

        # Define orbital parameters used for lat/lon calculation
        orb_par = {'projection_longitude': 0.0,
                   'projection_latitude': 0.0,
                   'projection_altitude': 35784012.5,
                   'satellite_nominal_latitude': 0.0,
                   'satellite_nominal_longitude': 0.0,
                   'satellite_nominal_altitude': 36000000.0}
        # Create test Scene
        scn = Scene()
        targ_area = create_area_def("source_area",
                                    "EPSG:4326",
                                    area_extent=(-50, -50, 50, 50),
                                    width=10,
                                    height=10)
        for band in band_names:
            wvl_mock = mock.MagicMock()
            wvl_mock.central = band
            scn[band] = xr.DataArray(da.from_array(np.zeros((10, 10))),
                                     coords={'y': np.zeros(10), 'x': np.zeros(10)},
                                     attrs={'start_time': datetime.utcnow(),
                                            'area': targ_area,
                                            'wavelength': wvl_mock,
                                            'orbital_parameters': orb_par})

        for band in angle_names.keys():
            scn[angle_names[band]] = xr.DataArray(da.from_array(np.zeros((10, 10))),
                                                  coords={'y': np.zeros(10), 'x': np.zeros(10)},
                                                  attrs={'start_time': datetime.utcnow()})
        return scn

    def setup_method(self):
        self.good_bnd_names = ['VIS006']
        self.missing_angle_names = {'sza': 'solar_zenith_angle',
                                    'vza': 'satellite_zenith_angle',
                                    'saa': 'solar_azimuth_angle'}

        self.good_angle_names = self.missing_angle_names.copy()
        self.good_angle_names['vaa'] = 'satellite_azimuth_angle'

        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        self.scn_missing = self.create_test_scene(self.good_bnd_names, self.missing_angle_names)

    @staticmethod
    def _get_lonlats():
        return np.array([10, 11]), np.array([5, 6])

    def test_badbools(self):
        """Test some of the boolean options"""
        # Bad BRDF
        with pytest.raises(ValueError):
            PyCoxMunk(None, None, do_brdf='potato')

        # Bad masking
        with pytest.raises(ValueError):
            PyCoxMunk(None, None, mask_bad='tomato')

        # Bad deletion
        with pytest.raises(ValueError):
            PyCoxMunk(None, None, delete_when_done='parsnip')

    def test_badbands(self):
        """Test errors are raised if bad band / angle combinations are used."""
        # Bad band names
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, band_names='potato')
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, band_names=[])
        with pytest.raises(KeyError):
            PyCoxMunk(self.scn_good, band_names=['VIS008'])

        # Bad angle names
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names=1.0)
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names={})
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names={'sza': 1.})
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names={'sza': 1., 'vza': 1., })
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names={'sza': 1., 'vza': 1., 'saa': 1., })
        # Missing dataset in Scene
        with pytest.raises(KeyError):
            PyCoxMunk(self.scn_missing, self.good_bnd_names, angle_names=self.good_angle_names)

        # Good angle names
        PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names=self.good_angle_names)

    def test_ocean_color_init(self):
        """Test that ocean color options are parsed correctly."""
        # Bad ocean color dir
        with pytest.raises(ValueError):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names=self.good_angle_names, oc_dir=1.0)

        # No ocean color dir
        pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names=self.good_angle_names, oc_dir=None)
        assert not pcm.use_occci

        # Good ocean color dir
        ocdir = '/media/test_oc_dir/'
        pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names=self.good_angle_names, oc_dir=ocdir)
        assert pcm.use_occci
        assert pcm.oc_dir == ocdir

    def test_calcangles(self):
        """Check that angle calculation is initialised correctly."""
        mocker = mock.MagicMock()
        mocker_lalo = mock.MagicMock()
        mocker_lalo.attrs['area'].get_lonlats = self._get_lonlats
        tmp_dict = {'VIS006': mocker_lalo,
                    'solar_zenith_angle': np.array([10, 10]),
                    'satellite_zenith_angle': np.array([10, 10]),
                    'solar_azimuth_angle': np.array([10, 10]),
                    'satellite_azimuth_angle': np.array([10, 10])}
        mocker.return_value = tmp_dict
        with mock.patch('pycoxmunk.CM_Main.cm_calcangles', mocker):
            PyCoxMunk(self.scn_good, self.good_bnd_names, angle_names='calc')
            mocker.assert_called()

    def test_misc(self):
        """Misc tests for creation of the CM class."""
        # Bad Scene
        with pytest.raises(ValueError):
            PyCoxMunk(None, None)

    def test_wind_setup(self):
        """Test that winds are set up correctly."""
        mocker = mock.MagicMock()
        mocker.return_value = True
        with mock.patch('pycoxmunk.CM_Main.CMSharedWind', mocker):
            with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
                pcm = PyCoxMunk(self.scn_good, self.good_bnd_names)
            pcm.setup_wind(13.4, -12.3)
            assert pcm.shared_wind

    def _make_mocked_cm_refl(self):
        cmr_mocker = mock.MagicMock()

        self.tmp_cmr_array = np.random.uniform(0, 10, (10, 10))
        self.tmp_cmr_array[0, 0] = -10.

        cmr_mocker.rho = self.tmp_cmr_array.copy()
        cmr_mocker.rho_0d = self.tmp_cmr_array.copy() * 2.
        cmr_mocker.rho_0v = self.tmp_cmr_array.copy() * 4.
        cmr_mocker.rho_dv = self.tmp_cmr_array.copy() * 6.
        cmr_mocker.rho_dd = self.tmp_cmr_array.copy() * 8.

        cmr_mocker.rhoul = True
        cmr_mocker.rhogl = True
        cmr_mocker.rhowc = True

        return cmr_mocker

    def test_pixmask(self):
        """Test that pixmask is correctly initialised in the class."""
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names)
        assert pcm.pixmask is None
        pcm.setup_pixmask(np.array([0, 0, 1]), np.array([2., 0, 0]), np.array([0., 0, 0]), np.array([1., 0, 1.]))
        np.testing.assert_allclose(pcm.pixmask.mask, np.array([1, 0, 1]))

    @mock.patch('pycoxmunk.CM_Main.calc_coxmunk_wrapper')
    def test_retr_cm(self, mock_cmr_func):
        """Tests for the main routine that retrieves reflectance."""

        cm_band_name = f'cox_munk_refl_{self.good_bnd_names[0]}'
        cm_rho0v_name = f'cox_munk_rho0v_{self.good_bnd_names[0]}'

        # This section tests that cox_munk reflectance is computed, but not BRDF
        # Recreate scene in case previous tests have messed with it.
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=False, do_brdf=False)
        pcm.shared_wind = ['testval']
        pcm.retr_coxmunk_refl()
        np.testing.assert_allclose(pcm.scn[cm_band_name].data, self.tmp_cmr_array)
        assert cm_rho0v_name not in pcm.scn
        with pytest.raises(AttributeError):
            pcm.cm_refl.rhoul

        # Now we test for case when we also want the BRDF
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=False, do_brdf=True)
        pcm.retr_coxmunk_refl()
        np.testing.assert_allclose(pcm.scn[cm_band_name].data, self.tmp_cmr_array)
        assert cm_rho0v_name in pcm.scn
        np.testing.assert_allclose(pcm.scn[cm_rho0v_name].data, self.tmp_cmr_array * 4)

        # Finally, test case where we mask bad values
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=True, do_brdf=True)
        pcm.retr_coxmunk_refl()
        new_arr = self.tmp_cmr_array.copy()
        new_arr[0, 0] = np.nan
        np.testing.assert_allclose(pcm.scn[cm_band_name].data, new_arr)
        assert cm_rho0v_name in pcm.scn
        np.testing.assert_allclose(pcm.scn[cm_rho0v_name].data, new_arr * 4)

        # Now set up secondary pixmask class
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        cmask = np.zeros((10, 10))
        cmask[5, 5] = 1

        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=True, do_brdf=True)
        pcm.setup_pixmask(cloud_mask=cmask)
        pcm.retr_coxmunk_refl()
        new_arr = self.tmp_cmr_array.copy()
        new_arr[0, 0] = np.nan
        new_arr[5, 5] = np.nan
        np.testing.assert_allclose(pcm.scn[cm_band_name].data, new_arr)

    @mock.patch('pycoxmunk.CM_Main.calc_coxmunk_wrapper')
    def test_deleter(self, mock_cmr_func):
        """Test that deletions are done correctly."""

        tmp_arr = np.array([1, 4., 64])

        # Test for case where we don't delete intermediate data
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=False, do_brdf=False, delete_when_done=False)
        pcm.shared_wind = tmp_arr
        pcm.geometry = 'geotest'
        pcm.cm_refl.rhogl = True
        pcm._run_delete()
        assert pcm.geometry == 'geotest'
        assert pcm.cm_refl.rhogl
        np.testing.assert_allclose(pcm.shared_wind, tmp_arr)

        # Test for case where we delete intermediate data
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=False, do_brdf=False, delete_when_done=True)
        pcm.shared_wind = tmp_arr
        pcm.geometry = 'geotest'
        pcm._run_delete()
        assert not hasattr(pcm, 'geometry')
        assert not hasattr(pcm.cm_refl, 'rhogl')
        assert not hasattr(pcm, 'shared_wind')

        # Test for case where we delete intermediate data with some missing
        self.scn_good = self.create_test_scene(self.good_bnd_names, self.good_angle_names)
        mock_cmr_func.return_value = self._make_mocked_cm_refl()
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            pcm = PyCoxMunk(self.scn_good, self.good_bnd_names, mask_bad=False, do_brdf=False, delete_when_done=True)
        del pcm.shared_wind
        del pcm.geometry
        del pcm.cm_refl.rhogl
        del pcm.cm_refl.rhoul
        del pcm.cm_refl.rhowc
        pcm._run_delete()
