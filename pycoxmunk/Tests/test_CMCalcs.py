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

from pycoxmunk.CM_Constants import CM_DATA_DICT, WaterData
from pycoxmunk.CM_Shared_Wind import CMSharedWind
from pycoxmunk.CM_SceneGeom import CMSceneGeom
import pycoxmunk.CM_Calcs as CMCalcs
from copy import deepcopy
from unittest import mock
import numpy as np


class TestCMCalcs:
    def setup_method(self):
        """Set up some common variables for the tests."""
        self.sza = np.array([58.655, 167.412, 3.173, 63.215])
        self.vza = np.array([57.733, 21.248, 57.777, 99.768])
        self.saa = np.array([285.129, 133.130, 57.125, 223.555])
        self.vaa = np.array([38.505, 185.793, 133.017, 341.330])
        self.raa = np.array([113.376, 52.663, 75.892, 117.775])
        self.lats = np.array([-5.247, -58.473, 17.666, 48.237])
        self.lons = np.array([-116.184, 73.047, 120.744, -121.522])

        self.refl = CMCalcs.CMReflectance(cwvl=0.84, rho=np.array([1., 10., 100.]))
        self.watprops = deepcopy(CM_DATA_DICT[0.47])
        self.watprops.chlabs = 0.0109369
        self.watprops.chlbsc = 0.01805958
        with pytest.warns(UserWarning, match="Some solar zenith values out of range. Clipping."):
            self.geom = CMSceneGeom(self.sza, self.saa, self.vza, self.vaa, self.lats, self.lons)
        self.wind = CMSharedWind(self.geom, -5., 1.2)

    def test_cmrefl(self):
        """Test the CM_Reflectance class."""
        with pytest.raises(ValueError):
            CMCalcs.CMReflectance(cwvl='1.3')
        with pytest.raises(ValueError):
            CMCalcs.CMReflectance(cwvl=1)
        with pytest.raises(ValueError):
            CMCalcs.CMReflectance(cwvl=[1., 2.])

        cmref = CMCalcs.CMReflectance(cwvl=1.3)
        assert 1.3 == cmref.cwvl
        assert cmref.rhowc is None

        cmref = CMCalcs.CMReflectance(cwvl=1.3, rho=1.54, rhowc=12.3)
        assert 1.54 == cmref.rho
        assert 12.3 == cmref.rhowc

    def test_compute_bands(self):
        """Test code that selects Cox-Munk bands."""
        blist = [0.22, 0.47, 0.65, 0.83, 1.4, 1.92, 2.56, 3.4, 4.1]
        ex_res = [(None, 0.47), (None, 0.47), (0.65, 0.87), (0.65, 0.87),
                  (1.375, 1.6), (1.6, 2.13), (2.13, 3.7), (2.13, 3.7), (3.7, None)]
        with pytest.warns(UserWarning,
                          match="Warning: Band wavelength 0.22 is less than PyCoxMunk minimum wavelength 0.47"):
            assert CMCalcs._compute_bands_to_use(blist[0]), ex_res[0]
        with pytest.warns(UserWarning,
                          match="Warning: Band wavelength 0.47 is less than PyCoxMunk minimum wavelength 0.47"):
            assert CMCalcs._compute_bands_to_use(blist[1]), ex_res[1]
        for i in range(2, len(blist) - 1):
            assert CMCalcs._compute_bands_to_use(blist[i]), ex_res[i]
        with pytest.warns(UserWarning,
                          match="Warning: Band wavelength 4.1 is greater than PyCoxMunk maximum wavelength 3.7"):
            assert CMCalcs._compute_bands_to_use(blist[-1]), ex_res[-1]
        with pytest.raises(ValueError):
            CMCalcs._compute_bands_to_use(0.18)
        with pytest.raises(ValueError):
            CMCalcs._compute_bands_to_use(5.18)

    def test_interp_frac(self):
        """Test interpolation of wavelengths."""
        w1 = 0.45
        w2 = 0.65
        cw = [0.45, 0.47, 0.53, 0.59, 0.64]
        omins = [(1.0, 0.0), (0.9, 0.1), (0.6, 0.4), (0.3, 0.7), (0.05, 0.95)]
        for i in range(0, len(cw)):
            retr = CMCalcs._get_interp_frac(w1, w2, cw[i])
            np.testing.assert_almost_equal(retr[0], omins[i][0])
            np.testing.assert_almost_equal(retr[1], omins[i][1])

        with pytest.raises(ValueError):
            CMCalcs._get_interp_frac(1, 2, 0.5)
        with pytest.raises(ValueError):
            CMCalcs._get_interp_frac(1, 2, 2.5)

    @mock.patch('pycoxmunk.CM_Calcs._compute_bands_to_use')
    @mock.patch('pycoxmunk.CM_Calcs._get_interp_frac')
    def test_compute_water_wvl(self, interp, band):
        """Test selection of wavelength used for water properties."""
        band_mock = mock.MagicMock()
        interp_mock = mock.MagicMock()

        band.return_value = band_mock
        interp.return_value = interp_mock
        # Check error is raised if bad interpolation method is requested
        with pytest.raises(ValueError):
            CMCalcs.compute_wavelength_specific_water_props(0.8, meth='random')

        # Check error is raised if no previous band is available
        with pytest.raises(ValueError):
            band.return_value = (None, 0.9)
            CMCalcs.compute_wavelength_specific_water_props(0.8, meth='prev')
        # Check error is raised if no next band is available
        with pytest.raises(ValueError):
            band.return_value = (0.9, None)
            CMCalcs.compute_wavelength_specific_water_props(0.8, meth='next')

        # Check warning is raised if switching from interp to prev/next method
        with pytest.warns(UserWarning):
            band.return_value = (0.47, None)
            CMCalcs.compute_wavelength_specific_water_props(0.8, meth='interp')
        with pytest.warns(UserWarning):
            band.return_value = (None, 1.375)
            CMCalcs.compute_wavelength_specific_water_props(0.8, meth='interp')

        # Check correct data is returned for prev method.
        ret_wvl = 0.47
        band.return_value = (ret_wvl, None)
        ret_data = CMCalcs.compute_wavelength_specific_water_props(0.8, meth='prev')
        assert type(ret_data) is WaterData
        assert ret_data.whitecap_refl, CM_DATA_DICT[ret_wvl].whitecap_refl
        assert ret_data.wavelength, CM_DATA_DICT[ret_wvl].wavelength

        # Check correct data is returned for next method.
        ret_wvl = 1.375
        band.return_value = (None, ret_wvl)
        ret_data = CMCalcs.compute_wavelength_specific_water_props(0.8, meth='next')
        assert type(ret_data) is WaterData
        assert ret_data.whitecap_refl, CM_DATA_DICT[ret_wvl].whitecap_refl
        assert ret_data.wavelength, CM_DATA_DICT[ret_wvl].wavelength

        # Check interp works correctly
        frac1 = 0.3
        frac2 = 0.7
        band1 = 0.65
        band2 = 0.87
        interp.return_value = (frac1, frac2)
        band.return_value = (band1, band2)
        ret_data = CMCalcs.compute_wavelength_specific_water_props(0.8, meth='interp')
        assert ret_data.wavelength == 0.8
        assert ret_data.base_abs == CM_DATA_DICT[band1].base_abs * frac1 + CM_DATA_DICT[band2].base_abs * frac2
        assert ret_data.whitecap_refl == (CM_DATA_DICT[band1].whitecap_refl * frac1
                                          + CM_DATA_DICT[band2].whitecap_refl * frac2)

    def test_oceancolor(self):
        """Test the ocean color calculations. Functionality currently very basic in the main code."""
        from copy import deepcopy

        # Test we raise error if OC data is supplied, as we don't support it yet.
        with pytest.raises(NotImplementedError):
            CMCalcs.run_oceancolor(None, 'a string')

        # Test that correct data is returned for chlorophyll.
        expected = (0.0109369, 0.01805958)
        in_data = deepcopy(CM_DATA_DICT[0.47])
        in_data.whitecap_refl = 50
        out_data = CMCalcs.run_oceancolor(in_data)
        np.testing.assert_almost_equal(out_data.chlabs, expected[0])
        np.testing.assert_almost_equal(out_data.chlbsc, expected[1])

    def test_get_abcd(self):
        """Test computation of the a, b, c and d parameters."""
        firstvals = np.array((np.deg2rad(20), np.deg2rad(0.5), np.deg2rad(-20), np.deg2rad(88)))
        secondvals = np.array((np.deg2rad(1), np.deg2rad(120), np.deg2rad(15), np.deg2rad(-90)))

        expected_a = np.array([0.32556815, -0.8703557, -0.57357644, 0.0348995])
        expected_b = np.array([0.35836795, 0.86162916, -0.08715574, -0.0348995])
        expected_c = np.array([0.34432761, 1.76749402, -0.70020754, -0.03492077])
        expected_d = np.array([0.38386404, -1.69766312, -0.08748866, -0.03492077])

        result = CMCalcs._compute_abcd(firstvals, secondvals)

        np.testing.assert_almost_equal(result[0], expected_a)
        np.testing.assert_almost_equal(result[1], expected_b)
        np.testing.assert_almost_equal(result[2], expected_c)
        np.testing.assert_almost_equal(result[3], expected_d)

    def test_calc_coxmunk_brdf(self):
        """Test calculation of the Cox-Munk BRDF parameters."""
        mock_wind = mock.MagicMock()
        mock_wind.u10 = np.array([10.])
        mock_wind.v10 = np.array([-2.1])
        geom = CMSceneGeom(10., np.array([[165.2, 123.1, 22.6], [34, 21.2, 170.4], [0.45, 128.89, 44.22]]),
                           58.23, 9.3, 12., -120., )

        exp_dv = np.array([[1.702098, 1.673689, 1.715214],
                           [1.706389, 1.715966, 1.707012],
                           [1.714822, 1.675567, 1.696584]])
        exp_0d = np.array([[2.389888, 2.282412, 2.424333],
                           [2.401838, 2.426135, 2.403517],
                           [2.423389, 2.29209, 2.373374]])
        exp_dd = np.array([[1.97325, 1.945448, 1.980272],
                           [1.975696, 1.980641, 1.976037],
                           [1.98008, 1.948862, 1.969804]])
        res = CMCalcs.calc_cox_munk_brdf_terms(self.refl, 0.8, geom, mock_wind, None)
        np.testing.assert_allclose(self.refl.rho, res.rho_0v)
        np.testing.assert_allclose(exp_dd, res.rho_dd, atol=1e-6)
        np.testing.assert_allclose(exp_0d, res.rho_0d, atol=1e-6)
        np.testing.assert_allclose(exp_dv, res.rho_dv, atol=1e-6)

    @mock.patch('pycoxmunk.CM_Calcs.compute_wavelength_specific_water_props')
    @mock.patch('pycoxmunk.CM_Calcs.run_oceancolor')
    @mock.patch('pycoxmunk.CM_Calcs._compute_abcd')
    def test_calc_coxmunk(self, abcd, run_oc, watprops):
        """Test the main Cox-Munk calculation."""

        watprops.return_value = self.watprops
        run_oc.return_value = self.watprops
        abcd.return_value = (0.32556815, 0.35836795, 0.34432761, 0.38386404)

        exp_rho = np.array([0.8256082, 0.8417819, 0.8079787, 0.8281935])
        exp_wc = np.array([0.0414373, 0.0414373, 0.0414373, 0.0414373])
        exp_gl = np.array([1.55336676e-10, 0.0, 2.42690648e-05, 0.0])
        exp_ul = np.array([0.78490883, 0.80109773, 0.76723845, 0.78749654])

        ref_data = CMCalcs.calc_cox_munk(0.47, self.geom, self.wind)

        np.testing.assert_almost_equal(exp_rho, ref_data.rho)
        np.testing.assert_almost_equal(exp_wc, ref_data.rhowc)
        np.testing.assert_almost_equal(exp_gl, ref_data.rhogl)
        np.testing.assert_almost_equal(exp_ul, ref_data.rhoul)

    @mock.patch('pycoxmunk.CM_Calcs.calc_cox_munk')
    @mock.patch('pycoxmunk.CM_Calcs.calc_cox_munk_brdf_terms')
    def test_cm_wrapper(self, ccm_brdf, ccm):
        """Test the wrapper function to ensure it makes calls correctly."""

        # Case with no BRDF
        CMCalcs.calc_coxmunk_wrapper(None, None, None, None, False)
        assert ccm.called
        assert not ccm_brdf.called
        # Case with BRDF
        CMCalcs.calc_coxmunk_wrapper(None, None, None, None, True)
        assert ccm.called
        assert ccm_brdf.called
