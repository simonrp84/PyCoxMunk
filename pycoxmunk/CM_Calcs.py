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
"""Class for the Scene wind and sea surface information."""

from pycoxmunk.CM_Constants import CM_DATA_DICT, WaterData, chlconc, n_air, dither_more, \
    cm_min_wvl, cm_max_wvl, n_quad_theta, n_quad_phi
from pycoxmunk.CM_Shared_Wind import CMSharedWind
from pycoxmunk.CM_SceneGeom import CMSceneGeom
from pycoxmunk.CM_Utils import gauss_leg_quadx
from copy import deepcopy
import dask.array as da
from numba import jit
import numpy as np

import warnings
import copy


class CMReflectance:
    def __init__(self, cwvl=None, rho=None, rhowc=None, rhogl=None, rhoul=None,
                 rho_0v=None, rho_0d=None, rho_dv=None, rho_dd=None):
        """Class for holding the output Cox-Munk reflectances and BRDF."""
        if type(cwvl) is float or cwvl is None:
            self.cwvl = cwvl  # The central wavelength of the band reflectances
        else:
            raise ValueError("Central wavelength must be a float!")
        self.rho = rho  # The reflectance
        self.rhowc = rhowc  # White cap reflectance
        self.rhogl = rhogl  # Glint reflectance
        self.rhoul = rhoul  # Underlight reflectance
        self.rho_0v = rho_0v  # The solar beam to satellite view reflectances
        self.rho_0d = rho_0d  # The solar beam to diffuse reflectances
        self.rho_dv = rho_dv  # The diffuse to satellite view reflectances
        self.rho_dd = rho_dd  # The diffuse to diffuse reflectances


def _compute_bands_to_use(cwvl):
    """Find which bands from the constants to use in processing.
    Inputs:
      - cwvl: Float, band central wavelength in micron.
    Returns:
      - bnames: Tuple of dict keys specifying previous/next wavelengths to use.
    """
    wvl_keys = list(CM_DATA_DICT.keys())
    wvl_keys.sort()
    if cwvl <= cm_min_wvl or cwvl > cm_max_wvl:
        raise ValueError(f"Central wavelength must be between {cm_min_wvl} and {cm_max_wvl}.")
    if cwvl <= wvl_keys[0]:
        warnings.warn(f"Warning: Band wavelength {cwvl} is less than PyCoxMunk minimum wavelength {wvl_keys[0]}")
        return None, wvl_keys[0]
    elif cwvl >= wvl_keys[-1]:
        warnings.warn(f"Warning: Band wavelength {cwvl} is greater than PyCoxMunk maximum wavelength {wvl_keys[-1]}")
        return wvl_keys[-1], None
    else:
        for i in range(0, len(wvl_keys) - 1):
            if wvl_keys[i] <= cwvl < wvl_keys[i + 1]:
                return wvl_keys[i], wvl_keys[i + 1]


def _get_interp_frac(prev_b, next_b, centre):
    """Compute relative contributions of two wavelengths to centre wavelength interpolation.
    Inputs:
      - prev_b: Float, previous wavelength to use.
      - next_b: Float, next wavelength to use.
      - centre: Float, requested central wavelength.
    Returns:
      - tuple, previous and next fractional contributions.
    """
    if centre < prev_b:
        raise ValueError("Central wavelength must be greater than previous wavelength.")
    if centre > next_b:
        raise ValueError("Central wavelength must be less than next wavelength.")

    w_diff = next_b - prev_b
    p_frac = 1 - (centre - prev_b) / w_diff
    n_frac = 1 - (next_b - centre) / w_diff

    return p_frac, n_frac


def compute_wavelength_specific_water_props(cwvl, meth='interp'):
    """Compute the water properties for a specific central wavelength using one of several methods.
    Inputs:
      - cwvl: Float, the wavelength of the desired band.
      - meth: Str, computation method. 'interp' for linear interpolation prev/next Cox-Munk bands,
             'prev' for using only the previous (lower wavelength) band,
             and 'next' for using only the next (higher wavelength) band.
             Default: 'interp'
    Returns:
      - WaterData, containing the required water information on the specified wavelength.
    """
    meths = ['interp', 'prev', 'next']
    if meth not in meths:
        raise ValueError("Wavelength selection method must be one of 'interp', 'prev' or 'next'.")

    # Compute which Cox-Munk bands we should use
    p_band, n_band = _compute_bands_to_use(cwvl)

    # Some sanity checks to prevent bad data
    if meth == 'prev' and p_band is None:
        raise ValueError("Cannot use previous wavelength as user-requested wavelength too low.")
    if meth == 'next' and n_band is None:
        raise ValueError("Cannot use next wavelength as user-requested wavelength too high.")
    if meth == 'interp' and p_band is None:
        warnings.warn("Cannot use wavelength interpolation as user-requested wavelength too low. Using next.",
                      UserWarning)
        meth = 'next'
    if meth == 'interp' and n_band is None:
        warnings.warn("Cannot use wavelength interpolation as user-requested wavelength too high. Using previous.",
                      UserWarning)
        meth = 'prev'

    # This is simple if we're only using previous / next wavelength
    if meth == 'prev':
        return deepcopy(CM_DATA_DICT[p_band])
    elif meth == 'next':
        return deepcopy(CM_DATA_DICT[n_band])

    # More complex calculations needed for interpolation. First, find fraction of each wavelength to use.
    p_frac, n_frac = _get_interp_frac(p_band, n_band, cwvl)

    new_wat_data = WaterData(wavelength=cwvl,
                             refrac_real=(p_frac * CM_DATA_DICT[p_band].refrac_real +
                                          n_frac * CM_DATA_DICT[n_band].refrac_real),
                             refrac_imag=(p_frac * CM_DATA_DICT[p_band].refrac_imag +
                                          n_frac * CM_DATA_DICT[n_band].refrac_imag),
                             base_abs=(p_frac * CM_DATA_DICT[p_band].base_abs +
                                       n_frac * CM_DATA_DICT[n_band].base_abs),
                             base_backscat=(p_frac * CM_DATA_DICT[p_band].base_backscat +
                                            n_frac * CM_DATA_DICT[n_band].base_backscat),
                             whitecap_refl=(p_frac * CM_DATA_DICT[p_band].whitecap_refl
                                            + n_frac * CM_DATA_DICT[n_band].whitecap_refl),
                             chl_a_coef=[p_frac * CM_DATA_DICT[p_band].chl_a_coef[0] +
                                         n_frac * CM_DATA_DICT[n_band].chl_a_coef[0],
                                         p_frac * CM_DATA_DICT[p_band].chl_a_coef[1] +
                                         n_frac * CM_DATA_DICT[n_band].chl_a_coef[1]],
                             total_abs=(p_frac * CM_DATA_DICT[p_band].total_abs +
                                        n_frac * CM_DATA_DICT[n_band].total_abs),
                             total_backscat=(p_frac * CM_DATA_DICT[p_band].total_backscat +
                                             n_frac * CM_DATA_DICT[n_band].total_backscat), )
    return new_wat_data


def run_oceancolor(water_data, oc_data=None):
    """Placeholder for adding Ocean Color CCI support."""
    # Interim values for Chlorophyll-a absorption and backscatter.
    # OC-CCI data contains actual values, so when that's implemented these two lines can be removed.
    if oc_data is None:
        water_data.chlabs = water_data.chl_a_coef[0] * (1.0 - np.exp(-1.61 * chlconc)) + \
                            water_data.chl_a_coef[1] * chlconc

        water_data.chlbsc = 0.02 * (0.5 - 0.25 * np.log10(chlconc)) * 0.55 / water_data.wavelength + 0.002

        water_data.total_backscat = water_data.base_backscat + water_data.chlbsc
        water_data.total_abs = water_data.base_abs + water_data.chlabs
    else:
        raise NotImplementedError("PyCoxMunk currently does not support Ocean Color data!")

    return water_data


def _compute_abcd(first, second):
    """Compute coefficients for Fresnel's equation."""
    a1 = np.sin(first - second)
    b1 = np.sin(first + second)
    c1 = np.tan(first - second)
    d1 = np.tan(first + second)

    # Catch very small values to prevent division by zero
    b1 = da.where(np.abs(b1) < dither_more, np.nan, b1)
    d1 = da.where(np.abs(d1) < dither_more, np.nan, d1)

    return a1, b1, c1, d1


def _set_arr_bytype(inarr, outarr, idx1, idx2, idx3=None):
    """Add data to an array based on type."""

    if idx3 is None:
        if type(inarr) == float:
            outarr[:, :, idx1, idx2] = inarr
        else:
            outarr[:, :, idx1, idx2] = inarr[:, :]
    else:
        if type(inarr) == float:
            outarr[:, :, idx1, idx2, idx3] = inarr
        else:
            outarr[:, :, idx1, idx2, idx3] = inarr[:, :]


def _calc_0d_dv(in_refl: CMReflectance,
                band_wvl: float,
                geom_info: CMSceneGeom,
                wind_info: CMSharedWind,
                water_data,
                oc_cci_data=None, ) -> CMReflectance:
    """Compute the 0d and dv terms."""
    qx_theta, qw_theta = gauss_leg_quadx(n_quad_theta, 0., np.pi / 2.)
    qx_phi, qw_phi = gauss_leg_quadx(n_quad_phi, 0., 2. * np.pi)
    qx_cos_sin_qw_theta = np.cos(qx_theta) * np.sin(qx_theta) * qw_theta

    # Need to initialise various values first
    tmp_geom = copy.deepcopy(geom_info)
    tmp_geom2 = copy.deepcopy(geom_info)

    tmp_saa = da.zeros((in_refl.rho.shape[0], in_refl.rho.shape[0], n_quad_theta, n_quad_phi))
    tmp_vaa = da.zeros_like(tmp_saa)
    tmp_sza = da.zeros_like(tmp_saa)
    tmp_vza = da.zeros_like(tmp_saa)
    tmp_raa = da.zeros_like(tmp_saa)
    tmp_sza2 = da.zeros_like(tmp_saa)
    tmp_vza2 = da.zeros_like(tmp_saa)

    tmp_u10 = da.zeros_like(tmp_saa)
    tmp_v10 = da.zeros_like(tmp_saa)

    for j in range(0, n_quad_theta):
        for k in range(0, n_quad_phi):
            _set_arr_bytype(geom_info.saa, tmp_saa, j, k)
            _set_arr_bytype(geom_info.vaa, tmp_vaa, j, k)
            _set_arr_bytype(geom_info.sza, tmp_sza, j, k)
            _set_arr_bytype(geom_info.vza, tmp_vza, j, k)

            tmp_raa[:, :, j, k] = da.rad2deg(qx_phi[k])

            tmp_vza[:, :, j, k] = da.rad2deg(qx_theta[j])

            tmp_sza2[:, :, j, k] = da.rad2deg(qx_theta[j])

            tmp_u10[:, :, j, k] = wind_info.u10
            tmp_v10[:, :, j, k] = wind_info.v10

    tmp_geom.sza = tmp_sza
    tmp_geom.vza = tmp_vza
    tmp_geom.vaa = tmp_vaa
    tmp_geom.saa = tmp_saa
    tmp_geom.raa = tmp_raa

    tmp_geom2.sza = tmp_sza2
    tmp_geom2.vza = tmp_vza2
    tmp_geom2.vaa = tmp_vaa
    tmp_geom2.saa = tmp_saa
    tmp_geom2.raa = tmp_raa

    tmp_geom.compute_additional()
    tmp_geom2.compute_additional()

    tmp_wind = CMSharedWind(tmp_geom, tmp_u10, tmp_v10)
    tmp_wind2 = CMSharedWind(tmp_geom2, tmp_u10, tmp_v10)

    # Now do the direct-diffuse and diffuse-direct terms.
    # Need to recompute some angles.

    # Compute reflectances
    tmp_cm_refl = calc_cox_munk(band_wvl, tmp_geom, tmp_wind, oc_cci_data, water_data)
    tmp_cm_refl2 = calc_cox_munk(band_wvl, tmp_geom2, tmp_wind2, oc_cci_data, water_data)

    tmp_rho = da.map_blocks(_loop_od_ov,
                            tmp_cm_refl.rho, tmp_cm_refl2.rho,
                            qw_phi, qx_cos_sin_qw_theta,
                            n_quad_theta,
                            dtype=np.float32,
                            chunks=(tmp_cm_refl.rho.chunks[0],
                                    tmp_cm_refl.rho.chunks[1],
                                    2,),
                            drop_axis=3)

    tmp_rho_0d = np.squeeze(tmp_rho[:, :, 0])
    tmp_rho_dv = np.squeeze(tmp_rho[:, :, 1])
    in_refl.rho_0d = tmp_rho_0d / np.pi
    in_refl.rho_dv = tmp_rho_dv / np.pi

    return in_refl


@jit(parallel=True, nopython=True)
def _loop_od_ov(rho: np.ndarray, rho2: np.ndarray,
                qw_phi: np.ndarray, qx_cos_sin_qw_theta: np.ndarray,
                nq_theta: int = 4, nq_phi: int = 4) -> np.ndarray:
    """Loop over quad theta and phi."""

    out_0d = np.zeros_like(rho[:, :, 0, 0])
    out_dv = np.zeros_like(rho[:, :, 0, 0])

    for j in range(0, nq_theta):
        aa1 = np.zeros_like(rho[:, :, 0, 0])
        aa2 = np.zeros_like(rho[:, :, 0, 0])
        for k in range(0, nq_phi):
            aa1 = aa1 + rho[:, :, j, k] * qw_phi[k]
            aa2 = aa2 + rho2[:, :, j, k] * qw_phi[k]
        out_0d = out_0d + aa1 * qx_cos_sin_qw_theta[j]
        out_dv = out_dv + aa2 * qx_cos_sin_qw_theta[j]

    tmp_rho = np.dstack((out_0d, out_dv))
    return tmp_rho


def _calc_dd(in_refl: CMReflectance,
             band_wvl: float,
             geom_info: CMSceneGeom,
             wind_info: CMSharedWind,
             water_data,
             oc_cci_data=None, ) -> CMReflectance:
    """Compute the dd term."""
    qx_theta, qw_theta = gauss_leg_quadx(n_quad_theta, 0., np.pi / 2.)
    qx_phi, qw_phi = gauss_leg_quadx(n_quad_phi, 0., 2. * np.pi)
    qx_cos_sin_qw_theta = np.cos(qx_theta) * np.sin(qx_theta) * qw_theta

    # Need to initialise various values first
    tmp_geom = copy.deepcopy(geom_info)

    tmp_saa = da.zeros((in_refl.rho.shape[0], in_refl.rho.shape[0], n_quad_theta, n_quad_theta, n_quad_phi))
    tmp_vaa = da.zeros_like(tmp_saa)
    tmp_sza = da.zeros_like(tmp_saa)
    tmp_vza = da.zeros_like(tmp_saa)
    tmp_raa = da.zeros_like(tmp_saa)

    if type(wind_info.u10) is not float:
        tmp_u10 = da.zeros((in_refl.rho.shape[0], in_refl.rho.shape[0], n_quad_theta, n_quad_theta, n_quad_phi))
    if type(wind_info.v10) is not float:
        tmp_v10 = da.zeros((in_refl.rho.shape[0], in_refl.rho.shape[0], n_quad_theta, n_quad_theta, n_quad_phi))

    for i in range(0, n_quad_theta):
        for j in range(0, n_quad_theta):
            for k in range(0, n_quad_phi):
                _set_arr_bytype(geom_info.saa, tmp_saa, i, j, k)
                _set_arr_bytype(geom_info.vaa, tmp_vaa, i, j, k)
                _set_arr_bytype(geom_info.sza, tmp_sza, i, j, k)
                _set_arr_bytype(geom_info.vza, tmp_vza, i, j, k)

                tmp_raa[:, :, i, j, k] = da.rad2deg(qx_phi[k])
                tmp_vza[:, :, i, j, k] = da.rad2deg(qx_theta[j])
                tmp_sza[:, :, i, j, k] = da.rad2deg(qx_theta[i])

                tmp_u10[:, :, i, j, k] = wind_info.u10
                tmp_v10[:, :, i, j, k] = wind_info.v10

    tmp_geom.sza = tmp_sza
    tmp_geom.vza = tmp_vza
    tmp_geom.vaa = tmp_vaa
    tmp_geom.saa = tmp_saa
    tmp_geom.raa = tmp_raa

    tmp_geom.compute_additional()

    tmp_wind = CMSharedWind(tmp_geom, tmp_u10, tmp_v10)

    # Now do the direct-diffuse and diffuse-direct terms.
    # Need to recompute some angles.

    # Compute reflectances
    tmp_cm_refl = calc_cox_munk(band_wvl, tmp_geom, tmp_wind, oc_cci_data, water_data)
    tmp_rho_dd = da.map_blocks(_loop_dd,
                               tmp_cm_refl.rho,
                               qw_phi, qx_cos_sin_qw_theta,
                               n_quad_theta, n_quad_phi,
                               dtype='f2',
                               chunks=(tmp_cm_refl.rho.chunks[0],
                                       tmp_cm_refl.rho.chunks[1]),
                               drop_axis=[2, 3, 4, ])

    in_refl.rho_dd = tmp_rho_dd * 2. / np.pi

    return in_refl


@jit(parallel=True, nopython=True)
def _loop_dd(rho: np.ndarray,
             qw_phi: np.ndarray, qx_cos_sin_qw_theta: np.ndarray,
             nq_theta: int = 4, nq_phi: int = 4) -> np.ndarray:
    """Loop over quad theta and phi."""

    out_dd = np.zeros_like(rho[:, :, 0, 0, 0])

    for i in range(0, nq_theta):
        aa1 = np.zeros_like(rho[:, :, 0, 0, 0])
        for j in range(0, nq_theta):
            aa2 = np.zeros_like(rho[:, :, 0, 0, 0])
            for k in range(0, nq_phi):
                aa2 = aa2 + rho[:, :, i, j, k] * qw_phi[k]
            aa1 = aa1 + aa2 * qx_cos_sin_qw_theta[j]
        out_dd = out_dd + aa1 * qx_cos_sin_qw_theta[i]

    return out_dd


def calc_cox_munk_brdf_terms(in_refl: CMReflectance,
                             band_wvl: float,
                             geom_info: CMSceneGeom,
                             wind_info: CMSharedWind,
                             oc_cci_data=None) -> CMReflectance:
    """Compute the BRDF terms for the sea surface reflectance.
    These terms are:
     - rho_0v: Solar beam to satellite view reflectances
     - rho_0d: Solar beam to diffuse reflectances
     - rho_dv: Diffuse to satellite view reflectances
     - rho_dd: Diffuse to diffuse reflectances
    Inputs:
      - band_wvl: Float, central wavelength of the band being processed.
      - geom_info: CMSceneGeom, the scene geometry.
      - wind_info: CMSharedWind, the wind-based information for the scene.
      - pix_mask: CMPixMask, the pixel masks for the Scene
      - oc_cci_data: None, placeholder for when this data is used.
    Returns:
      - coxmunk_data: CM_Reflectance, output reflectances and BRDF components.
     """

    # Initialise remaining reflectances
    in_refl.rho_0d = da.zeros_like(in_refl.rho)
    in_refl.rho_dv = da.zeros_like(in_refl.rho)
    in_refl.rho_dd = da.zeros_like(in_refl.rho)

    # Get specific information about water properties at wavelength
    water_data = compute_wavelength_specific_water_props(band_wvl, meth='interp')
    # Apply ocean color CCI data (currently not available)
    water_data = run_oceancolor(water_data, oc_cci_data)

    # The direct-direct term is simply the sea surface reflectance
    # That has already been computed by `calc_cox_munk`
    in_refl.rho_0v = in_refl.rho.copy()

    in_refl = _calc_0d_dv(in_refl, band_wvl, geom_info, wind_info, water_data, oc_cci_data)

    # Lastly, do the diffuse-diffuse term.
    in_refl = _calc_dd(in_refl, band_wvl, geom_info, wind_info, water_data, oc_cci_data)

    return in_refl


def calc_cox_munk(band_wvl: float,
                  geom_info: CMSceneGeom,
                  wind_info: CMSharedWind,
                  oc_cci_data=None,
                  water_data=None) -> CMReflectance:
    """Compute the bidirectional reflectance from scene information using the Cox-Munk approach.
    Currently, this uses only the band central wavelength, not spectral response, and doesn't include the
    ability to process with Ocean Color CCI (or other ocean color) datasets.
    This is typically called from within the PyCoxMunk class but can be called directly by external code.
    Inputs:
      - band_wvl: Float, central wavelength of the band being processed.
      - geom_info: CMSceneGeom, the scene geometry.
      - wind_info: CMSharedWind, the wind-based information for the scene.
      - pix_mask: CMPixMask, the pixel masks for the Scene
      - oc_cci_data: None, placeholder for when this data is used.
      - water_data: None or WaterData, information about the water body spectrum.
    Returns:
      - coxmunk_data: CM_Reflectance, output reflectances and, if requested, BRDF components.
    """

    # Get specific information about water properties at wavelength
    if water_data is None:
        water_data = compute_wavelength_specific_water_props(band_wvl, meth='interp')
        # Apply ocean color CCI data (currently not available)
        water_data = run_oceancolor(water_data, oc_cci_data)

    # Calculate white cap reflectance
    rhowc = wind_info.wcfrac * water_data.whitecap_refl

    # Reflectance of water body using total backscatter and absorption, ideally derived from Ocean Color CCI
    # data but currently this is not implemented.
    eta_oc = 0.5 * water_data.base_backscat / water_data.total_backscat

    # Coefficient of R
    f = 0.6279 - (0.2227 * eta_oc) - (0.00513 * eta_oc) + (0.2465 * eta_oc - 0.3119) * geom_info.cos_sza

    # Water body reflectance
    refl_water = f * water_data.total_backscat / water_data.total_abs

    # Use Fresnel equation and Snell's law to determine how much light enters water body from surface.
    # Upward transmission and reflectance assumed constant for wavelengths where underlight is significant
    t_u = 0.52
    r_u = 1.0 - t_u

    w = np.deg2rad(geom_info.vza)
    wprime = np.arcsin(n_air * geom_info.sin_sza / water_data.refrac_real)
    a1, b1, c1, d1 = _compute_abcd(w, wprime)

    t_d = 1.0 - 0.5 * ((a1 * a1) / (b1 * b1) + (c1 * c1) / (d1 * d1))
    # Deal with the NaN values we potentially introduced
    t_d = da.where(np.isnan(t_d), 0, t_d)

    # Combine surface transmission with underlight to give total underlight contribution
    rhoul = (t_u + t_d + refl_water) / (1.0 - r_u * refl_water)

    # Now we can move on to calculate the  reflectance, which is the actual Cox-Munk step
    wprime = np.arcsin((n_air * wind_info.sin_w) / water_data.refrac_real)
    a1, b1, c1, d1 = _compute_abcd(wind_info.w, wprime)

    r_sf = 1.0 - 0.5 * ((a1 * a1) / (b1 * b1) + (c1 * c1) / (d1 * d1))
    # Deal with the NaN values we potentially introduced
    r_sf = da.where(da.isnan(r_sf), 0, r_sf)

    # Calculate glint contribution
    rhogl = da.where(da.abs(wind_info.a) > dither_more,
                     np.pi * wind_info.p * r_sf / wind_info.a, 0)

    # Calculate overall reflectance
    rho = rhowc + (1 - wind_info.wcfrac) * (rhogl + rhoul)

    # Add results to the output class
    coxmunk_data = CMReflectance(band_wvl, rho, rhowc, rhogl, rhoul, None, None, None, None)

    return coxmunk_data


def calc_coxmunk_wrapper(band_wvl: float,
                         geom_info: CMSceneGeom,
                         wind_info: CMSharedWind,
                         oc_cci_data=None,
                         do_brdf: bool = False) -> CMReflectance:
    """Wrapper for the above functions that combines the C-M and BRDF calcs."""

    # If we want to compute BRDF terms then call the correct function now.
    coxmunk_data = calc_cox_munk(band_wvl, geom_info, wind_info, oc_cci_data)
    if do_brdf:
        coxmunk_data = calc_cox_munk_brdf_terms(coxmunk_data, band_wvl, geom_info, wind_info, oc_cci_data)

    return coxmunk_data
