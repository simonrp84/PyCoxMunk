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

from pycoxmunk.CM_SceneGeom import CMSceneGeom, cm_calcangles
from pycoxmunk.CM_Calcs import calc_cox_munk, CM_Reflectance
from pycoxmunk.CM_Shared_Wind import CMSharedWind
from satpy import Scene
import numpy as np


class PyCoxMunk:
    """The main class for the library, sets up and runs processing."""

    def __init__(self, scn, band_names, oc_dir=None, angle_names=None,
                 do_brdf=False, mask_bad=True, delete_when_done=True):
        """Initialise the class.
        Inputs:
        - scn: Satpy Scene, the scene containing data and angles.
        - band_names: List, if supplied lists all band names to process.
        - oc_dir: (optional) String, if supplied points to location of Ocean Color CCI data.
        - angle_names: (optional) Dict, if supplied this should list all
          solar/satellite angle dataset names. If not supplied, these are assumed.
        - do_brdf: Bool, if true then PyCoxMunk will also compute BRDF coefficients.
        - mask_bad: Bool, if true then pixels with bad data (f.ex zenith too high) will be set to np.nan
        """

        # Check types and set up variables
        if type(do_brdf) is bool:
            self.do_brdf = do_brdf
        else:
            raise ValueError("do_brdf variable must be boolean!")

        if type(mask_bad) is bool:
            self.mask_bad = mask_bad
        else:
            raise ValueError("mask_bad variable must be boolean!")

        if type(delete_when_done) is bool:
            self.delete_when_done = delete_when_done
        else:
            raise ValueError("delete_when_done variable must be boolean!")

        if type(scn) is not Scene:
            raise ValueError("input scene variable must be a satpy Scene!")
        else:
            self.scn = scn

        if type(band_names) is not list:
            raise ValueError("band_names variable must be a list!")
        elif len(band_names) < 1:
            raise ValueError("band_names must contain at least one band!")
        else:
            self.band_names = band_names

        if angle_names is not None and type(angle_names) is not dict and type(angle_names) is not str:
            raise ValueError("angle_names variable must be None, a dict or a string!")
        elif angle_names is None or angle_names == 'calc':
            self.angle_names = {'sza': 'solar_zenith_angle',
                                'vza': 'satellite_zenith_angle',
                                'saa': 'solar_azimuth_angle',
                                'vaa': 'satellite_azimuth_angle'}
        elif 'sza' not in angle_names.keys():
            raise ValueError("angle_names dict must contain 'sza' key!")
        elif 'vza' not in angle_names.keys():
            raise ValueError("angle_names dict must contain 'vza' key!")
        elif 'saa' not in angle_names.keys():
            raise ValueError("angle_names dict must contain 'saa' key!")
        elif 'vaa' not in angle_names.keys():
            raise ValueError("angle_names dict must contain 'vaa' key!")
        else:
            self.angle_names = angle_names
        if angle_names == 'calc':
            self.scn = cm_calcangles(scn, refband=band_names[0])

        if oc_dir is None:
            self.use_occci = False
        elif type(oc_dir) is not str:
            raise ValueError("oc_dir variable must be None or a string!")
        else:
            self.oc_dir = oc_dir
            self.use_occci = True

        # If angle and band names are supplied then check these are actually in the scene
        for bname in self.band_names:
            if bname not in self.scn:
                raise KeyError(f"User-supplied band dataset {bname} not in input scene!")

        if self.angle_names is not None:
            for kname in self.angle_names.keys():
                bname = self.angle_names[kname]
                if bname not in self.scn:
                    print(bname)
                    raise KeyError(f"User-supplied angle dataset {bname} not in input scene!")

        # Cox Munk output class, used later
        self.cm_refl = CM_Reflectance()

        # Compute longitudes and latitudes from first selected band. Assumes all bands are same
        # dimensions! This means that high res bands (such as Himawari B03) must be resampled to
        # lower res band resolution before passing to PyCoxMunk. If you wish to process both
        # high and low res bands at native resolution then you must call PyCoxMunk multiple times.
        lons, lats = self.scn[self.band_names[0]].attrs['area'].get_lonlats()
        lons = np.array(lons)
        lats = np.array(lats)
        self.geometry = CMSceneGeom(np.array(self.scn[self.angle_names['sza']]),
                                    np.array(self.scn[self.angle_names['saa']]),
                                    np.array(self.scn[self.angle_names['vza']]),
                                    np.array(self.scn[self.angle_names['vaa']]),
                                    lats, lons)

        # Initialise shared winds, not yet loaded
        self.shared_wind = None

    def setup_wind(self, u10, v10):
        """Set up the fields that depend on wind speed.
        Inputs:
          - u10: Float or np.ndarray, u-direction wind speed at 10m in m/s.
          - v10: Float or np.ndarray, v-direction wind speed at 10m in m/s.
        Returns:
          - self.shared_wind: CM_Shared_Wind class.
        """

        self.shared_wind = CMSharedWind(self.geometry, u10, v10)

    def retr_coxmunk_refl(self):
        for band_id in self.band_names:
            out_band_id = f'cox_munk_refl_{band_id}'
            self.scn[out_band_id] = self.scn[band_id].copy()
            self.cm_refl = calc_cox_munk(self.scn[band_id].attrs['wavelength'].central,
                                         self.geometry,
                                         self.shared_wind)

            self.scn[out_band_id].data = self.cm_refl.rho

            if self.do_brdf:
                out_band_id = f'cox_munk_rho0d_{band_id}'
                self.scn[out_band_id] = self.scn[band_id].copy()
                self.scn[out_band_id].data = self.cm_refl.rho_0d

                out_band_id = f'cox_munk_rho0v_{band_id}'
                self.scn[out_band_id] = self.scn[band_id].copy()
                self.scn[out_band_id].data = self.cm_refl.rho_0v

                out_band_id = f'cox_munk_rhodv_{band_id}'
                self.scn[out_band_id] = self.scn[band_id].copy()
                self.scn[out_band_id].data = self.cm_refl.rho_dv

                out_band_id = f'cox_munk_rhodd_{band_id}'
                self.scn[out_band_id] = self.scn[band_id].copy()
                self.scn[out_band_id].data = self.cm_refl.rhodd

            # Mask bad pixels
            if self.mask_bad:
                masker_rho = np.argwhere(self.cm_refl.rho < - 0.5)
                masker_sza = np.argwhere(self.geometry.sza >= self.geometry.zenith_max)
                masker_vza = np.argwhere(self.geometry.vza >= self.geometry.zenith_max)

                self.cm_refl.rho[masker_rho] = np.nan
                self.cm_refl.rho[masker_sza] = np.nan
                self.cm_refl.rho[masker_vza] = np.nan
                if self.do_brdf:
                    self.cm_refl.rho_0d[masker_rho] = np.nan
                    self.cm_refl.rho_0v[masker_rho] = np.nan
                    self.cm_refl.rho_dd[masker_rho] = np.nan
                    self.cm_refl.rho_dv[masker_rho] = np.nan
                    self.cm_refl.rho_0d[masker_sza] = np.nan
                    self.cm_refl.rho_0v[masker_sza] = np.nan
                    self.cm_refl.rho_dd[masker_sza] = np.nan
                    self.cm_refl.rho_dv[masker_sza] = np.nan
                    self.cm_refl.rho_0d[masker_vza] = np.nan
                    self.cm_refl.rho_0v[masker_vza] = np.nan
                    self.cm_refl.rho_dd[masker_vza] = np.nan
                    self.cm_refl.rho_dv[masker_vza] = np.nan

        if self.delete_when_done:
            del self.shared_wind
            del self.geometry
            del self.cm_refl.rhoul
            del self.cm_refl.rhogl
            del self.cm_refl.rhowc
