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
from pycoxmunk.CM_Calcs import calc_coxmunk_wrapper, CM_Reflectance
from pycoxmunk.CM_Shared_Wind import CMSharedWind
from pycoxmunk.CM_PixMask import CMPixMask
from satpy import Scene
import dask.array as da
import numpy as np


class PyCoxMunk:
    """The main class for the library, sets up and runs processing.

    This code will estimate the sea surface reflectance as-seen from a given satellite
    instrument. It is closely integrated with the `satpy` library, which is used for reading
    the original satellite data and to store the reflectance output.

    To use, first load your data using a satpy Scene. Then, call this class giving the scene
    as an argument together with a list of band names to be processed.
    Optionally, you can provide an `oc_dir` that contains ocean color CCI files for use in estimating the
    reflectance.
    `angle_names` is optional and if used should be a dict of satellite and solar angle datasets within
    the Scene. The form is:
    ```
    {'sza': 'sza_name', 'vza': 'vza_name', 'saa': 'saa_name', 'vaa': 'vaa_name'}
    ```
    Where `sza_name` is the name of the solar zenith angle dataset within the Scene, `vza_name` is for
    the viewing zenith angle and `saa_name` plus `vaa_name` are the names of the solar and viewing azimuth
    angle datasets respectively.
    `do_brdf` is an optional argument to specify whether BRDF properties are calculated. This is
    `False` by default.
    `mask_bad` is optional and enables/disables masking data of bad quality, typically reflectance values
    that are unphysically low. It is `True` by default.
    `delete_when_done` is an optional boolean to specify whether intermediate datasets are deleted once
    the sea surface reflectance is calculated. This can help reduce memory use and is `True` by default.

    Once you have initialised this class, you can optionally add data masks with `setup_pixmask` and add
    wind information, which improves accuracy of the estimated reflectance, with `setup_wind`.

    The actual reflectances can then be computed with `retr_coxmunk_refl`. Results will then be available
    in the `scn` variable.
    Example:
    ```
    my_band = 'VIS800'

    my_scn = satpy.Scene(some_file)
    my_scn.load([my_band])
    my_pcm = pycoxmunk.PyCoxMunk(my_scn, [my_band], do_brdf=True)
    my_pcm.setup_wind(my_u10_wind, my_v10_wind)
    my_pcm.retr_coxmunk_refl()

    my_sea_refl = my_pcm.scn['cox_munk_refl_VIS800'].data
    """

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
                    raise KeyError(f"User-supplied angle dataset {bname} not in input scene!")

        # Cox Munk output class, used later
        self.cm_refl = CM_Reflectance()

        # Compute longitudes and latitudes from first selected band. Assumes all bands are same
        # dimensions! This means that high res bands (such as Himawari B03) must be resampled to
        # lower res band resolution before passing to PyCoxMunk. If you wish to process both
        # high and low res bands at native resolution then you must call PyCoxMunk multiple times.
        lons, lats = self.scn[self.band_names[0]].attrs['area'].get_lonlats_dask()

        self.geometry = CMSceneGeom(da.array(self.scn[self.angle_names['sza']]),
                                    da.array(self.scn[self.angle_names['saa']]),
                                    da.array(self.scn[self.angle_names['vza']]),
                                    da.array(self.scn[self.angle_names['vaa']]),
                                    da.array(lats),
                                    da.array(lons))

        # Initialise shared winds and mask, not yet loaded
        self.shared_wind = None
        self.pixmask = None

    def setup_pixmask(self, cloud_mask=None, land_mask=None, sol_zen_mask=None, sat_zen_mask=None):
        """Add a pixel mask to the class.
        Inputs (all optional):
        - cloud_mask: An array for masking clouds, any values > 0 are assumed cloudy and not processed.
        - land_mask: An array for masking land pixels, any values > 0 are assumed cloudy and not processed.
        - sol_zen_mask: An array for masking high solar zeniths, any values > 0 are assumed cloudy and not processed.
        - sat_zen_mask: An array for masking high viewing zeniths, any values > 0 are assumed cloudy and not processed.
        Zenith angle masks can be added post-creation by calling the `cut_high_zen` function in the PyCoxMunk.pixmask.
        """
        self.pixmask = CMPixMask(cloud_mask, land_mask,
                                 sol_zen_mask, sat_zen_mask)

    def setup_wind(self, u10, v10):
        """Set up the fields that depend on wind speed.
        Inputs:
          - u10: Float or darray, u-direction wind speed at 10m in m/s.
          - v10: Float or array, v-direction wind speed at 10m in m/s.
        Returns:
          - self.shared_wind: CM_Shared_Wind class.
        """

        self.shared_wind = CMSharedWind(self.geometry, u10, v10)

    def retr_coxmunk_refl(self):
        """Main function for computing the sea surface reflectance.

        This uses data previously loaded including, if available, wind and ocean color
        data to compute the reflectance using the approach described in Sayer (2010).
        DOI: 10.5194/amt-3-813-2010

        This accounts for white caps and Chlorophyll content (via Ocean Color). However,
        it is designed for use over open oceans and may give poor results in coastal regions,
        lakes, and other water types.
        """
        for band_id in self.band_names:
            out_band_id = f'cox_munk_refl_{band_id}'
            self.scn[out_band_id] = self.scn[band_id].copy()
            self.cm_refl = calc_coxmunk_wrapper(self.scn[band_id].attrs['wavelength'].central,
                                                self.geometry,
                                                self.shared_wind,
                                                oc_cci_data=None,
                                                do_brdf=self.do_brdf)
            # Mask bad pixels
            mlist = []
            if self.mask_bad:
                masker_rho = da.where(self.cm_refl.rho < - 0.5, np.nan, 1)
                mlist = [masker_rho]

            if self.pixmask is not None:
                pmask = da.where(self.pixmask.mask >= 1, np.nan, 1)
                mlist.append(pmask)
            for masker in mlist:
                self.cm_refl.rho = self.cm_refl.rho * masker
                if self.do_brdf:
                    self.cm_refl.rho_0d = self.cm_refl.rho_0d * masker
                    self.cm_refl.rho_0v = self.cm_refl.rho_0v * masker
                    self.cm_refl.rho_dd = self.cm_refl.rho_dd * masker
                    self.cm_refl.rho_dv = self.cm_refl.rho_dv * masker

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
                self.scn[out_band_id].data = self.cm_refl.rho_dd

        if self.delete_when_done:
            try:
                del self.shared_wind
            except NameError:
                pass
            try:
                del self.geometry
            except NameError:
                pass
            try:
                del self.cm_refl.rhoul
            except NameError:
                pass
            try:
                del self.cm_refl.rhogl
            except NameError:
                pass
            try:
                del self.cm_refl.rhowc
            except NameError:
                pass
