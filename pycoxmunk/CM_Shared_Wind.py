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

from pycoxmunk.CM_Utils import check_and_reshape, check_type
from pycoxmunk.CM_SceneGeom import CMSceneGeom
from pycoxmunk.CM_Constants import dither_more
import dask.array as da
import numpy as np


class CMSharedWind:
    """A class for Cox-Munk calculations involving wind speed.
    These are typically shared across the Cox-Munk code and are
    computed from scene geometry and wind speed."""

    @staticmethod
    def zeisse_ba(theta, cos_theta, ws):
        """Applies Zeisse correction for zenith angles greater than 70.
        This replaces cos(vza) with the area of the so-called 'ergodic cap'.
        See Zeisse (1995) for more details.
        Inputs:
          - theta: Float or np.ndarray, satellite zenith angles in radians.
          - cos_theta: Float or np.ndarray, cosine of satellite zenith angles.
          - ws: Float or np.ndarray, the absolute 10m wind speed in m/s
        Returns:
          - float or np.ndarray, cos(vza) or ergodic cap area, depending on vza.
        """
        from pycoxmunk.CM_Constants import zeisse_coef as zc

        tp = (theta - da.deg2rad(70)) / 5.
        wp = 4. * da.log10(ws) / 1.30103

        delta = (zc[0, 0] + (zc[1, 0] + zc[2, 0] * wp) * wp) + \
                ((zc[0, 1] + (zc[1, 1] + zc[2, 1] * wp) * wp) +
                 (zc[0, 2] + (zc[1, 2] + zc[2, 2] * wp) * wp) * tp) * tp

        ba = da.where((theta >= da.deg2rad(70)) & (ws > 1.), delta + cos_theta, cos_theta)

        return ba

    def __init__(self,
                 scenegeom: CMSceneGeom,
                 u10,
                 v10):
        """A class for Cox-Munk calculations involving wind speed.
        Inputs:
          - scenegeom: CMSceneGeom, class containing scene geometry information.
          - u10: Float or np.ndarray, u-direction wind speed at 10m in m/s.
          - v10: Float or np.ndarray, v-direction wind speed at 10m in m/s.
        Returns:
          - Class containing calculated values.
        """
        # Wind speeds
        self.u10 = check_type(u10, "U-direction 10m wind speed")
        self.v10 = check_type(v10, "V-direction 10m wind speed")

        # Here we assume that the scene geometry already has correctly sized arrays
        self.u10 = check_and_reshape(self.u10, scenegeom.sza.shape)
        self.v10 = check_and_reshape(self.v10, scenegeom.sza.shape)

        # Calculate parameters
        # First 10m absolute wind speed
        self.ws = da.sqrt(self.u10 * self.u10 + self.v10 * self.v10)

        # Wind direction
        self.wd = da.arccos(self.v10 / self.ws)
        # We want wind direction in -180 -> 180 range
        self.wd = da.where(self.wd < 0, -self.wd, self.wd)

        # Consistency check
        self.wd = da.where(self.ws < 0, 0, self.wd)
        self.ws = da.where(self.ws < 0, 0, self.ws)

        # White-cap fraction
        self.wcfrac = 2.951e-6 * da.power(self.ws, 3.52)
        self.wcfrac = da.where(self.wcfrac > 1, 1, self.wcfrac)

        # Convert wind direction to be relative to solar azimuth
        self.wd = self.wd - da.deg2rad(scenegeom.saa)
        self.wd = da.where(self.wd < 0, 2 * np.pi + self.wd, self.wd)

        # Angle between incident light and instrument with respect to
        # the sloping sea surface.
        self.cosomega = (scenegeom.cos_sza * scenegeom.cos_vza) + \
                        (scenegeom.sin_sza * scenegeom.sin_vza * scenegeom.cos_raa)

        self.cosbeta = da.where(da.abs(self.cosomega) + 1. > dither_more,
                                (scenegeom.cos_sza + scenegeom.cos_vza) / (da.sqrt(2 + 2 * self.cosomega)),
                                0.)

        self.w = 0.5 * da.arccos(self.cosomega)
        self.sin_w = da.sin(self.w)

        # Zeisse correction
        self.ergodic = self.zeisse_ba(da.deg2rad(scenegeom.vza), scenegeom.cos_vza, self.ws)
        self.a = 4 * scenegeom.cos_sza * self.ergodic * da.power(self.cosbeta, 4)

        # Surface slopes
        self.dangle = scenegeom.cos_sza + scenegeom.cos_vza

        self.Zx = da.where(self.dangle > dither_more,
                           -1. * (scenegeom.sin_vza * scenegeom.sin_raa) / self.dangle,
                           0.)
        self.Zy = da.where(self.dangle > dither_more,
                           -1. * (scenegeom.sin_sza + scenegeom.sin_vza * scenegeom.cos_raa) / self.dangle,
                           0.)

        self.cos_wd = da.cos(self.wd)
        self.sin_wd = da.sin(self.wd)

        self.Zxprime = self.cos_wd * self.Zx + self.sin_wd * self.Zy
        self.Zyprime = -self.sin_wd * self.Zx + self.cos_wd * self.Zy

        # Coefficients for Cox-Munk
        self.sigx = da.sqrt(0.003 + 0.00192 * self.ws)
        self.sigy = da.sqrt(0.00316 * self.ws)

        self.zeta = self.Zxprime / self.sigx
        self.eta = self.Zyprime / self.sigy

        self.p = da.exp(-0.5 * (self.zeta * self.zeta + self.eta * self.eta)) / (2. * np.pi * self.sigx * self.sigy)
