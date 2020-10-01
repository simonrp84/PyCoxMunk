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
"""Class for the Scene geometry information."""

import numpy as np


class CMSceneGeom:
    """A class for data about the scene geometry.

    This stores the viewing and solar zenith and azimuth angles, plus the relative azimuth,
    latitude and longitude for each pixel / input point.
    All input data should be numpy arrays of the same size or should be a single float.
    If a single value is passed while all other input are numpy arrays, the single value
    will be converted into a numpy array of constant values with size equal to the other inputs.

    ***Notes on the geometry:***
    All angles should be supplied in degrees.

    Zenith angles should be in the form where:
        - 0.0 zenith means overhead
        - 90.0 zenith means on the horizon
        - 180.0 is the extreme case of below the horizon (on other side of Earth)

    Solar and Viewing azimuth angles are measured clockwise from North:
        - 0.0 azimuth means sun or sensor is due North of pixel
        - 90.0 azimuth means sun or sensor is due East of pixel
        - 180.0 azimuth means sun or sensor is due South of pixel
        - 270.0 azimuth means sun or sensor is due West of pixel

    Relative azimuth angles are defined where:
        - 0.0 RAA means sun and sensor are on same side of pixel (backscatter)
        - 180.0 RAA means sun and sensor are on opposite sides of pixel (forward scatter)

    """

    # Specify some hard-coded limits
    zenith_min = 0.
    zenith_max = 180.

    azimuth_min = 0.
    azimuth_max = 360.

    relazi_min = 0.
    relazi_max = 180.

    lat_min = -90.
    lat_max = 90.

    lon_min = -180.
    lon_max = 180.

    @staticmethod
    def _check_and_reshape(arr, good_shape):
        """If array is single value then scale to match other arrays."""
        if arr.shape == (1,):
            arr = np.full(good_shape, arr[0])
        return arr

    @staticmethod
    def _check_type(in_val, var_typ):
        """Check that input variable is correct type.
        All inputs should be numpy array or a float."""
        if type(in_val) == float:
            return np.array([in_val])
        elif type(in_val) == np.ndarray:
            return in_val
        else:
            raise TypeError(var_typ + " must be a single float or numpy array!")

    def __init__(self, sza, saa, vza, vaa, lats, lons, raa=None):
        # Solar zenith angles
        self.sza = self._check_type(sza, "Solar zenith angles")
        # Solar azimuth angles
        self.saa = self._check_type(saa, "Solar azimuth angles")
        # Viewing zenith angles
        self.vza = self._check_type(vza, "Viewing zenith angles")
        # Viewing azimuth angles
        self.vaa = self._check_type(vaa, "Viewing azimuth angles")
        # Latitudes
        self.lats = self._check_type(lats, "Latitudes")
        # Longitudes
        self.lons = self._check_type(lons, "Longitudes")

        # Check arrays are equal size or are single values
        self.check_array_shapes()

        # Relative azimuth angle (range: 0 -> 180)
        # Can be provided by user or calculated internally
        if raa is not None:
            raa = self._check_type(raa, "Relative azimuth angles")
            self.raa = self._check_and_reshape(raa, self.sza.shape)
        else:
            self.raa = self.calc_relazi()

        # Run checks to ensure values are within bounds
        self.check_angle_bounds()

    def check_array_shapes(self):
        """Ensure that input arrays have same shape as each other."""
        # First, find which variables are single values
        in_arr = [self.sza, self.saa, self.vza, self.vaa, self.lats, self.lons]
        arr_shapes = [in_arr[0].shape == (1,), in_arr[1].shape == (1,),
                      in_arr[2].shape == (1,), in_arr[3].shape == (1,),
                      in_arr[4].shape == (1,), in_arr[5].shape == (1,)]
        # If everything is a single value, this is OK
        if all(arr_shapes):
            return
        # If everything is an array, also OK but we need to check shapes
        elif not any(arr_shapes):
            # Assume SZA shape is the reference
            good_shp = self.sza.shape
            if (self.saa.shape != good_shp or
                    self.vza.shape != good_shp or
                    self.vaa.shape != good_shp or
                    self.lats.shape != good_shp or
                    self.lons.shape != good_shp):
                raise ValueError("All input arrays must have identical shape or be single values.")
        # If we have a mixture of single values and larger arrays
        else:
            good_shp = None
            for i in range(0, len(in_arr)):
                if not arr_shapes[i]:
                    good_shp = in_arr[i].shape
                    break
            self.saa = self._check_and_reshape(self.saa, good_shp)
            self.vza = self._check_and_reshape(self.vza, good_shp)
            self.vaa = self._check_and_reshape(self.vaa, good_shp)
            self.lats = self._check_and_reshape(self.lats, good_shp)
            self.lons = self._check_and_reshape(self.lons, good_shp)

    def check_angle_bounds(self):
        """Ensure that zenith angles are in acceptable range."""
        # Solar zenith
        if np.any(self.sza < self.zenith_min) or np.any(self.sza >= self.zenith_max):
            raise ValueError("All solar zenith angles must be in the range " +
                             str(self.zenith_min) + ' to ' + str(self.zenith_max))
        # Viewing zenith
        if np.any(self.vza < self.zenith_min) or np.any(self.vza >= self.zenith_max):
            raise ValueError("All viewing zenith angles must be in the range " +
                             str(self.zenith_min) + ' to ' + str(self.zenith_max))
        # Solar azimuth
        if np.any(self.saa < self.azimuth_min) or np.any(self.saa >= self.azimuth_max):
            raise ValueError("All solar azimuth angles must be in the range " +
                             str(self.azimuth_min) + ' to ' + str(self.azimuth_max))
        # Viewing azimuth
        if np.any(self.vaa < self.azimuth_min) or np.any(self.vaa >= self.azimuth_max):
            raise ValueError("All viewing azimuth angles must be in the range " +
                             str(self.azimuth_min) + ' to ' + str(self.azimuth_max))
        # Latitudes
        if np.any(self.lats < self.lat_min) or np.any(self.lats >= self.lat_max):
            raise ValueError("All latitudes must be in the range " +
                             str(self.lat_min) + ' to ' + str(self.lat_max))
        # Longitudes
        if np.any(self.lons < self.lon_min) or np.any(self.lons >= self.lon_max):
            raise ValueError("All longitudes must be in the range " +
                             str(self.lon_min) + ' to ' + str(self.lon_max))
        # Relative azimuth
        if np.any(self.raa < self.relazi_min) or np.any(self.raa >= self.relazi_max):
            raise ValueError("All Relative azimuth angles must be in the range " +
                             str(self.azimuth_min) + ' to ' + str(self.relazi_max))

    def calc_relazi(self):
        """Compute the relative azimuth angles from solar and viewing azimuths."""
        raa = np.abs(self.vaa - self.saa)

        # If any RAA values are greater than 180, invert them.
        raa = np.where(raa >= 180., 360. - raa, raa)
        return raa
