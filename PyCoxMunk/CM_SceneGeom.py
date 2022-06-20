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

from PyCoxMunk.CM_Utils import check_and_reshape, check_type
import numpy as np
import warnings


def cm_calcangles(inscn, refband):
    """Calculate satellite and solar angles from a given dataset.
    NOTE: This only works for geostationary satellites!
    Inputs:
      - inscn: Scene, containing datasets.
      - refband: Str, a reference dataset name used to get area definition.
    Returns:
      - Modified scene containing angles.
    """
    from satpy.modifiers.angles import _get_sensor_angles, _get_sun_angles

    try:
        suna, sunz = _get_sun_angles(inscn[refband])
        sata, satz = _get_sensor_angles(inscn[refband])

        inscn['satellite_azimuth_angle'] = inscn[refband].copy()
        inscn['solar_azimuth_angle'] = inscn[refband].copy()
        inscn['satellite_zenith_angle'] = inscn[refband].copy()
        inscn['solar_zenith_angle'] = inscn[refband].copy()

        inscn['satellite_azimuth_angle'].data = sata
        inscn['solar_azimuth_angle'].data = suna
        inscn['satellite_zenith_angle'].data = satz
        inscn['solar_zenith_angle'].data = sunz
    except KeyError:
        raise KeyError("Input scene does not contain reference dataset, please check inputs.")
    except ValueError:
        raise ValueError("Cannot retrieve solar and satellite angles, please compute manually and specify.")
    return inscn


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

    def __init__(self, sza, saa, vza, vaa, lats, lons, raa=None,
                 zenith_min=0., zenith_max=85.,
                 azimuth_min=0., azimuth_max=360.,
                 relazi_min=0., relazi_max=180.,
                 lat_min=-90., lat_max=90.,
                 lon_min=-180., lon_max=180.,
                 check_raa=False, fix_angs=True):
        # Solar zenith angles
        self.sza = check_type(sza, "Solar zenith angles")
        # Solar azimuth angles
        self.saa = check_type(saa, "Solar azimuth angles")
        # Viewing zenith angles
        self.vza = check_type(vza, "Viewing zenith angles")
        # Viewing azimuth angles
        self.vaa = check_type(vaa, "Viewing azimuth angles")
        # Latitudes
        self.lats = check_type(lats, "Latitudes")
        # Longitudes
        self.lons = check_type(lons, "Longitudes")

        # Set some limits for the angles and geolocation
        self.zenith_min = zenith_min
        self.zenith_max = zenith_max

        self.azimuth_min = azimuth_min
        self.azimuth_max = azimuth_max

        self.relazi_min = relazi_min
        self.relazi_max = relazi_max

        self.lat_min = lat_min
        self.lat_max = lat_max

        self.lon_min = lon_min
        self.lon_max = lon_max

        # Check arrays are equal size or are single values
        self.check_array_shapes()

        # Run checks to ensure values are within bounds
        self.check_angle_bounds(check_raa=False, fix=fix_angs)

        # Relative azimuth angle (range: 0 -> 180)
        # Can be provided by user or calculated internally
        if raa is not None:
            raa = check_type(raa, "Relative azimuth angles")
            self.raa = check_and_reshape(raa, self.sza.shape)
        else:
            self.raa = self.calc_relazi()

        # Check again with RAA
        if check_raa:
            self.check_angle_bounds(check_raa=True, fix=fix_angs)

        # Compute additional parameters
        self.cos_sza = np.cos(np.deg2rad(self.sza))
        self.cos_vza = np.cos(np.deg2rad(self.vza))

        self.sin_sza = np.sin(np.deg2rad(self.sza))
        self.sin_vza = np.sin(np.deg2rad(self.vza))

        self.cos_raa = np.cos(np.deg2rad(self.raa))
        self.sin_raa = np.sin(np.deg2rad(self.raa))

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
            self.saa = check_and_reshape(self.saa, good_shp)
            self.vza = check_and_reshape(self.vza, good_shp)
            self.vaa = check_and_reshape(self.vaa, good_shp)
            self.lats = check_and_reshape(self.lats, good_shp)
            self.lons = check_and_reshape(self.lons, good_shp)

    def check_angle_bounds(self, fix=True, check_raa=True):
        """Ensure that zenith angles are in acceptable range.
        The optional fix argument allows user to explicitly set whether to
        attempt to fix bad angles (by normalising their range) or to raise an error."""
        # Solar zenith
        self.sza = np.where(np.isinf(self.sza), np.nan, self.sza)
        if np.any(self.sza < self.zenith_min) or np.any(self.sza >= self.zenith_max):
            if fix:
                warnings.warn("Some solar zenith values out of range. Clipping.")
                self.sza = np.where(self.sza < self.zenith_min, self.zenith_min, self.sza)
                self.sza = np.where(self.sza > self.zenith_max, self.zenith_max, self.sza)
            else:
                raise ValueError("All solar zenith angles must be in the range " +
                                 str(self.zenith_min) + ' to ' + str(self.zenith_max))
        # Viewing zenith
        self.vza = np.where(np.isinf(self.vza), np.nan, self.vza)
        if np.any(self.vza < self.zenith_min) or np.any(self.vza >= self.zenith_max):
            if fix:
                warnings.warn("Some satellite zenith values out of range. Clipping.")
                self.vza = np.where(self.vza < self.zenith_min, self.zenith_min, self.vza)
                self.vza = np.where(self.vza > self.zenith_max, self.zenith_max, self.vza)
            else:
                raise ValueError("All viewing zenith angles must be in the range " +
                                 str(self.zenith_min) + ' to ' + str(self.zenith_max))
        # Solar azimuth
        self.saa = np.where(np.isinf(self.saa), np.nan, self.saa)
        if np.any(self.saa < self.azimuth_min) or np.any(self.saa >= self.azimuth_max):
            if fix:
                warnings.warn("Some solar azimuth values out of range. Scaling.")
                self.saa = np.where(self.saa > self.azimuth_max, self.saa - self.azimuth_max, self.saa)
                self.saa = np.where(self.saa < self.azimuth_min, self.azimuth_max + self.saa, self.saa)
            else:
                raise ValueError("All solar azimuth angles must be in the range " +
                                 str(self.azimuth_min) + ' to ' + str(self.azimuth_max))
        # Viewing azimuth
        self.vaa = np.where(np.isinf(self.vaa), np.nan, self.vaa)
        if np.any(self.vaa < self.azimuth_min) or np.any(self.vaa >= self.azimuth_max):
            if fix:
                warnings.warn("Some satellite azimuth values out of range. Scaling.")
                self.vaa = np.where(self.vaa > self.azimuth_max, self.vaa - self.azimuth_max, self.vaa)
                self.vaa = np.where(self.vaa < self.azimuth_min, self.azimuth_max + self.vaa, self.vaa)
            else:
                raise ValueError("All viewing azimuth angles must be in the range " +
                                 str(self.azimuth_min) + ' to ' + str(self.azimuth_max))
        # Latitudes
        self.lats = np.where(np.isinf(self.lats), np.nan, self.lats)
        if np.any(self.lats < self.lat_min) or np.any(self.lats >= self.lat_max):
            raise ValueError("All latitudes must be in the range " +
                             str(self.lat_min) + ' to ' + str(self.lat_max))
        # Longitudes
        self.lons = np.where(np.isinf(self.lons), np.nan, self.lons)
        if np.any(self.lons < self.lon_min) or np.any(self.lons >= self.lon_max):
            raise ValueError("All longitudes must be in the range " +
                             str(self.lon_min) + ' to ' + str(self.lon_max))
        # Relative azimuth
        if check_raa:
            self.raa = np.where(np.isinf(self.raa), np.nan, self.raa)
            if np.any(self.raa < self.relazi_min) or np.any(self.raa >= self.relazi_max):
                raise ValueError("All Relative azimuth angles must be in the range " +
                                 str(self.azimuth_min) + ' to ' + str(self.relazi_max))

    def calc_relazi(self):
        """Compute the relative azimuth angles from solar and viewing azimuths."""
        raa = np.abs(self.vaa - self.saa)

        # If any RAA values are greater than 180, invert them.
        raa = np.where(raa >= 180., 360. - raa, raa)
        return raa
