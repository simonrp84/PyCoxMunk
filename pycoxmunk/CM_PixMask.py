#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2022 Simon R Proud
#
# This file is part of PyCoxMunk.
#
# PyCoxMunk is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, version 3.
#
# PyCoxMunk is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PyCoxMunk.  If not, see <http://www.gnu.org/licenses/>.
"""Class for the pixel mask information."""

from pycoxmunk.CM_Utils import check_type
import dask.array as da
import numpy as np


class CMPixMask:
    """A class for data about the pixel masking.

    This can combine multiple masks into one final pixel mask that
    is used in the Cox Munk calculations. Any masked pixel will be
    set to `np.nan` when the mask is applied.

    At present, the following masks are supported:
    - cloud_mask: A mask specifying cloudy and non-cloudy pixels.
                    Cloudy pixels are set to fill.
    - land_mask: A mask specifying land and sea pixels.
                    Land pixels are set to fill
    - sol_zen_mask: A mask to cut off pixels with high solar zenith angles.
                    High zenith pixels (`1`) are cut off.
    - sat_zen_mask: A mask to cut off pixels with high sat zenith angles.
                    High zenith pixels (`1`) are cut off.

    A helper function is available for creating solar and satellite zenith masks:
    cut_high_zen(zeniths, threshold=80)

    This sets up a suitable mask for all pixels greater than threshold.
    Data can be in degrees or radians, but both arguments must be consistent
    and the default for the threshold is in degrees.

    In the case of all masks, pixels to be processed should have a
    value of 0 and pixels to be masked should have a value of >=1.
    """
    def __init__(self, cloud_mask=None, land_mask=None, sol_zen_mask=None, sat_zen_mask=None):

        self.mask = None
        self._check_and_add_mask(cloud_mask, "Cloud mask")
        self._check_and_add_mask(land_mask, "Land mask")
        self._check_and_add_mask(sol_zen_mask, "Solar zenith mask")
        self._check_and_add_mask(sat_zen_mask, "Satellite zenith mask")
        if self.mask is not None:
            self.mask = da.where(self.mask > 0, 1, 0)

    def _check_and_add_mask(self, mask, mtype):
        if mask is not None:
            mask = check_type(mask, mtype)
            if self.mask is None:
                self.mask = mask
            else:
                self.mask = self.mask + mask

    def cut_high_zen(self, zeniths, threshold=80):
        """Create a mask for high solar or satellite zenith angles."""

        # We take the absolute value here to simplify things.
        # Some datasets have zeniths in -90 to 90 range, others in 0 - 90+
        tmp_mask = da.where(np.abs(zeniths) > threshold, 1, 0)
        if self.mask is None:
            self.mask = tmp_mask
        else:
            self.mask = self.mask + tmp_mask
        if self.mask is not None:
            self.mask = da.where(self.mask > 0, 1, 0)
