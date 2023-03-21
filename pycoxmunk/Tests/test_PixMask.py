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
"""Test the pixel masking module."""

import numpy as np
from pycoxmunk.CM_PixMask import CMPixMask


class TestCMPixMask:
    """Test the various pixel masking methods."""
    def setup_method(self):
        """Set up the initial data."""
        self.zens = np.array([0., 0.01, 12., 15, 45.35, 79.99, 84.23, 91.1, 112.12, -10.])
        self.cloud_mask = np.array([0, 1, 1, 0, 0, 0, 0, 1, 1, 1])
        self.sol_mask = np.array([1, 1, 1, 1, 0, 0, 1, 0, 0, 1])
        self.sat_mask = np.array([0, 0, 1, 0, 1, 0, 0, 0, 0, 1])
        self.land_mask = np.array([1, 0, 0, 0, 1, 1, 1, 1, 1, 0])

        self.test_mask = np.array([1, 1, 1, 1, 0, 0, 1, 1, 1, 1])

        self.test_zenmask = np.array([0, 1, 1, 0, 0, 0, 1, 1, 1, 1])
        self.test_zenmask2 = np.array([0, 0, 0, 0, 0, 1, 1, 1, 1, 0])

    def test_init(self):
        """Test the masks are combined correctly."""
        # No masks
        main_mask = CMPixMask()
        assert main_mask.mask is None

        main_mask = CMPixMask(cloud_mask=self.cloud_mask)
        np.testing.assert_allclose(main_mask.mask, self.cloud_mask)

        main_mask = CMPixMask(cloud_mask=self.cloud_mask, sol_zen_mask=self.sol_mask)
        np.testing.assert_allclose(main_mask.mask, self.test_mask)

    def test_zencutter(self):
        """Test the utility function for cutting high zeniths."""
        main_mask = CMPixMask(cloud_mask=self.cloud_mask)
        main_mask.cut_high_zen(self.zens)
        np.testing.assert_allclose(main_mask.mask, self.test_zenmask)

        main_mask = CMPixMask()
        main_mask.cut_high_zen(self.zens, threshold=70)
        np.testing.assert_allclose(main_mask.mask, self.test_zenmask2)
