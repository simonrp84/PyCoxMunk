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
"""Test the constants module."""

from pycoxmunk import CM_Constants as CM_C
import pytest


class TestWaterClass:

    def test_exceptions(self):
        """Check that passing wrong types to the WaterData class raises an exception.

        Correct type for the chl_a_coef value is a 2-element list"""
        # Single float
        with pytest.raises(TypeError):
            CM_C.WaterData(None, None, None, None, None, None, 1.0, None, None)
        # Single item list
        with pytest.raises(TypeError):
            CM_C.WaterData(None, None, None, None, None, None, [1.0], None, None)
        # Three item list
        with pytest.raises(TypeError):
            CM_C.WaterData(None, None, None, None, None, None, [1.0, 1.0, 1.0], None, None)

    def test_types(self):
        # Check constants are floats."""
        assert type(CM_C.chl_a_conc) == float
        assert type(CM_C.n_air) == float

    def test_waterdata(self):
        """Check WaterData class works correctly."""
        test_data = CM_C.WaterData(0.5, 1.0, 0.02, 0.1, 0.2, 0.8, [0.1, 0.2], 0.4, 0.6)
        assert isinstance(test_data, CM_C.WaterData)
        assert test_data.wavelength == 0.5
        assert test_data.total_backscat == 0.6
        assert test_data.chl_a_coef[0] == 0.1

    def test_default_channel_data(self):
        """Check that default water data is correct."""
        assert CM_C.CM_DATA_DICT[0.47].chl_a_coef[0] == 3.460e-02
        assert CM_C.CM_DATA_DICT[1.6].refrac_real == 1.323e+00
        assert CM_C.CM_DATA_DICT[3.7].total_abs == 1.223e+04
        assert CM_C.CM_DATA_DICT[0.55].wavelength == 0.55
        assert CM_C.CM_DATA_DICT[0.87].base_abs == 5.365e+00
