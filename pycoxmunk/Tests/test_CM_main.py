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
from pycoxmunk.CM_Main import PyCoxMunk
from unittest import mock
from satpy import Scene
import unittest


class TestCMCalcs(unittest.TestCase):
    """Tests for the overarching CM_Main module."""

    def test_init(self):
        """Test creation of the CM class."""

        # Check we get an error if we don't pass a bool to various options.
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, do_brdf='potato')
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, mask_bad='potato')
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, delete_when_done='potato')
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None)
        with self.assertRaises(ValueError):
            PyCoxMunk(None, band_names='potato')
        with self.assertRaises(ValueError):
            PyCoxMunk(None, band_names=[])
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, angle_names=1.0)
        with self.assertRaises(ValueError):
            PyCoxMunk(None, None, oc_dir=1.0)
