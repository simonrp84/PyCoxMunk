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

from PyCoxMunk.Utils import CM_Constants as CM_C
import unittest


class TestWaterClass(unittest.TestCase):
    def test_exceptions(self):
        # Check single value is raise das exception
        with self.assertRaises(TypeError):
            test_cls = CM_C.WaterData(None, None, None, None, None, None, 1.0, None, None)
        # Check single item list is raised as exception
        with self.assertRaises(TypeError):
            test_cls = CM_C.WaterData(None, None, None, None, None, None, [1.0], None, None)
        # Check three item list is raised as exception
        with self.assertRaises(TypeError):
            test_cls = CM_C.WaterData(None, None, None, None, None, None, [1.0, 1.0, 1.0], None, None)


if __name__ == '__main__':
    unittest.main()
