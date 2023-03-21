# !/usr/bin/env python
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
"""Utility functions used in multiple places in the code."""
import dask.array as da_ar
import xarray as xr
import numpy as np


def _gauss_leg_point(n, i, ia, ib, da):
    """Compute a given Gauss-Legendre point."""

    db = (4 * i + 3) / ia * np.pi
    y1 = da_ar.cos(db + da / da_ar.tan(db))

    p2 = 0

    for looper in range(0, 100000):
        p2 = 0
        p1 = 1
        for j in range(0, int(da_ar.floor(n))):
            p3 = p2
            p2 = p1

            db = y1 * p2
            p1 = j / (j + 1) * (db - p3) + db

        db = 1 - y1 * y1
        p1p = n * (p2 - y1 * p1) / db

        p1pp = (2 * y1 * p1p - ib * p1) / db

        y2 = y1
        db = p1 / p1p
        y1 = y2 - db * (1 + db * p1pp / (2 * p1p))

        if da_ar.abs(y1 - y2) < 3e-14:
            break

    return y1, p2


def gauss_leg_quadx(n, x1, x2):
    """Calculate the abscissas and weights of the Gauss-Legendre n-point quadrature formula.

    Uses the method described in:
    Philip J. Davis and Philip Rabinowitz. Methods of Numerical Integration. Dover
    Publications Inc., 31 East 2nd Street, Mineola, Ney York 11501, second edition,
    1984. ISBN 0486453391.

    Inputs:
     - n: Int, the number of points to compute.
     - x1: Float, lower limit for computation.
     - x2: Float, upper limit for computation.
    Returns:
     - x: array, the abscissas.
     - w: array, the weights.
    """

    if n <= 0:
        raise ValueError("Gauss-Leg number of points must be greater than zero.")

    x = da_ar.zeros(int(np.floor(n)))
    w = da_ar.zeros(int(np.floor(n)))

    nn = (n + 1) / 2.

    ia = 4 * n + 2
    ib = n * (n + 1)

    da = (n - 1) / (8. * n ** 3)

    db = (x1 + x2) / 2.
    dc = (x2 - x1) / 2.

    ii = int(da_ar.floor(n - 1))
    for i in range(0, int(da_ar.floor(nn))):
        y1, p2 = _gauss_leg_point(n, i, ia, ib, da)

        dd = dc * y1

        x[i] = db - dd
        x[ii] = db + dd

        de = n * p2
        de = 2 / (de * de)
        w[i] = (dc - dd * y1) * de

        w[ii] = w[i]

        ii = ii - 1

    return x, w


def check_and_reshape(arr, good_shape):
    """If array is single value then scale to match other arrays."""
    if type(arr) is float:
        return arr
    if arr.shape == (1,):
        arr = da_ar.full(good_shape, arr[0])
    elif arr.shape != good_shape:
        raise ValueError("Cannot resize array and sizes do not match.")
    return arr


def check_type(in_val, var_typ):
    """Check that input variable is correct type.
    All inputs should be numpy array or a float."""
    if type(in_val) == float or type(in_val) == np.float64 or type(in_val) == float:
        return in_val
    elif isinstance(in_val, np.ndarray):
        return da_ar.array(in_val)
    elif isinstance(in_val, xr.DataArray):
        return da_ar.array(in_val)
    elif isinstance(in_val, da_ar.Array):
        return in_val
    else:
        raise TypeError(f'{var_typ} must be a single float or numpy array! Got: {type(in_val)}')


def _write_gdal(fname, datas):  # pragma: no cover
    from osgeo import gdal
    driver = gdal.GetDriverByName("GTiff")
    shp = datas.shape
    dst_ds = driver.Create(fname,
                           shp[1],
                           shp[0],
                           1,
                           gdal.GDT_Float32)
    dst_ds.GetRasterBand(1).WriteArray(datas)
    del dst_ds
