.. _api_cmutils:

CM_Utils Functions
==================

The `pycoxmunk.CM_Utils` module contains utility functions used elsewhere in the code. In general, these should not
need to be used by the user but they are documented here:


The `gauss_leg_quadx` function
------------------------------

This function computes the abscissas and weights of the Gauss-Legendre n-point quadrature. This utility function is
included primarily for use within the Cox-Munk BRDF calculations, for computing the diffuse terms, but is available for
users if needed. It is called with three arguments:

 - `n`: Integer, the number of points to compute. Higher `n` is more accurate in downstream computations, but slower.
 - `x1`: Float, the lower limit for computation. Often `0`.
 - `x2`: Float, the upper limit for computation. Often `π/2`.

Two variables are returned:

 - `x`: An array of length `n` giving the abscissa values.
 - `w`: An array of length `n` giving the weights.

The `check_and_reshape` function
--------------------------------

This utility function checks that a given variable is either a float or an array of the same shape as another,
reference, array. It is used during processing to ensure that arrays of different sizes are not passed around the code,
it a singular array is given as input, it will be converted to the same shape as the reference. It is called with two
arguments:

 - `arr`: An input variable of any type.
 - `good_shape`: Tuple giving the shape of a reference array.

The function returns:
 - A `float` if `arr` is a float
 - An array of shape `good_shape` with all values equalling `arr` if `arr` is a singular array.
 - An unmodified copy of `arr` if it is an array of the same shape as `good_shape`.

A `ValueError` is raised if `arr` is neither a float nor an array, or if `arr` does not have a compatible shape.

The `check_type` function
-------------------------

This function ensures that a given input variable is of the correct type. It expects the input to be: Any float, a numpy
`ndarray`, an xarray `DataArray` or a dask array. Numpy and xarray types are converted into dask arrays. Any other type
of input raises a `TypeError`.

Inputs:
 - `in_val`: The input variable of any type.
 - `var_typ`: The name of the input variable, `solar_zenith_angle` for example. Used only when raising `TypeError`.

Returns:
 - If `in_val` is a float or dask array, an unmodified `in_val` is returned.
 - If `in_val` is a numpy array or xarray, a dask version of `in_val` is returned.

The `_write_gdal` function
--------------------------

This function is intended to be private and is used during testing to save intermediate output as TIFF images. It takes
two arguments, `fname`, a string specifying the filename to save an image to, and `datas`, a 2d array containing the
data to be saved.

The `_gauss_leg_point` function
-------------------------------

This function is not intended to be called by users, it is called internally by `gauss_leg_quadx`.