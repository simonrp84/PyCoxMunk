.. _api_cmscenegeom:

The CM_SceneGeom Module
=======================

The `CM_SceneGeom` module defines the `CMSceneGeom` class that stores solar and satellite geometry information that is
used throughout the `pycoxmunk` code. In a typical processing workflow, it will be initialised during the initialisation
of the main `PyCoxMunk` class and does not need to be initialised directly by the user. The documentation below
describes the class to assist with advanced use cases or use outside `pycoxmunk`'s core processing.


The `CMSceneGeom` class
------------------------

This class houses the solar and satellite angles required by processing, namely the solar zenith and azimuth angles, and
the satellite zenith and azimuth angles. A variety of other geometry information is computed within the module,
either during initialisation or by using the helper functions. In addition, the class checks that input values are
within the expected range, as `pycoxmunk` behaviour is undefined for non-compliant angles.

Notes on angles
^^^^^^^^^^^^^^^

Scene geometry has been defined in many different ways between various sensors, platforms and codebases. In `pycoxmunk`
the following definitions are used:

Zenith angles:
 - Should be 0° for the case when the sun/satellite is directly above the pixel.
 - Should be 90° for the case when the sun/satellite is on the horizon.
 - Should be 180° for the extreme case of the sun/satellite being on the opposite side of the Earth.

Azimuth angles:
 - Should be 0° when the sun/satellite is directly North of the pixel.
 - Should be 90° when the sun/satellite is directly East of the pixel.
 - Should be 180° when the sun/satellite is directly South of the pixel.
 - Should be 270° when the sun/satellite is directly West of the pixel.

Relative azimuth angles:
 - Should be 0° when the sun and satellite are aligned on the same side of the pixel (backscatter).
 - Should less than 90° when the sun and satellite are on the same side of the pixel but not aligned (backscatter).
 - Should greater than 90° when the sun/satellite are on opposite sides of the pixel but not aligned (forward scatter).
 - Should be 180° when the sun/satellite are aligned  on opposite sides of the pixel (forward scatter).

Latitudes:
 - Should be between 0° and 90° for pixels North of the Equator
 - Should be between 0° and -90° for pixels South of the Equator

Longitudes:
 - Should be between 0° and 180° for pixels East of the prime meridian (Greenwich).
 - Should be between 0° and -180° for pixels West of the prime meridian (Greenwich).

Initialisation of the class
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Initialisation of the `CMSceneGeom` class requires 6 arguments and can also accept various optional arguments:

Mandatory arguments:

 - `sza`: Float or array, the solar zenith angle in degrees.
 - `saa`: Float or array, the solar azimuth angle in degrees.
 - `vza`: Float or array, the satellite zenith angle in degrees.
 - `vaa`: Float or array, the satellite azimuth angle in degrees.
 - `lats`: Float or array, the pixel latitude in degrees.
 - `lons`: Float or array, the pixel longitude in degrees.

Optional arguments:

 - `raa`: Float or array, the relative azimuth angle in degrees. Calculated internally if not supplied.
 - `zenith_min`: Float, the minimum zenith angle in degrees. Zeniths lower than this will be clipped. Default = 0
 - `zenith_max`: Float, the maximum zenith angle in degrees. Zeniths higher than this will be clipped. Default = 85
 - `azimuth_min`: Float, the minimum azimuth angle in degrees. Azimuths lower than this will be clipped. Default = 0
 - `azimuth_max`: Float, the maximum azimuth angle in degrees. Azimuths higher than this will be clipped. Default = 360
 - `relazi_min`: Float, the minimum rel azimuth angle in degrees. RAAs lower than this will be clipped. Default = 0
 - `relazi_max`: Float, the maximum rel azimuth angle in degrees. RAAs higher than this will be clipped. Default = 180
 - `lat_min`: Float, the maximum latitude in degrees. Latitudes higher than this will be clipped. Default = -90
 - `lat_max`: Float, the maximum latitude in degrees. Latitudes higher than this will be clipped. Default = 90
 - `lon_min`: Float, the maximum longitude in degrees. Longitides higher than this will be clipped. Default = -180
 - `lon_max`: Float, the maximum longitude in degrees. Longitides higher than this will be clipped. Default = 180
 - `check_raa`: Bool, switch whether to check the RAA range against `relazi_min` and `relazi_max`. Default: False
 - `fix_angs`: Bool, switch whether to clip 'bad' angles into acceptable range. Default: True

The `compute_additional` function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To ease the workload on the user, many angles are calculated internally rather than being mandatory arguments. The
`compute_additional` function computes the cosine and sine of the solar and viewing zenith angles, and of the relative
azimuth angle. Results are stored in `self.cos_{val}` and `self.sin_{val}` where val is `sza`, `vza` or `raa` for the
solar zenith, satellite zenith and relative azimuth angles respectively. The function takes no arguments and returns
no variables.

The `check_angle_bounds` function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The `pycoxmunk` code requires angles within a certain range, described above. This function ensures that angles are
within this range and, optionally, clips angles outside the range to the minimum or maximum acceptable value.

The function takes two arguments, both optional:

 - `fix`: Boolean, specifies whether to clip out-of-range angles. Default: True
 - `check_raa`: Boolean, specifies whether to include the relative azimuth angle in checks. Default: True

The function returns no values, but modifies class instance variables if out of range and `fix=True`. When out-of-range
values are found a warning is raised if `fix=True` and a `ValueError` otherwise. The latitude, longitude and relative
azimuth angle checks will always raise an error and are not clipped if out-of-range.

The `calc_relazi` function
^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a helper function to calculate the relative azimuth angle in the case that it is not supplied by the user.
The function calculates:

.. math::
    {\phi_r} = abs(\phi_{sat} - \phi_{sun})

Values greater than 180° are scaled via:

.. math::
    {\phi_r} = 360 - {\phi_r}

The `cm_calcangles` function
----------------------------

This function is separate to the `CMSceneGeom` class, and is used to calculate satellite and solar angles for a given
satellite Scene. This uses the `satpy` functions `_get_sensor_angles` and `_get_sun_angles`.

Arguments:

 - `inscn`: Scene, the input satellite data Scene.
 - `refband`: String, the band name used for computing the angles. Required as multiple resolutions of data may be
   present in the input scene.

Returns:
 - `inscn`: Scene, a modified version of the scene that contains the following angles datasets: `solar_azimuth_angle`,
   `solar_zenith_angle`, `satellite_azimuth_angle`, `satellite_zenith_angle`.

If the `refband` is not present then a `KeyError` will be raised. If the angles cannot be computed, which can happen
for some satellites - although it should work successfully for most, then a `ValueError` is raised.

