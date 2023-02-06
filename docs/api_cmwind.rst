.. _api_cmwind:

The CM_Shared_Wind Module
=========================

The `CM_Shared_Wind` module defines the `CMSharedWind` class that is used py `pycoxmunk` for storing and computing
data associated with wind speed and its effect upon sea surface reflectance. The `CMSharedWind` class contains numerous
variables and helper functions, which are described here. In general, the user should only have need of initialising
the class itself, with the helper functions being used internally by `pycoxmunk`.

The `CMSharedWind` class
------------------------

Initialisiation
^^^^^^^^^^^^^^^

This class is initialised by the user with three arguments: Firstly, the scene geometry must be supplied via the
`CMSceneGeom <./api_cmscenegeom.html>`_ class. The user must also supply `u10` and `v10`, which define the U and V
components of wind speed at 10m above the surface. These data are commonly output from numerical weather models such as
NOAA and ECMWF's operational forecasts, and are also available via widely-used reanalysis datasets such as ERA5. These
two wind speeds can be provided each as a single float, which will result in a constant wind speed being used across the
satellite image, or as arrays of the same shape as the satellite image, in which case the wind is assumed to be
per-pixel.

The class also contains numerous other variables, but these are calculated internally and are not required to be
provided by the user.

The simplest way to initialise `CMSharedWind` is to use the helper function from the main `pycoxmunk.PyCoxMunk` class:

.. code-block:: python

    from pycoxmunk import PyCoxMunk

    # Initialise pycoxmunk
    pcm = PyCoxMunk(scn, bnames)  # Where scn and bnames are defined elsewhere in your code

    # Some function to retrieve wind
    wind_u10, wind_v10 = get_wind()

    # Initialise winds
    pcm.setup_wind(wind_u10, wind_v10 )


You can also initialise the class manually:

.. code-block:: python

    from pycoxmunk.CM_SceneGeom import CMSceneGeom
    from pycoxmunk.CM_Shared_Wind import CMSharedWind

    # Some function to retrieve wind
    wind_u10, wind_v10 = get_wind()

    # Some function to set up the scene geometry
    scn_geom = CMSceneGeom()

    # Initialise winds
    my_winds = CMSharedWind(scn_geom, wind_u10, wind_v10 )


The `zeisse_ba` function
^^^^^^^^^^^^^^^^^^^^^^^^

The `zeisse_ba` function is used internally for computing a correction factor to the sea surface reflectance that is
used at zenith angles greater than 70°. This process is described in
`Zeisse's 1995 paper <https://doi.org/10.1364/JOSAA.12.002022>`_.

Arguments:
 - `theta`: Float or array, the zenith angle in radians.
 - `cos_theta`: Float or array, the cosine of the zenith angle.
 - `ws`: Float or array, the wind speed in m/s.

Returns:
 - `ba`: Float or array (depending on input variables), the ergodic cap area when the zenith angle is greater than 70°
   and wind speed is greater than 1m/s. Otherwise, the cosine of the zenith angle.

The `calc_wind` function
^^^^^^^^^^^^^^^^^^^^^^^^

A helper function to calculate variables required by the `pycoxmunk` code that are not directly provided by the user.
This function is called internally when the `CMSharedWind` class is initialised.
Argument:
- `scenegeom`: A `CMSceneGeom` type that contains the various satellite and solar geometry variables.

Returns:
 - Nothing

Nothing is returned by this function, but various internal variables relating to sea surface roughness and reflectance
are calculated and stored within the `CMSharedWind` class.