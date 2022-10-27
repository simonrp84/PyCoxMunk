.. _api_cmscenegeom:

PyCoxMunk Class
===============

The core of this library is the `PyCoxMunk` class, from which all data can be loaded, processed and retrieved. This
class can be initialised with:

.. code-block:: python

    import pycoxmunk
    pcm = pycoxmunk.PyCoxMunk(scn,
                              band_names,
                              oc_dir,
                              angle_names,
                              do_brdf,
                              mask_bad,
                              delete_when_done)

The first two arguments are mandatory:

 - `scn`: A `satpy.Scene` class containing loaded satellite data.
 - `band_names`: A list of string specifying which band / channel names to process.

The other arguments are optional:

 - `oc_dir`: The directory containing ocean color data. This is currently unused and is a placeholder for when ocean
    color support is added at a later date.
 - `angle_names`: This can be a string, 'calc' to request that `pycoxmunk` calculates all required sun and satellite
    angles, which is not possible for all satellites but can typically be done for geostationary platforms.
    Alternatively, it can be a dict with the following keys: 'sza', 'saa', 'vza', 'vaa'. In this case, the key values
    should be strings specifying the solar zenith and azimuth and satellite zenith and azimuth datasets respectively,
    giving their name as loaded via satpy. The default option is `calc`. See also note below.
 - `do_brdf`: A boolean value to set whether `pycoxmunk` calculates bidirectional reflectance values in addition to the
    reflectance along the sun-satellite path. Calculating these data is of value to some applications, such as radiative
    transfer, but substantially increased runtime and compute (CPU + memory) requirements. By default it is `False`.
 - `mask_bad`: A boolean to specify whether unrealistic reflectance values (such as those below zero) are masked. By
    default this is `True`.
 - `delete_when_done`: A boolean to specify whether ancillary data produced during processing is deleted once the
    reflectance values have been calculated. `pycoxmunk` produces some extra data - such as water underlight
    reflectance - that is typically not of interest to the user. This can be deleted to save memory.
    Default value is `True`. Data can also be manually deleted by calling the `_run_delete` function within `PyCoxMunk`.

NOTE: Zenith and azimuth angles must be specified in the correct range! `pycoxmunk` expects zenith angles >0°, with 0°
indicating that the sun or satellite is directly overhead a given pixel. A value of 90° indicates that the sun or
satellite is on the horizon and values of >90° specify that the sun or satellite is below the horizon. For azimuths,
if the sun or satellite is due North then the azimuth should be 0° while due South the azimuth must be 180°. For due
East, the azimuth is 90° and for due West it is 270°. Different satellites, and indeed different datasets for the same
satellite, may use alternative definitions of these angles. It is up to the user to ensure that the angles supplied to
`pycoxmunk` are in the correct format.

Using the PyCoxMunk class
-------------------------

The `PyCoxMunk` class contains various functions that provide the core functionality for this library:

 - `setup_pixmask`: Initialises a pixel-based mask to specify which pixels to process and which should be set to a fill
    value.
 - `setup_wind`: Initialises the wind data used to compute the sea surface reflectance.
 - `retr_coxmunk_refl`: Runs the sea surface reflectance code to retrieve both sun-satellite reflectance and, if
    requested, also the bidirectional reflectance.

In addition, the `_run_delete` function is not designed for direct use, instead the `delete_when_done` argument to
`PyCoxMunk` should be used, but this can be called to delete unneeded ancillary information to reduce memory use.

The `setup_pixmask` function
----------------------------

This function computes a pixel-level mask of which regions should be processed. It takes four arguments, all of which
are optional:

 - `cloud_mask`: An array of equal size to the satellite data, specifying where clouds are located. A value of 0
    indicates cloud-free while 1 or greater indicates cloudy pixels that are not to be processed.
 - `land_mask`: An array of equal size to the satellite data, specifying where land is located. A value of 0
    indicates water while 1 or greater indicates land pixels that are not to be processed.
 - `sol_zen_mask`: An array of equal size to the satellite data, specifying where pixel solar zenith angle exceeds a
    given threshold value.
 - `sat_zen_mask`: An array of equal size to the satellite data, specifying where pixel solar zenith angle exceeds a
    given threshold value.

This function populates the `pixmask` in the `PyCoxMunk` class. More details are given in the `CMPixMask` documentation.

The `setup_wind` function
_________________________

This function adds wind information to the `PyCoxMunk` class and is called with two arguments:

 - `u10`: A single value or array of same shape as the satellite data, specifying the U (East-West) wind speed in
    meters per second. If a single value, this is assumed to be the wind for each pixel in the satellite image. If an
    array, it is expected that each pixel has a defined wind speed.
 - `v10`: A single value or array of same shape as the satellite data, specifying the V (North-South) wind speed. The
    same assumptions apply as for `u10`.

This function populates the `shared_wind` in the `PyCoxMunk` class. More details are given in the `CMSharedWind`
documentation.

The `retr_coxmunk_refl` function
________________________________

This function calls the computation of sea surface reflectance for all input data channels (given via the `band_names`
argument to `PyCoxMunk`. It will also compute the bidirectional reflectance if requested via `do_brdf` in `PyCoxMunk`.
The function takes no arguments and returns no data. Results are stored inside the `PyCoxMunk.scn` class as datasets.
For each channel in `band_names`, the resulting reflectance can be retrieved via:

`PyCoxMunk.scn[f'cox_munk_refl_{band_name}']`

Similarly, the bidirectional components can be retrieved with:

`PyCoxMunk.scn[f'cox_munk_rho0d_{band_name}']`
`PyCoxMunk.scn[f'cox_munk_rho0v_{band_name}']`
`PyCoxMunk.scn[f'cox_munk_rhodv_{band_name}']`
`PyCoxMunk.scn[f'cox_munk_rhodd_{band_name}']`
