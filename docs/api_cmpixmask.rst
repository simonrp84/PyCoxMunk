.. _api_cmpixmask:

The CM_PixMask Module
=====================

The `CM_PixMask` module defines the `CMPixMask` class, which implements a pixel-based masking scheme for specifying
which regions of a given satellite image to process.

The CMPixMask class
-------------------

The `CMPixMask` can be initialised with four masks, each of which should use a value of `0` to indicate a pixel to be
processed and a value of greater than zero for a pixel to mask. The four masks are nominally defined as:

 - `cloud_mask`: An array containing pixels to be masked due to cloud contamination.
 - `land_mask`: An array containing pixels to be masked as they are not water.
 - `sol_zen_mask`: An array containing pixels to be masked due to their high solar zenith angle.
 - `sat_zen_mask`: An array containing pixels to be masked due to their high satellite zenith angle.

These four masks are optional, and any - or none - can be specified. The core of `CMPixMask`, once initialised, is the
`CMPixMask.mask` array, which contains a combined mask utilising all of the other masks, if supplied. If no masks are
supplied, all pixels are processed.

The `cut_high_zen` helper function
----------------------------------

The `CMPixMask` class also contains a helper function, `cut_high_zen` that can compute a mask based on a user-supplied
zenith angle array and threshold value:

.. code-block:: python

    from pycoxmunk.CM_PixMask import CMPixMask
    cm_mask = CMPixMask()
    my_solar_zeniths = load_solzens()  # where this is some function to load zenith angles, if not done via satpy
    cm_mask.cut_high_zen(my_solar_zeniths, threshold=70)

The above code will set all pixels in `cm_mask.mask` to `1` where `my_solar_zeniths` is greater than 70. Be aware of
confusion between degrees and radians, the function is agnostic but requires both to be in the same units. Typically,
zenith (and azimuth) angles for satellites are specified in degrees, but this is not always the case.