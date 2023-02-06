.. _PCM_Quickstart:

==========
Quickstart
==========

Pycoxmunk operates in tandem with the  `satpy library <https://github.com/pytroll/satpy>`_, and it is recommended that
users new to working with satellite data first familiarise themselves with the basic operation of `satpy` before
beginning to use `pycoxmunk`.

This page gives a brief overview of how to use `pycoxmunk`. Further information can be found in the
`API documentation <./pcm_api.html>`_ and in the example code / notebooks that are part of the `pycoxmunk`
`github repository <https://github.com/simonrp84/PyCoxMunk/tree/main/Examples>`_.

Initialising pycoxmunk
======================

Before using _pycoxmunk_, you must load satellite data with _satpy_ or your own code that emulates satpy xarray
`datasets`:

.. code-block:: python

    # Load libraries
    from satpy import Scene
    from glob import glob

    # Set input dir and find files
    my_input_dir = '/some_data/directory/'
    my_files = glob(f'{my_input_dir}/*.nc')

    # Create Scene
    scn = Scene(my_files, reader='my_reader')  # Where 'my_reader' is a satpy reader such as 'abi_l1b'.

    # Which bands do you want to load
    my_bands = ['my_band1', 'my_band2'] # For 'abi_l1b' this could be, for example, ['C01', 'C04']

    # Load the data
    scn.load(my_bands)

We can now load _pycoxmunk_ and set up our processing

.. code-block:: python

    # Load library
    import pycoxmunk

    # Initialise the main class
    pcm = pycoxmunk.PyCoxMunk(scn, my_bands)


This example assumes the default options, so please view the `API documentation <./pcm_api.html>`_ for more advanced usage.
_pycoxmunk_ requires geometry information: The sun and satellite zenith and azimuth angles relative to the pixels in the
data. These can be supplied via a dict or - for some satellites - can be calculated directly.

Next, we can (optionally) load wind and pixel masks. The wind data should be on the same grid as the satellite data, so
may require resampling if it is, for example, output from an NWP model. Alternatively, you can use `float` values for
the wind speeds, which will set constant wind over the entire scene. If you do not set wind speeds then _pycoxmunk_ will
assume zero wind.

Pixel masks can mask out areas that you do not wish to process, such as land pixels or pixels filled with cloud. These
are optional and by default _pycoxmunk_ will process all pixels in an image. Each mask should contain zeros for 'clear'
regions where you want to retrieve Cox-Munk reflectances and be greater than 0 for all regions where you wish to skip
processing.

.. code-block:: python

    my_u10_wind, my_u10_wind = load_some_wind() # Resulting in numpy or xarray of same shape as satellite data
    pcm.setup_wind(my_u10_wind, my_u10_wind)

    my_cloudmask = get_cloudmask() # Resulting in array of same shape as satellite data. Sometimes this data is
                                   # available via satpy

    pcm.setup_pixmask(cloud_mask = my_cloudmask)



Retrieving results
==================

Now, finally, we can calculate the reflectances themselves:

.. code-block:: python

    pcm.retr_coxmunk_refl()

The results will be stored as new datasets in the `scn` variable within the PyCoxMunk class, such as
`pcm.scn['cox_munk_refl_my_band']`. These datasets can then either be saved to disk via satpy, or worked on further by
accessing the underlying data arrays.

.. code-block:: python

    # Save the data via satpy
    import numpy as np
    pcm.save_dataset('cox_munk_refl_my_band', enhance=False, dtype=np.float32, fill_value=0, filename='/out_dir/file.tif')

    # Continue processing the data
    my_reflectance_diff = pcm.scn['cox_munk_refl_my_band'] - pcm.scn['my_band']

Using dask
==========

The `dask` library is used by `pycoxmunk`, as this enables large images - such as whole disk geostationary satellite
data - to be processed on machines with limited memory availability. However, dask can require careful tuning to ensure
it works optimally on machines with limited compute. The `complex_seviri.ipynb` example notebook gives an example of
tuning dask to run on a standard desktop machine.

In general, the user should set the maximum number of workers to be used by dask and the chunk size, as follows:

.. code-block:: python

    import dask
    dask.config.set({'array.chunk-size': '32M', 'num_workers': 4})

Some experimentation may be necessary to find the best values for these, and the above work well to process SEVIRI data
on a machine with 16Gb memory and 8 cores. For most processing, this may be unneeded as the `pycoxmunk` code is not
memory intensive. However, the optional BRDF calculations