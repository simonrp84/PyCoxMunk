
[![PyCoxMunk Tests](https://github.com/simonrp84/PyCoxMunk/workflows/Tests/badge.svg)](https://github.com/simonrp84/PyCoxMunk/actions)
[![codecov](https://codecov.io/gh/simonrp84/PyCoxMunk/branch/main/graph/badge.svg?token=4GMLURHA5V)](https://codecov.io/gh/simonrp84/PyCoxMunk)

PyCoxMunk
=========

### A python library for simulating satellite-viewed sea surface reflectances.

The pycoxmunk library computes the sea surface reflectance, using the 
[Cox-Munk method](https://doi.org/10.1364/JOSA.44.000838), as expected to be seen from space based on the satellite and 
weather conditions at the sea surface.

The main purpose of this library is to enable an easy method for computing sea surface reflectance that is applicable
to a wide variety of low-earth and geostationary orbit satellites and sensors. The library is closely linked to the 
[satpy library](https://github.com/pytroll/satpy) that loads, calibrates and produces projection information for many 
satellites. Under the hood, pycoxmunk uses [dask](https://github.com/dask/dask) and 
[xarray](https://github.com/pydata/xarray) to store and process data. Results are made available via the same `satpy` 
scene used to provide the input satellite images.

As well as computing the sea surface reflectance, pycoxmunk can also compute the four BRDF terms, which can be useful in
radiative transfer, cloud and aerosol properties retrieval codes, and other uses such as cloud or sea ice detection.

Installation
------------

pycoxmunk can be installed from PyPI with pip:

```bash
    pip install pycoxmunk
```

It is also available from `conda-forge` for conda installations:

```bash
    conda install -c conda-forge pycoxmunk
```

Credits
-------

This code was written by Simon Proud and is based on a fortran implementation of the Cox-Munk algorithm written by 
Greg McGarragh as part of the [ORAC algorithm](https://github.com/ORAC-CC/orac).

Feedback and Contribution
-------------------------

Feedback, suggestions, bug reports or any other type of contribution is welcome.

If you encounter any problems with this code or the documentation then please file an 
[issue](https://github.com/simonrp84/PyCoxMunk/issues).
It may help in debugging any problems to enable satpy's debug mode:


```python
    from satpy import debug_on
    debug_on()
```

This will print additional log and diagnostic information.

Suggestions for new features are welcome, but may not always be possible for me to code due to limited time. You can
also submit your own [pull requests](https://github.com/simonrp84/PyCoxMunk/pulls) that add features of fix bugs. This is the recommended way to change the library
code, rather than emailing me your updates. By submitting a pull request your changes are documented and your 
contribution to the code is clear.
