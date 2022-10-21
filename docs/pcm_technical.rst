.. _PCM_Technical:
Details of the pycoxmunk algorithm
=================

Introduction
------------
`pycoxmunk` uses a  modified version of the mathematical equations described by Cox and Munk in their two 1954 papers
describing the `roughness <https://doi.org/10.1364/JOSA.44.000838>`_
and `statistics <https://images.peabody.yale.edu/publications/jmr/jmr13-02-04.pdf>`_ of the sea surface derived
from photographs of sunglint. A more detailed description of the method as-applied to satellite imagery is given in
`Sayer et al, 2010 <https://doi.org/10.5194/amt-3-813-2010>`_

The algorithm described here computes both the per-wavelength sea surface reflectance, :math:`\rho`, and the four
bidirectional reflectance terms:
 - :math:`\rho_{0v}`: Solar beam to satellite view reflectances
 - :math:`\rho_{0d}`: Solar beam to diffuse reflectances
 - :math:`\rho_{dv}`: Diffuse to satellite view reflectances
 - :math:`\rho_{dd}`: Diffuse to diffuse reflectances

To estimate the reflectances, `pycoxmunk` requires knowledge of the sun-satellite viewing geometry through the following
angles:
- :math:`\theta_v`: The satellite zenith angle, defined as 0° when the satellite is directly overhead and 90° when the
  satellite is at the horizon.

- :math:`\phi_v`: The satellite azimuth angle, defined as 0° when the satellite is directly North of the pixel, and 90°
  when the satellite is directly East, 180° when the satellite is directly South and 270° when the satellite is directly
  West.

- :math:`\theta_s`: The solar zenith angle, defined as 0° when the sun is directly overhead and 90° when the sun is at
  the horizon.

- :math:`\phi_s`: The satellite azimuth angle, defined as 0° when the sun is directly North of the pixel, and 90° when
  the sun is directly East, 180° when the sun is directly South and 270° when the sun is directly West.

In addition, satellite sensors typically produce images for a series of `channels`, where each channel has a defined
range of wavelengths over which it is sensitive. This is typically called the `spectral response function` (SRF) of a
sensor.

At present, `pycoxmunk` does not take into account this SRF. Instead, the `central wavelength` for each channel is
used. The spectral properties of water at various wavelengths are defined within `pycoxmunk` in the `CM_Constants.py`
file, and for a given central wavelength the two closest wavelengths in the spectral library are interpolated to
provide the wavelength-specific water properties for each channel. If a channel is selected that is outside the range
of wavelengths supported by the library (~450nm - 3700nm) then the closest library wavelength will be used instead.

Finally, in order to compute the surface properties a source of near-surface wind information is needed. This could be
from observations, numerical weather prediction models or a reanalysis such as the ECMWF/Copernicus `ERA5` products.
The 10 meter wind speeds in the ERA5 dataset have been found to produce best results.


Calculation of water properties
_______________________________

Before performing the reflectance calculations we first compute the properties of the water being observed. Currently,
this is done using fixed values from the literature (described in Sayer, 2010) but in future we plan to implement
support for `ESA's Ocean Color CCI <https://climate.esa.int/en/projects/ocean-colour/>`_ data, which will enable more
a dynamic computation of water properties.

Chlorophyll-A absorption is computed via:
.. math::

    {chl_{abs}} = coef_0 * ( 1 - e^{-1.61 * chl_{conc}})