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
 - :math:`\rho_0_v`: Solar beam to satellite view reflectances
 - :math:`\rho_0_d`: Solar beam to diffuse reflectances
 - :math:`\rho_d_v`: Diffuse to satellite view reflectances
 - :math:`\rho_d_d`: Diffuse to diffuse reflectances

To estimate the reflectances, `pycoxmunk` requires knowledge of the sun-satellite viewing geometry through the following
angles:
 - :math:`\theta_v`: The satellite zenith angle, defined as 0 :math:`^\circ` when the satellite is directly
overhead and 90 :math:`^\circ` when the satellite is at the horizon.
 - :math:`\phi_v`: The satellite azimuth angle, defined as 0 :math:`^\circ` when the satellite is directly
North of the pixel, and 90 :math:`^\circ` when the satellite is directly East, 180 :math:`^\circ` when the satellite is
directly South and 270 :math:`^\circ` when the satellite is directly West.
 - :math:`\theta_s`: The solar zenith angle, defined as 0 :math:`^\circ` when the sun is directly
overhead and 90 :math:`^\circ` when the sun is at the horizon.
 - :math:`\phi_s`: The satellite azimuth angle, defined as 0 :math:`^\circ` when the sun is directly
North of the pixel, and 90 :math:`^\circ` when the sun is directly East, 180 :math:`^\circ` when the sun is
directly South and 270 :math:`^\circ` when the sun is directly West.