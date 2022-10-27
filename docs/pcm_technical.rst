.. _PCM_Technical:
==================================
Details of the pycoxmunk algorithm
==================================

Introduction
============
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
from observations, numerical weather prediction models or a reanalysis such as the
`ECMWF/Copernicus ERA5 <https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5>`_ products. The 10 meter
wind speeds in the ERA5 dataset have been found to produce best results.


Calculation of water properties
===============================
Before performing the reflectance calculations we first compute the properties of the water being observed. Currently,
this is done using fixed values from the literature (described in Sayer, 2010) but in future we plan to implement
support for `ESA's Ocean Color CCI <https://climate.esa.int/en/projects/ocean-colour/>`_ data, which will enable more
a dynamic computation of water properties.

Chlorophyll-A absorption and backscattering are computed via:

.. math::
    {chl_{abs}} = coef_0 \cdot ( 1 - e^{-1.61 \cdot chl_{conc}}) + coef_1 \cdot chl_{conc}
.. math::
    {chl_{bsc}} = \frac{0.02 \cdot ( 0.5 - 0.25 \cdot log(chl_{conc})) \cdot 0.55}{\lambda} + 0.002

where :math:`coef_n` are predefined coefficients extracted from the spectral properties library by wavelength
and :math:`chl_{conc}` is the Chlorophyll concentration, which is currently a fixed value of 0.18.

The total absorption and backscatter are then defined by:

.. math::
    {tot_{abs}} = base_{abs} + chl_{abs}
.. math::
    {tot_{bsc}} = base_{bsc} + chl_{bsc}

Where :math:`base_{abs}` and :math:`base_{bsc}` are extracted from the spectral properties library by wavelength.

The white cap fraction is defined by:

.. math::
    {wc_{frac}} = 2.951x10^{-6} \cdot v_{wind}^{3.52}


Calculation of water reflectance
================================

Now, the water body reflectance is calculated via:

.. math::
    \eta_{oc} = 0.5 \cdot \frac{base_{bsc}}{tot_{bsc}}
.. math::
    f = 0.6279 - (0.2227 \cdot \eta_{oc}) - (0.00513 \cdot \eta_{oc}) + (0.2465 \cdot \eta_{oc} - 0.3119) \cdot cos(\theta_s)
.. math::
    \rho_{water} = f \cdot \frac{tot_{bsc}}{tot_{abs}}

The underlight contribution to total water reflectance is given by:

.. math::
    \rho_{ul} = \frac{t_u + t_d + \rho_{water}}{1 - r_u \cdot \rho_{water}}

where :math:`t_u` is the upward transmission, defined as 0.52, :math:`t_d` is computed from the scene geometry (see
Sayer, 2010 for more details) and :math:`r_u` is the upward reflectance, defined as 0.48.

The sunglint reflectance is calculated by:

.. math::
    \rho_{gl} = \frac{\pi \cdot r_{sf} \cdot P_{slo}}{\beta}


Where:

    :math:`P{slo}` is the wave slope distribution

    :math:`r_{sf}` is the Fresnel reflection coefficient and

    :math:`\beta` is the facet tilt defined by:

.. math::
    cos(\beta) = \frac{cos(\theta_s) + cos(\theta_v)}{\sqrt{2 + 2 \cdot cos(2\cdot\Theta)}}


with :math:`\Theta` defined via:

.. math::
    cos(2\cdot\Theta) = cos(\theta_V)cos(\theta_s) + sin(\theta_v)sin(\theta_s)cos(\phi_r)`

Finally, the total reflectance is calculated using:

.. math::
    \rho = \rho_{wc} + (1 - wc_{frac}) \cdot (\rho_{gl} + \rho_{ul})



Calculation of bidirectional reflectance
========================================
In addition to calculating the reflectance along the sun-surface-satellite path, `pycoxmunk` can also calculate the
bidirectional reflectance terms. These comprise of the:

- :math:`\rho_{0v}`: Solar beam to satellite view reflectances

- :math:`\rho_{0d}`: Solar beam to diffuse reflectances

- :math:`\rho_{dv}`: Diffuse to satellite view reflectances

- :math:`\rho_{dd}`: Diffuse to diffuse reflectances

Calculation of these terms is enabled optionally by the user at runtime and is disabled by default. The BRDF terms
require significant additional computation and hence extend both the processing time and memory requirements
for `pycoxmunk`.

Internally, these terms are calculated by simulating the surface reflectance across a range of angles and summing the
results to gain the direct and diffuse terms listed above. The surface reflectance is simulated using the terms
described in the preceeding sections.

The angles used in the calculation are computed using the
`Gauss-Legendre n-point quadrature formula <https://doi.org/10.1016/C2013-0-10566-1>`_ and by default `pycoxmunk` uses
four terms for both :math:`\theta` and :math:`\phi`. This can be altered within the `pycoxmunk` code by editing
the ``n_quad_theta`` and ``n_quad_phi`` terms in ``CM_Constants.py``. In the future these will be user-configurable options.
Note, however, that increasing either of these values will significantly increase runtime and memory use.

The :math:`\rho_{0v}` term is simply the standard Cox-Munk sea surface reflectance, so no additional calculations are
performed for this term.

The :math:`\rho_{0d}` and :math:`\rho_{dv}` terms use the actual solar or viewing geometry (respectively) and then
simulate diffuse radiation by substituting the viewing or solar geometry with the Gauss-Legendre terms described above.
The following python code offers a simplified workflow for computing these values:

.. code-block:: python

    def get_rho_od_dv():
        qx_qw_sincos = np.cos(gauss_leg_theta_abscissas) * np.sin(gauss_leg_theta_abscissas) * gauss_leg_theta_weights

        # Loop over zeniths
        for i in range(0, n_quad_theta):
            cur_gl_zen = gauss_leg_theta_abscissas[i]
            tmp_0d = 0
            tmp_dv = 0
            # Loop over azimuths
            for j in range(0, n_quad_phi):
                cur_gl_azi = gauss_leg_phi_abscissas[j]
                cur_gl_phi_w = gauss_leg_phi_weights[j]
                # Compute the reflectances
                tmp_0d = tmp_0d + calc_cox_munk_refl(sol_zen, sol_az, cur_gl_zen, cur_gl_azi) * cur_gl_phi_w
                tmp_dv = tmp_dv + calc_cox_munk_refl(cur_gl_zen, cur_gl_azi, sat_zen, sat_azi) * cur_gl_phi_w

            rho_0d = rho_0d + tmp_0d * qx_qw_sincos[i]
            rho_dv = rho_dv + tmp_dv * qx_qw_sincos[i]

        rho_0d = rho_0d / np.pi
        rho_dv = rho_dv / np.pi



Finally, the :math:`\rho_{dd}` term is calculated in a similar way, but substitutes both the solar and viewing zenith
angles as opposed to just one or the other:

.. code-block:: python

    def get_rho_dd():
        qx_qw_sincos = np.cos(gauss_leg_theta_abscissas) * np.sin(gauss_leg_theta_abscissas) * gauss_leg_theta_weights

        # Initialise direct_diffuse reflectance
        rho_dd = 0

        # Loop over solar zeniths
        for i in range(0, n_quad_theta):
            outer_tmp_dd = 0
            cur_zen1 = gauss_leg_theta_abscissas[i]
            # Loop over satellite zeniths
            for j in range(0, n_quad_theta):
                inner_tmp_dd = 0
                cur_zen2 = gauss_leg_theta_abscissas[j]
                # Loop over relative azimuths
                for k in range(0, n_quad_phi):
                    cur_relazi = gauss_leg_phi_abscissas[k]
                    cur_phi_w = gauss_leg_phi_weights[k]
                    # Compute the reflectances
                    inner_tmp_dd = inner_tmp_dd + calc_cox_munk_refl(cur_zen1, .0, cur_zen2, cur_relazi) * cur_phi_w

                outer_tmp_dd = outer_tmp_dd + inner_tmp_dd * qx_qw_sincos[j]
            rho_dd = rho_dd + outer_tmp_dd * qx_qw_sincos[i]

        rho_dd = rho_dd / np.pi