.. _api_cmconsts:

The CM_Constants module
=======================

This module stores a variety of constants and other data used during processing. It also defines the `WaterData` class
that stores the wavelength-dependant water properties.

The `WaterData` class
---------------------

This class stores per-wavelength water properties, used within the code to compute - for example - the white cap
reflectance.

Mandatory class attributes:
 - `wavelength`: Float, the wavelength, in micrometers, at which the remaining attributes are defined.
 - `refrac_real`: Float, he real component of the refractive index.
 - `refrac_imag`: Float, he imaginary component of the refractive index.
 - `base_abs`: Float, the base water absorption coefficient.
 - `base_backscat`: Float, the base water backscatter coefficient.
 - `whitecap_refl`: Float, the whitecap reflectance of water at this wavelength.
 - `chl_a_coef`: List of length 2, the Chlorophyll-A absorption coefficient defined as
 - `total_abs`: Float, the total water absorption.
 - `total_backscat`: Float, the total water backscatter.
 - `clhabs`: Float, the Chlorophyll-A absorption.
 - `clhbsc`: Float, the Chlorophyll-A backscatter.

The two values in `chl_a_coef` are used to compute absorption using the following equation:

.. math::
    abs_{chl} = coef_0 * ( 1 - e^{-1.61 * conc_{chl}}) + coef_1 * conc_{chl}

Where :math:`conc_{chl}` is the Chlorophyll concentration in :math:`mg/m^3`.

At present, `pycoxmunk` contains predefined `WaterData` information for the following wavelengths: 0.47, 0.55, 0.65,
0.87, 1.24, 1.375, 1.6, 2.13, and 3.7 :math:`\mu m`

The constants
-------------

Also defined in this module are several constants used throughout the processing:

 - `dither_more`: A float specifying a very small number, used to test if variables are significantly greater than zero.
 - `chl_a_conc`: The default Chlorophyll-A concentration, 0.18 :math:`mg/m^3`. In the future this will be replaced by
   data from ESA's Ocean Color CCI.
 - `n_air`: The refractive index of air, defined as 1.00029
 - `cm_min_wvl`: The minimum wavelength for processing via `pycoxmunk`. Defined as 0.2 :math:`\mu m`
 - `cm_max_wvl`: The maximum wavelength for processing via `pycoxmunk`. Defined as 4.5 :math:`\mu m`
 - `zeisse_coef`: Coefficients used for computing the Zeisse correction, which is applied at high zenith angles.
 - `n_quad_theta`: The number of quadrature points for computing zenith component of the diffuse BRDF.
 - `n_quad_phi`: The number of quadrature points for computing azimuth component of the diffuse BRDF.

The default value of both `n_quad_theta` and `n_quad_phi` is 4. Increasing either or both of these beyond 4 will result
in more accurate calculation of the diffuse BRDF terms but will have a significant impact on processing time and memory
requirements.

Lastly, `CM_DATA_DICT` is a dictionary that stores `WaterData` information at all of the wavelengths specified earlier.
