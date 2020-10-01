#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2020 Simon R Proud
#
# This file is part of PyCoxMunk.
#
# PyCoxMunk is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# PyCoxMunk is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# PyCoxMunk.  If not, see <http://www.gnu.org/licenses/>.
"""Some constants used throughout the rest of the code."""


class WaterData:
    """A class for data about water at a given wavelength."""

    def __init__(self,
                 wavelength,       # Wavelength of given information
                 refrac_real,      # Real part of refractive index
                 refrac_imag,      # Imaginary part of refractive index
                 base_abs,         # Base water absorption
                 base_backscat,    # Base water backscatter coefficient
                 whitecap_refl,    # Reflectance of whitecaps
                 chl_a_coef,       # Chlorophyll-a absorption coefficient part
                 total_abs,        # Total water absorption
                 total_backscat):  # Total water backscatter

        # Check the chl-a coefficients are a 2-element list
        if type(chl_a_coef) is not list:
            raise TypeError("Chlorophyll-a absorption coefficient should be a 2-element list.")
        if len(chl_a_coef) != 2:
            raise TypeError("Chlorophyll-a absorption coefficient should be a 2-element list.")
        self.wavelength = wavelength
        self.refrac_real = refrac_real
        self.refrac_imag = refrac_imag
        self.base_abs = base_abs
        self.base_backscat = base_backscat
        self.whitecap_refl = whitecap_refl
        self.chl_a_coef = chl_a_coef
        self.total_abs = total_abs
        self.total_backscat = total_backscat


# Approximate median chlorophyll - A concentration from GlobCOLOUR (mg / m3)
chl_a_conc = 0.18

# Refractive index of air
n_air = 1.00029

# Wavelength dependent coefficients.
# Values are taken and / or computed from the sources indicated in Sayer et al., 2010.
CM_DATA_DICT = {0.47: WaterData(wavelength=0.47,
                                refrac_real=1.345e+00,
                                refrac_imag=9.520e-10,
                                base_abs=1.600e-02,
                                base_backscat=3.780e-03,
                                whitecap_refl=4.408e-01,
                                chl_a_coef=[3.460e-02, 1.240e-02],
                                total_abs=2.694e-02,
                                total_backscat=3.761e-03),
                0.55: WaterData(wavelength=0.55,
                                refrac_real=1.341e+00,
                                refrac_imag=1.960e-09,
                                base_abs=6.400e-02,
                                base_backscat=1.930e-03,
                                whitecap_refl=4.024e-01,
                                chl_a_coef=[2.790e-03, 6.400e-03],
                                total_abs=6.585e-02,
                                total_backscat=2.594e-03),
                0.65: WaterData(wavelength=0.65,
                                refrac_real=1.338e+00,
                                refrac_imag=1.640e-08,
                                base_abs=3.500e-01,
                                base_backscat=9.379e-04,
                                whitecap_refl=3.544e-01,
                                chl_a_coef=[3.286e-03, 5.200e-03],
                                total_abs=3.518e-01,
                                total_backscat=1.879e-03),
                0.87: WaterData(wavelength=0.87,
                                refrac_real=1.334e+00,
                                refrac_imag=3.714e-07,
                                base_abs=5.365e+00,
                                base_backscat=2.662e-04,
                                whitecap_refl=2.488e-01,
                                chl_a_coef=[0.000e+00, 0.000e+00],
                                total_abs=5.365e+00,
                                total_backscat=1.239e-03),
                1.24: WaterData(wavelength=1.24,
                                refrac_real=1.327e+00,
                                refrac_imag=3.551e-05,
                                base_abs=3.599e+02,
                                base_backscat=5.759e-05,
                                whitecap_refl=7.120e-02,
                                chl_a_coef=[0.000e+00, 0.000e+00],
                                total_abs=3.599e+02,
                                total_backscat=8.667e-04),
                1.375: WaterData(wavelength=1.375,
                                 refrac_real=1.325e+00,
                                 refrac_imag=1.220e-04,
                                 base_abs=1.115e+03,
                                 base_backscat=3.685e-05,
                                 whitecap_refl=6.400e-03,
                                 chl_a_coef=[0.000e+00, 0.000e+00],
                                 total_abs=1.115e+03,
                                 total_backscat=7.944e-04),
                1.6: WaterData(wavelength=1.6,
                               refrac_real=1.323e+00,
                               refrac_imag=8.550e-05,
                               base_abs=6.715e+02,
                               base_backscat=1.915e-05,
                               whitecap_refl=0.000e+00,
                               chl_a_coef=[0.000e+00, 0.000e+00],
                               total_abs=6.715e+02,
                               total_backscat=7.056e-04),
                2.13: WaterData(wavelength=2.13,
                                refrac_real=1.313e+00,
                                refrac_imag=5.729e-04,
                                base_abs=3.380e+03,
                                base_backscat=5.563e-06,
                                whitecap_refl=0.000e+00,
                                chl_a_coef=[0.000e+00, 0.000e+00],
                                total_abs=3.380e+03,
                                total_backscat=5.771e-04),
                3.7: WaterData(wavelength=3.7,
                               refrac_real=1.374e+00,
                               refrac_imag=3.600e-03,
                               base_abs=1.223e+04,
                               base_backscat=5.120e-07,
                               whitecap_refl=0.000e+00,
                               chl_a_coef=[0.000e+00, 0.000e+00],
                               total_abs=1.223e+04,
                               total_backscat=4.188e-04)}
