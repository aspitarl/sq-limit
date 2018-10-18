# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

##Why are numbers not correct?

#working but off by a factor of 4??

#can't figure things out exactly with the factors but get the correct result setting the solid angle to pi. see the one note notebook.

# Just some constants for the upcoming math

import numpy as np
import scipy.constants as constants
import matplotlib.pyplot as plt
import matplotlib
import pandas as pd

font = {'size' : 14}

matplotlib.rc('font',**font)


c = constants.value('speed of light in vacuum')
h = constants.value('Planck constant in eV s')
e = constants.value('elementary charge')
k = constants.value('Boltzmann constant in eV/K')

pi = 3.1415

# Globals
Tcell = 300  # Kelvin
# Energy Gap
Egap = 1.1  #electron volts

r_earth = 6e6
r_sun = 6.95e8
d_sun = 1.50e11

def photons_above_bandgap(egap, spectrum):
    """
    Counts number of photons above given bandgap
    
    input: #/m^2*s*eV
    output: #/m^2*s
    """
    indexes = np.where(spectrum.index > egap)
    y = spectrum.iloc[indexes]
    x = y.index
    return np.trapz(y[::-1], x[::-1])
    
def recomb_rate(egap, spectrum, voltage):
    print( 'recomb rate')
    return e * rr_v(egap, spectrum, 0) * np.exp(voltage / (k * Tcell))

def rr_v(egap, spectrum, voltage):

    const = (2 * np.pi) / (c**2 * h**3)

    E = spectrum.index[::-1, ]  # in increasing order of bandgap energy
    egap_index = np.where(spectrum.index>= egap)
    numerator = E**2
    exponential_in = (E - voltage) / (k * Tcell)
    denominator = np.exp(exponential_in) - 1
    integrand = numerator / denominator

    integral = np.trapz(integrand[egap_index], E[egap_index])

    result = const * integral
    return result


#Was here before, errors when trying to remove
def rr0(egap, spectrum):

    const = (2 * np.pi) / (c**2 * h**3)

    E = spectrum.index[::-1, ]  # in increasing order of bandgap energy
    egap_index = np.where(spectrum.index>= egap)
    numerator = E**2
    exponential_in = E / (k * Tcell)
    denominator = np.exp(exponential_in) - 1
    integrand = numerator / denominator

    integral = np.trapz(integrand[egap_index], E[egap_index])

    result = const * integral
    return result

def recomb_rate_v(egap, spectrum, voltage):
    print( 'recomb rate with correct voltage ')
    return e * rr_v(egap, spectrum, voltage)

def current_density(egap, spectrum, voltage):
    PAG = photons_above_bandgap(egap, spectrum) 
    return e * (PAG- rr_v(egap, spectrum,0) * (np.exp(voltage / (k * Tcell)) -1) )
    
    
def jsc(egap, spectrum):
    return current_density(egap, spectrum, 0)

def voc(egap, spectrum):
    return (k * Tcell) * np.log((photons_above_bandgap(egap, spectrum) / rr_v(egap, spectrum,0)) +1)


def v_at_mpp(egap, spectrum):
    v_open = voc(egap, spectrum)
    v = np.linspace(0, v_open)
    index = np.where(v * current_density(egap, spectrum, v)==max(v * current_density(egap, spectrum, v)))
    return v[index][0]

def j_at_mpp(egap, spectrum):
    return max_power(egap, spectrum) / v_at_mpp(egap, spectrum)


def max_power(egap, spectrum):
    v_open = voc(egap, spectrum)
    v = np.linspace(0, v_open)
    #index = np.where(v * current_density(egap, spectrum, v)==max(v * current_density(egap, spectrum, v)))
    return max(v * current_density(egap, spectrum, v))

def int_irr(egap, spectrum):
    irradiance =  np.trapz(spectrum.iloc[::-1] * e * spectrum.index[::-1], spectrum.index[::-1])
    return irradiance

def max_eff(egap, spectrum):
    return max_power(egap, spectrum) / int_irr(egap, spectrum)

def max_eff_array(spectra_ph_all,index):
    
    max_eff_all_Si = pd.Series(index = index)
    max_eff_all_GaSb = pd.Series(index = index)
    max_pow_all_Si = pd.Series(index = index)
    max_pow_all_GaSb = pd.Series(index = index)

    i=0
    for idx in index:
        spectra_ph = spectra_ph_all[idx]

        max_eff_all_Si[idx] = max_eff(Egap, spectra_ph)*100
        max_eff_all_GaSb[idx] = max_eff(0.7, spectra_ph)*100
        
        max_pow_all_Si[idx] = max_power(Egap, spectra_ph)
        max_pow_all_GaSb[idx] = max_power(0.7, spectra_ph)
        i=i+1
    
    return max_eff_all_Si, max_eff_all_GaSb, max_pow_all_Si,max_pow_all_GaSb

