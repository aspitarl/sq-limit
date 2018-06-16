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
import pandas as pd

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



def gen_emissivity(e_low,e_high,E_cutoff,E_bb):
    emissivity = np.copy(E_bb)
    
    i=0
    for E in E_bb:
        if (E<E_cutoff):
            emissivity[i] = e_low
        else: 
            emissivity[i] = e_high
        i=i+1    
    
    return emissivity

def lor_emissivity(e_bg,e_high,E_cutoff,w,E_bb):
    emissivity = np.copy(E_bb)
    
    i=0
    for E in E_bb:
        emissivity[i] = e_bg + (e_high-e_bg)/(1 + ((E-E_cutoff)**2/(w**2))  )
        i=i+1   
    return emissivity


def stephan(T):
    result = (5.670367e-8)*(T**4)
    #units W/m^2
    return result

def solid_angle_sun (r_e, d_sun):
    return 4*(pi)*((pi*r_e**2)/(4*pi*d_sun**2))

def spect_rad(T,E_bb, emissivity, powfactor = 1):
    a = (2*(E_bb**3))/(h**3*c**2)
    b = 1/( np.exp((E_bb)/(k*T)) -1 )
    intensity = a*b*powfactor
    # units (eV/s)/m^2*eV*sr
    
    intensity = intensity*e 
     # units W/m^2*eV*sr
     
    intensity = intensity*emissivity
    spectra = np.transpose(np.stack((E_bb,intensity)))
    return spectra

def rad_to_rad(spectrum, solidangle ,emitterarea, absorberarea):
    spectra =  spectrum[:,1]
    spectra = spectra*emitterarea   
    # units W/eV*sr
    spectra = spectra*solidangle/4   ### fudge factor of 4 for now. see above.
    # units W/eV
    spectra = spectra/absorberarea  
    # units W/m^2*eV
    spectrum[:,1] = spectra
    
    return spectrum

def solar(T,E_bb,emissivity, powfactor = 1):
    BB = spect_rad(T,E_bb,emissivity, powfactor)
    BB = rad_to_rad(BB, solid_angle_sun(r_earth,d_sun), 4*pi*r_sun**2, (pi*r_earth**2))
    # units W/eV
    return BB

####Analysis


    
# convert to photons from energy
def power_to_photons(spectrum):
    converted = np.copy(spectrum)
    converted[:, 1] = converted[:, 1]/converted[:,0] * 1/e    
    # units #/eV    
    return converted

def photons_above_bandgap(egap, spectrum):
    """Counts number of photons above given bandgap"""
    indexes = np.where(spectrum[:, 0] > egap)
    y = spectrum[indexes, 1][0]
    x = spectrum[indexes, 0][0]
    return np.trapz(y[::-1], x[::-1])
    
def rr0(egap, spectrum):

    const = (2 * np.pi) / (c**2 * h**3)

    E = spectrum[::-1, ]  # in increasing order of bandgap energy
    egap_index = np.where(E[:, 0] >= egap)
    numerator = E[:, 0]**2
    exponential_in = E[:, 0] / (k * Tcell)
    denominator = np.exp(exponential_in) - 1
    integrand = numerator / denominator

    integral = np.trapz(integrand[egap_index], E[egap_index, 0])

    result = const * integral
    return result[0]

def recomb_rate(egap, spectrum, voltage):
    print( 'recomb rate')
    return e * rr0(egap, spectrum) * np.exp(voltage / (k * Tcell))

def rr_v(egap, spectrum, voltage):

    const = (2 * np.pi) / (c**2 * h**3)

    E = spectrum[::-1, ]  # in increasing order of bandgap energy
    egap_index = np.where(E[:, 0] >= egap)
    numerator = E[:, 0]**2
    exponential_in = (E[:, 0] - voltage) / (k * Tcell)
    denominator = np.exp(exponential_in) - 1
    integrand = numerator / denominator

    integral = np.trapz(integrand[egap_index], E[egap_index, 0])

    result = const * integral
    return result[0]

def recomb_rate_v(egap, spectrum, voltage):
    print( 'recomb rate with correct voltage ')
    return e * rr_v(egap, spectrum, voltage)

def current_density(egap, spectrum, voltage):
    return e * (photons_above_bandgap(egap, spectrum) - rr0(egap, spectrum) * (np.exp(voltage / (k * Tcell))) - 1)
    #return e * (photons_above_bandgap(egap, spectrum) - rr_v(egap, spectrum,voltage))
def jsc(egap, spectrum):
    return current_density(egap, spectrum, 0)

def voc(egap, spectrum):
    return (k * Tcell) * np.log(photons_above_bandgap(egap, spectrum) / rr0(egap, spectrum))

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
    irradiance =  np.trapz(spectrum[::-1, 1] * e * spectrum[::-1, 0], spectrum[::-1, 0])
    return irradiance

def max_eff(egap, spectrum):
    irradiance =  np.trapz(spectrum[::-1, 1] * e * spectrum[::-1, 0], spectrum[::-1, 0])
    return max_power(egap, spectrum) / irradiance

def max_eff_temp(E_bb,emissivity, Tmin,Tmax,dT,sourcetype = 0):
    # sourcetype 1 for sun and 0 for full angle
    Temps = pd.Series(np.arange(Tmin,Tmax,dT))
    spectra_ph_all = pd.Series(index = Temps,dtype=object)
    max_eff_all_Si = pd.Series(index = Temps)
    max_eff_all_GaSb = pd.Series(index = Temps)
    max_pow_all_Si = pd.Series(index = Temps)
    max_pow_all_GaSb = pd.Series(index = Temps)

    i=0
    for temp in Temps:
        if sourcetype:
            spectra = solar(temp,E_bb,emissivity)
        else:
            spectra = spect_rad(temp,E_bb,emissivity)
            spectra = rad_to_rad(spectra, solidangle = 4*pi, emitterarea = 1, absorberarea = 1)
        spectra_ph = power_to_photons(spectra)

        spectra_ph_all[temp] = spectra_ph
        max_eff_all_Si[temp] = max_eff(Egap, spectra_ph)*100
        max_eff_all_GaSb[temp] = max_eff(0.7, spectra_ph)*100

        max_pow_all_Si[temp] = max_power(Egap, spectra_ph)
        max_pow_all_GaSb[temp] = max_power(0.7, spectra_ph)
        i=i+1
    
    return max_eff_all_Si, max_eff_all_GaSb, max_pow_all_Si,max_pow_all_GaSb


def max_eff_power(E_bb,emissivity):
    
    powers = pd.Series(np.arange(10, 1000,10))
    spectra_ph_all = pd.Series(index = powers,dtype=object)
    max_eff_all_Si = pd.Series(index = powers)
    max_eff_all_GaSb = pd.Series(index = powers)
    max_pow_all_Si = pd.Series(index = powers)
    max_pow_all_GaSb = pd.Series(index = powers)

    i=0
    for power in powers:
        spectra = spect_rad(5500,E_bb,emissivity, powfactor = power )
        spectra = rad_to_rad(spectra, solidangle = 4*pi, emitterarea = 1, absorberarea = 1)
        spectra_ph = power_to_photons(spectra)
        
        spectra_ph_all[power] = spectra_ph
        max_eff_all_Si[power] = max_eff(Egap, spectra_ph)*100
        max_eff_all_GaSb[power] = max_eff(0.7, spectra_ph)*100
        
        max_pow_all_Si[power] = max_power(Egap, spectra_ph)
        max_pow_all_GaSb[power] = max_power(0.7, spectra_ph)
        i=i+1
    
    return max_eff_all_Si, max_eff_all_GaSb, max_pow_all_Si,max_pow_all_GaSb


###Plots
    
def photons_above_bandgap_plot(spectrum,E_gaps):
    """Plot of photons above bandgap as a function of bandgap"""
    PAG = np.copy(E_gaps)
    i=0
    for Eg in E_gaps:
        PAG[i] = photons_above_bandgap(Eg,spectrum)
        i=i+1   
    plt.plot(E_gaps,PAG)

    p_above_1_1 = photons_above_bandgap(Egap, spectrum)
    plt.plot([Egap], [p_above_1_1], 'ro')
    plt.text(Egap+0.05, p_above_1_1, '{}eV, {:.4}'.format(Egap, p_above_1_1))

    plt.xlabel('$E_{gap}$ (eV)')
    plt.ylabel('# Photons $m^{-2}s^{-1}$')
    plt.title('Number of above-bandgap \nphotons as a function of bandgap')
    plt.show()

def ideal_jsc_plot(spectrum,E_gaps):
    """Plot of photons above bandgap as a function of bandgap"""
#    a = np.copy(spectrum)
    Jscs = np.copy(E_gaps)
    i=0
    for Eg in E_gaps:
        Jscs[i] = jsc(Eg,spectrum)
        i=i+1   
    
    plt.plot(E_gaps, Jscs)
    e_gap = 1.1
    p_above_1_1 = jsc(e_gap, spectrum)
    plt.plot([e_gap], [p_above_1_1], 'ro')
    plt.text(e_gap+0.05, p_above_1_1, '{}eV, {:.4}'.format(e_gap, p_above_1_1))

    plt.xlabel('$E_{gap}$ (eV)')
    plt.ylabel('$J_{SC}$ $Am^{-2}$')
    plt.title('Ideal short-circuit current')


def ideal_voc_plot(spectrum, E_gaps):
    """Plot of the ideal open circuit voltage as a function of bandgap"""
    Vocs = np.copy(E_gaps)
    i=0
    for Eg in E_gaps:
        Vocs[i] = voc(Eg,spectrum)
        i=i+1
    
    plt.plot(E_gaps, Vocs)
    plt.plot(E_gaps, E_gaps)
    e_gap = 1.1
    p_above_1_1 = voc(e_gap, spectrum)
    plt.plot([e_gap], [p_above_1_1], 'ro')
    plt.text(e_gap+0.05, p_above_1_1, '{}eV, {:.4}'.format(e_gap, p_above_1_1))

    plt.xlabel('$E_{gap}$ (eV)')
    plt.ylabel('$V_{OC}$ (V)')
    plt.xlim((0,3.5))
    plt.ylim((0,3.5))
    plt.title('Ideal open-circuit voltage. Straight line is bandgap.')
    
    
        
def iv_curve_plot(egap, spectrum, power=False):
    """Plots the ideal IV curve, and the ideal power for a given material"""
    v_open = voc(egap, spectrum)
    v = np.linspace(0, v_open)

    fig, ax1 = plt.subplots()
    p =  v * current_density(egap, spectrum, v)
    i =  current_density(egap, spectrum, v)
    
    ax1.plot(v, i)
    ax1.set_xlabel('Voltage (V)')
    ax1.set_ylabel('Current density $J$ ($Am^{-2}$)')
    ax1.legend(['Current'], loc=3)
    
    if(power):
        ax2 = ax1.twinx()
        ax2.plot(v, p, color='orange')
        ax2.set_ylabel('Power generated ($W$)')
        ax2.legend(['Power'], loc=3)
    return

def sq_limit_plot(spectrum,E_gaps):
    # Plot the famous SQ limit
    SQlim = np.copy(E_gaps)
    i=0
    for Eg in E_gaps:
        SQlim[i] = max_eff(Eg,spectrum)
        i=i+1   
    plt.plot(E_gaps,SQlim)
    # Not plotting whole array becase some bad values happen
    #plt.plot(a[2:, 0], a[2:, 1])
    e_gap = Egap
    p_above_1_1 = max_eff(e_gap, spectrum)
    plt.plot([e_gap], [p_above_1_1], 'ro')
    plt.text(e_gap+0.05, p_above_1_1, '{}eV, {:.4}'.format(e_gap, p_above_1_1))

    plt.xlabel('$E_{gap}$ (eV)')
    plt.ylabel('Max efficiency')
    plt.title('SQ Limit')