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



def gen_emissivity(e_low,e_high,E_cutoff,E_ph):
    emissivity = np.copy(E_ph)
    
    i=0
    for E in E_ph:
        if (E<E_cutoff):
            emissivity[i] = e_low
        else: 
            emissivity[i] = e_high
        i=i+1    
    
    return emissivity

def lor_emissivity(e_bg,e_high,E_cutoff,w,E_ph):
    emissivity = np.copy(E_ph)
    
    i=0
    for E in E_ph:
        emissivity[i] = e_bg + (e_high-e_bg)/(1 + ((E-E_cutoff)**2/(w**2))  )
        i=i+1   
    return emissivity


def stephan(T):
    result = (5.670367e-8)*(T**4)
    #units W/m^2
    return result

def solid_angle_sun (r_e, d_sun):
    return 4*(pi)*((pi*r_e**2)/(4*pi*d_sun**2))

def spect_rad(E_ph,T, emissivity, powfactor = 1):
    a = (2*(E_ph**3))/(h**3*c**2)
    b = 1/( np.exp((E_ph)/(k*T)) -1 )
    intensity = a*b
    # units (eV/s)/m^2*eV*sr
    
    intensity = intensity*e 
     # units W/m^2*eV*sr
     
    intensity = intensity*emissivity
    spectra = np.transpose(np.stack((E_ph,intensity)))
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


def gen_spectrum(E_ph,constants):
    BB = spect_rad(E_ph,constants['Temp'],constants['emissivity'])
    BB = rad_to_rad(BB,constants['solidangle'] ,constants['emitterarea'], constants['absorberarea'])
    return BB

def gen_spectrum_lor(E_ph,I_bg,I_high,E_cutoff,w):
    intensity = np.copy(E_ph)  
    i=0
    for E in E_ph:
        intensity[i] = I_bg + (I_high-I_bg)/(1 + ((E-E_cutoff)**2/(w**2))  )
        i=i+1  
    
    spectra = np.transpose(np.stack((E_ph,intensity)))
    return spectra

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
    PAG = photons_above_bandgap(egap, spectrum) 
    return e * (PAG- rr0(egap, spectrum) * (np.exp(voltage / (k * Tcell)) - 1) )
#    j = np.copy(voltage)
#    for i in range(len(voltage)):
#        j[i] = e * (photons_above_bandgap(egap, spectrum) - rr_v(egap, spectrum,voltage[i]))
#    return j
    
    
def jsc(egap, spectrum):
    return current_density(egap, spectrum, 0)

def voc(egap, spectrum):
    return (k * Tcell) * np.log((photons_above_bandgap(egap, spectrum) / rr0(egap, spectrum)) +1)


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


###Plots
    
def em_ir_ph_plot(E_ph, constants, BB, BB_ph):
    w,h =plt.figaspect(1.5) 
    fig, ax = plt.subplots(3,1,figsize = (w,h) )
    ax[0].plot(E_ph,constants['emissivity'])
    ax[0].set_xlim(0,4)
    ax[0].set_ylabel('Emissivity')
    
    ax[1].plot(BB[:,0],BB[:,1] )
    ax[1].set_xlim(0,4)
    ax[1].set_ylabel('Spectral Irradiance \n ($Wm^{-2}eV^{-1}$)')
    
    ax[2].plot(BB_ph[:,0],BB_ph[:,1] )
    ax[2].set_xlim(0,4)
    ax[2].set_ylabel('Photon Irradiance\n (# $m^{-2}eV^{-1}$)')
    plt.xlabel('Photon Energy (eV)')
    
    plt.tight_layout()
    
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
#    v_open = 10
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