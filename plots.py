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
import detailedbalance as db

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



###Plots
    
def em_ir_ph_plot(E_ph, constants, BB, BB_ph):
    matplotlib.rc('font',**font)
    
    w,h =plt.figaspect(1.5) 
    fig, ax = plt.subplots(3,1,figsize = (w,h) )
    ax[0].plot(E_ph,constants['emissivity'])
    ax[0].set_xlim(0,3)
    ax[0].set_ylabel('Emissivity')
    
    ax[1].plot(BB[:,0],BB[:,1] )
    ax[1].set_xlim(0,3)
    ax[1].set_ylabel('Spectral Irradiance \n ($Wm^{-2}eV^{-1}$)')
    
    ax[2].plot(BB_ph[:,0],BB_ph[:,1] )
    ax[2].set_xlim(0,3)
    ax[2].set_ylabel('Photon Irradiance\n (# $m^{-2}eV^{-1}$)')
    plt.xlabel('Photon Energy (eV)')
    
    plt.tight_layout()
    
def photons_above_bandgap_plot(spectrum,E_gaps):
    """Plot of photons above bandgap as a function of bandgap"""
    PAG = np.copy(E_gaps)
    i=0
    for Eg in E_gaps:
        PAG[i] = db.photons_above_bandgap(Eg,spectrum)
        i=i+1   
    plt.plot(E_gaps,PAG)

    p_above_1_1 = db.photons_above_bandgap(Egap, spectrum)
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
        Jscs[i] = db.jsc(Eg,spectrum)
        i=i+1   
    
    plt.plot(E_gaps, Jscs)
    e_gap = 1.1
    p_above_1_1 = db.jsc(e_gap, spectrum)
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
        Vocs[i] = db.voc(Eg,spectrum)
        i=i+1
    
    plt.plot(E_gaps, Vocs)
    plt.plot(E_gaps, E_gaps)
    e_gap = 1.1
    p_above_1_1 = db.voc(e_gap, spectrum)
    plt.plot([e_gap], [p_above_1_1], 'ro')
    plt.text(e_gap+0.05, p_above_1_1, '{}eV, {:.4}'.format(e_gap, p_above_1_1))

    plt.xlabel('$E_{gap}$ (eV)')
    plt.ylabel('$V_{OC}$ (V)')
    plt.xlim((0,3.5))
    plt.ylim((0,3.5))
    plt.title('Ideal open-circuit voltage. Straight line is bandgap.')
    
    
        
def iv_curve_plot(egap, spectrum, power=False):
    """Plots the ideal IV curve, and the ideal power for a given material"""
    v_open = db.voc(egap, spectrum)
#    v_open = 10
    v = np.linspace(0, v_open)

    fig, ax1 = plt.subplots()
    p =  v * db.current_density(egap, spectrum, v)
    i =  db.current_density(egap, spectrum, v)
    
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
        SQlim[i] = db.max_eff(Eg,spectrum)
        i=i+1   
    plt.plot(E_gaps,SQlim)
    # Not plotting whole array becase some bad values happen
    #plt.plot(a[2:, 0], a[2:, 1])
    e_gap = Egap
    p_above_1_1 = db.max_eff(e_gap, spectrum)
    plt.plot([e_gap], [p_above_1_1], 'ro')
    plt.text(e_gap+0.05, p_above_1_1, '{}eV, {:.4}'.format(e_gap, p_above_1_1))

    plt.xlabel('$E_{gap}$ (eV)')
    plt.ylabel('Max efficiency')
    plt.title('SQ Limit')