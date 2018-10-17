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

###Todo understand better irradiance vs radiance, may be causing fudge factor of 4
def rad_blackbody(E_ph,constants): #Todo: switch to kwargs
    """High level function to generate a blackbody spectral radiance depending on temperature and geometrical constants"""
    BB = spect_irrad(E_ph,constants['Temp'])
    BB = irrad_to_rad(BB,constants['solidangle'] ,constants['emitterarea'], constants['absorberarea'])
    return BB

def spect_irrad(E_ph,T):
    """
    Generate spectral irradance for temperature T
    
    # units W/m^2*eV*sr
    """
    a = (2*(E_ph**3))/(h**3*c**2)
    b = 1/( np.exp((E_ph)/(k*T)) -1 )
    intensity = a*b
    # units (eV/s)/m^2*eV*sr
    
    intensity = intensity*e 
     
    spectra = pd.Series(intensity, index = E_ph)
    return spectra

def irrad_to_rad(spectra, solidangle ,emitterarea, absorberarea):
    """
    Convert spectral irradiance to spectral radiance.
    
    Calculate total energy emitting from body of area emitterarea over angle solidangle, then cast that energy on an absorber of area absorberarea. 
    # units W/m^2*eV
     """
    spectra = spectra*emitterarea   
    # units W/eV*sr
    spectra = spectra*solidangle/4   ### fudge factor of 4 for now. see above/notebook.
    # units W/eV
    spectra = spectra/absorberarea  
    
    return spectra


def power_to_photons(spectrum):
    """
    convert to photons from energy
    
    # units #/eV    
    """
    spectrum = spectrum/spectrum.index * 1/e    
    
    return spectrum

def solid_angle_sun (r_e, d_sun):
    """Calculate the solid angle of earth with respect to sun"""
    return 4*(pi)*((pi*r_e**2)/(4*pi*d_sun**2))

def stephan(T):
    """
    Return the intensity from the stephan boltzman equation 
    
    #units W/m^2
    """
    result = (5.670367e-8)*(T**4)
    
    return result


###basic math Function generation, used to make emissivity and filter arrays

def stepfn(e_low,high,E_0,E_ph):
    temp = np.copy(E_ph)
    
    for i in range(len(E_ph)):
        E = E_ph[i]
        if (E<E_0):
            temp[i] = e_low
        else: 
            temp[i] = high
    
    emissivity = pd.Series(temp,index= E_ph)
    return emissivity

def lorentzian(low,high,E_0,width,E_ph):
    temp = np.zeros(len(E_ph))
    for i in range(len(E_ph)):
        E = E_ph[i]
        temp[i] = low + (high-low)/(1 + ((E-E_0)**2/(width**2))  )   
    emissivity = pd.Series(temp,index= E_ph)
    return emissivity



if __name__ == '__main__':
    E_ph = np.arange(0.01, 10,0.001) 
    E_ph = np.flip(E_ph,0)

    E_gaps = np.arange(0.3, 5,0.01) 
    E_gaps = np.flip(E_gaps,0)

    constants = {}

    constants['Temp'] = 1800
    constants['solidangle'] = 2*pi
    constants['emitterarea'] = 1
    constants['absorberarea'] = 1
    
    constants['emissivity'] = stepfn(1,1,1.1,E_ph)