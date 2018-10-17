# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 21:15:19 2018

@author: aspit
"""

import detailedbalance as db
import numpy as np
import scipy.constants as constants
import matplotlib.pyplot as plt
import pandas as pd
import importlib

#importlib.reload(db)

E_ph = np.arange(0.01, 10,0.001) 
E_ph = np.flip(E_ph,0)

E_gaps = np.arange(0.3, 5,0.01) 
E_gaps = np.flip(E_gaps,0)

sourcetype = 1 

constants = {}

# sourcetype 1 for sun and 0 for full angle. This should perhaps be combined with the max_eff_temp function
if(sourcetype):
    #Sun
    constants['Temp'] = 5750
    constants['solidangle'] = db.solid_angle_sun(r_earth,d_sun)
    constants['emitterarea'] = 4*pi*r_sun**2
    constants['absorberarea'] = pi*r_earth**2
else:
    #Custom
    constants['Temp'] = 5500
    constants['solidangle'] = 2*pi
    constants['emitterarea'] = 1
    constants['absorberarea'] = 1
    
constants['emissivity'] = db.gen_emissivity(1,1,1.1,E_ph)
#constants['emissivity'] = db.lor_emissivity(0,1,3,0.00001,E_ph)
   
BB = db.gen_spectrum(E_ph,constants)

#BB = db.gen_spectrum_lor(E_ph,0,1,1.12,0.00001)


#check integrated irradiance is ~1kW/m^2
integrate = -np.trapz(BB[:,1],BB[:,0] )

#BB[:,1] = BB[:,1]/integrate
BB[:,1] = BB[:,1]*10**-17
#BB[:,1] = 0

BB_ph = db.power_to_photons(BB)

#print(db.stephan(5750)*4*pi*r_sun**2)

print(voc(1.1,BB_ph))
db.iv_curve_plot(1.1, BB_ph)
print(rr_v(1.1,BB_ph,1.2))