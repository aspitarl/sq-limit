# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import pandas as pd
import numpy as np

#Load filter data

def read_filter(filepath,refspect):
    """loads in filter data and interpolates data points to match reference speectrum"""
    filt = pd.read_csv(filepath, usecols = [0,1], index_col = 0)
    return reshape_spect(filt['% Transmission'],refspect)

#Load QE data

def read_QE(filepath,refspect):
    QEcsv = pd.read_csv(filepath, names = ['QE', 'Wavelength'], index_col=1)
    return reshape_spect(QEcsv['QE'],refspect)


def reshape_spect(spectrum, refspect):
    """converts a spectral array to be in eV and interpolates the data to match refspect"""
    #flip and switch to eV
    wavelength = spectrum.index
    ev = 1240/wavelength
    if (ev[0] > ev[-1]):
        ev = np.flip(ev,0)
        spectrum = np.flip(spectrum, 0)
    #interpolate
    spect_interp = np.interp(refspect.index,ev,spectrum)
    ev_interp = np.interp(refspect.index,ev,ev)
    spect = pd.Series(spect_interp, index = ev_interp)
    return spect

