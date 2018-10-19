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
    filt_data = filt['% Transmission']
    ev = 1240/filt.index
    filt_interp = np.interp(refspect.index,ev,filt_data)
    ev_interp = np.interp(refspect.index,ev,ev)
    filt = pd.Series(filt_interp, index = ev_interp)
    return filt

#Load QE data

def read_QE(filepath,refspect):
    QEcsv = pd.read_csv(filepath, names = ['QE', 'Wavelength'])
    ev = 1240/QEcsv['Wavelength']
    ev = np.flip(ev,0)
    QE = QEcsv['QE']
    QE = np.flip(QE,0)
    QE_interp = np.interp(refspect.index,ev,QE)
    ev_interp = np.interp(refspect.index,ev,ev)
    QE = pd.Series(QE_interp, index = ev_interp)
    return QE