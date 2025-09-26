# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 12:25:31 2024

@author: as836
"""

import os
import numpy as np
from dqc_kernel import *
from scipy.interpolate import interp1d 
import matplotlib.pyplot as plt
import seaborn as sns
import time

def scmm(data):    
    return (data - data.min()) / (data.max() - data.min())

def movav(data, n=3):
    _ndat = np.pad(data, (int((n-1)/2), int((n-1)/2)), 'constant', constant_values = (data[0], data[-1]))
    return np.convolve(_ndat, np.ones(n)/n, mode='valid')


# define some constants
_uB = 9.274*1E-24
_h = 6.626*1E-34
_g = 2.0023

# molecular parameters

## set P(r) 
_rrange = np.arange(1.5, 5.0, 0.05)
_mu1s, _sm1s = [2.5, 3.0, 4.0], [0.2, 0.36, 0.3] # in nm
_pr1 = (np.exp(-(_rrange - _mu1s[0])**2/(2*_sm1s[0]**2)) + 
        np.exp(-(_rrange - _mu1s[1])**2/(2*_sm1s[1]**2)) +
        np.exp(-(_rrange - _mu1s[2])**2/(2*_sm1s[2]**2)))
_pr1 /= _pr1.sum()
_prf = interp1d(_rrange, _pr1)

# set pulse amplitudes and lengths 
_wps = [125., 62.5, 31.25]
_tps = [4, 8, 16]

_tau1s = np.arange(10, 210, 5) # list of t1 values in ns

# run simulations
for k in range(len(_wps)):
    # other parameters
    _thtab = np.linspace(0, np.pi/2, 100)
    _pth = np.sin(_thtab)
    _pth /= _pth.sum()
    
    # pulse parameters
    _wid = 30 - k*10  
    _mult = 2
    _omtab = np.linspace(-_mult*_wid, _mult*_wid, 100)
    _pom = np.exp(-_omtab**2/(2*_wid**2))
    _pom /= _pom.sum()
        
    _wp = _wps[k] # pulse field in MHz
    _tp = _tps[k]*1E-3 # pulse length in us
    
    # pulse separations in time-domain
    _tm = 2.5
    _dt = 0.016
    _ttab = np.arange(0, _tm, _dt)
    
    # simulate the signal
    _mmax = 1*1E4 ## max number of iterations
    


    for j in range(len(_tau1s)):
        _resY = 0

        _t1 = _tau1s[j]*1E-3

        _m = 1
        while _m <= _mmax:
    
            _om1, _om2 = np.random.choice(_omtab, p = _pom), np.random.choice(_omtab, p = _pom)
        
            if _om1==_om2:
                _om2 += 1E-4
    
            _r = np.random.choice(_rrange, p = _pr1)
            _th = np.random.choice(_thtab, p = _pth)
        
            _dip = 2*np.pi*52.04/_r**3*(1-3*np.cos(_th)**2)
    
            _sig1 = sig(_om1, _om2, _dip, (_tm - _ttab)/2, _t1, (_tm + _ttab)/2, _wp, _tp, path = "1")
        
            _resY += _sig1.real
        
            _m += 1
        
        _resY /= _mmax
    

        # save the simulation
        _wpath = r'C:\Users\as836\Desktop\Tau_DQC_Analysis\Simulations'
    
        _fname = 'DQC_20A_pi-'+str(int(_tps[k]))+'ns_t1-'+ str(int(_tau1s[j]))+'ns.txt'
    
        np.savetxt(os.path.join(_wpath,_fname), np.c_[_ttab, _resY])