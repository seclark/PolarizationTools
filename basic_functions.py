from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import pyfits
import pylab
import time
import copy
import math
import os.path
import cPickle as pickle
from matplotlib import rc
from astropy.io import fits

vectorize = np.vectorize

@vectorize
def angle_residual(ang1, ang2, degrees = True):
    """
    Determines circular difference between two angles.
    Assumes radian inputs unless degrees = True
    """

    if degrees is True:
        ang1 = np.radians(ang1)
        ang2 = np.radians(ang2)

    dang_num = (np.sin(2*ang1)*np.cos(2*ang2) - np.cos(2*ang1)*np.sin(2*ang2))
    dang_denom = (np.cos(2*ang1)*np.cos(2*ang2) + np.sin(2*ang1)*np.sin(2*ang2))
    dang = 0.5*np.arctan2(dang_num, dang_denom)
    
    if degrees is True:
        dang = np.degrees(dang)
    
    return dang

@vectorize
def sigma_psi_P(Q, U, sig_QQ, sig_UU):
    """
    Output:: sigma_psi: error on polarization angle [degrees]
             sigma_P: error on polarized intensity
    """
    
    Psquared = Q**2 + U**2
    P = np.sqrt(Psquared)
    sig_P = np.sqrt((1/Psquared)*(Q**2*sig_QQ**2 + U**2*sig_UU**2))
    
    sig_psi = 28.65*np.sqrt((Q**2*sig_UU**2 + U**2*sig_QQ**2)/(Q**2*sig_QQ**2 + U**2*sig_UU**2))*(sig_P/P)
    
    return sig_psi, sig_P
    
    
    
    