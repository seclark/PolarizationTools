from __future__ import division
import numpy as np
import cPickle
import matplotlib
matplotlib.rcParams['backend'] = "Qt4Agg"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
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
def sigma_psi_P(Q, U, sig_QQ, sig_UU, degrees = True):
    """
    Output:: sigma_psi: uncertainty on polarization angle [degrees]
             sigma_P: uncertainty on polarized intensity
    """
    
    Psquared = Q**2 + U**2
    P = np.sqrt(Psquared)
    sig_P = np.sqrt((1/Psquared)*(Q**2*sig_QQ**2 + U**2*sig_UU**2))
    
    sig_psi = 28.65*np.sqrt((Q**2*sig_UU**2 + U**2*sig_QQ**2)/(Q**2*sig_QQ**2 + U**2*sig_UU**2))*(sig_P/P)
    
    if degrees is False:
        sig_psi = np.radians(sig_psi)
    
    return sig_psi, sig_P

@vectorize
def polarization_angle(Q, U, negU = False):
    """
    Returns polarization angle in radians.
    If negU is True, multiplies U by -1. Useful for converting from IAU <-> Planck standard.
    NB: returns polarization angle, *not* B-field angle.
    """
    
    if negU is True:
        pol_ang = np.mod(0.5*np.arctan2(-U, Q), np.pi)
    else:
        pol_ang = np.mod(0.5*np.arctan2(U, Q), np.pi)
    
    return pol_ang
    
@vectorize
def mod_halfpolar_center_0(angle):
    """
    Returns angle (in radians) on a domain from -pi/2 to +pi/2. 
    """
    
    angle = np.mod(angle, np.pi)
    
    if angle > np.pi/2:
        angle = angle - np.pi
        
    return angle
    