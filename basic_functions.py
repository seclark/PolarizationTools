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

@np.vectorize
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
    
