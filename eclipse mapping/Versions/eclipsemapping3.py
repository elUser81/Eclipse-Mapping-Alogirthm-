from astropy.table import Table
import numpy as np
from scipy.optimize import fsolve
import scipy.stats as st
from tqdm import tqdm
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import astropy.units as u
import os
cwd = os.getcwd()

def get_starobrits_shape(D, q, inclination):
    '''
    :param D: distance between the stars
    :param q: mass ratio (mass of second star/ mass of white dwarf...M2/M1)
    :return: a function giving the xyz coordinates for both stars
    '''

    def orbits_pos(phase):
        '''
        :param phase: the phasepoint of the orbit
        :return: the xyz coordinates for both stars
        '''
        x = D * np.sin(phase)
        y = D * np.cos(inclination) * np.cos(phase)
        z = D * np.sin(inclination) * np.cos(phase)

        x1 = -x/(1 + (1/q))
        y1 = -y / (1 + (1 / q))
        z1 = -z / (1 + (1 / q))

        x2 = x / (1 + q)
        y2 = y / (1 + q)
        z2 = z / (1 + q)

        wdcoord = [x1,y1,z1]
        rgcoord = [x2,y2,z2]

        return wdcoord, rgcoord
    return orbits_pos






