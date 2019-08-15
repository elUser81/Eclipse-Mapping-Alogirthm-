"""
An eclipsemapping and imaging algorithm for Catalclysmic Variable star systems.

Produces an image of the accretion disk around the primary star
given the mass ratio, eccentricity, and light curve of the system.

Code written by Nathan Smith, Thesis research at the University of Dallas
nsmith1@udallas.edu


Assumptions about the system:
    -- The radius of the two stars is equal to the radius of their roche lobe potiental
    -- The eccentricity is approximately zero, with the primary star at the center of the orbit
    -- the accretion disk is within the orbital plane
    -- there's one more... but its in the paper but I forget off hand

"""

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

class Star_System:
    '''
    A star system object containing all of the orbital pramaeters and data assigned by the Imager object,
    such as the grid and orbit-space conversions unique to this star system
    '''


    def __init__(self,mass_ratio, light_curve, distance, semi_major_axis, orbital_inclination,
                 eccentricity, period = None):

        #mass ratio q = M2/M1 where M1 is the primary white dwarf star
        self.mass_ratio = mass_ratio
        self.period = period
        self.light_curve = light_curve
        self.distance = distance
        self.a = semi_major_axis
        self.e = eccentricity
        self.i = orbital_inclination
        a, e = self.a, self.e
        self.b = np.sqrt((a**2)*(1-(e*e)))
        self.c = self.a * self.e




        """Radius of the secondary star.. found an approximation of RL1 
        (the distance between the secondary star and the first lagrange point) on wikipedia.... it assumes M2 << M1. 
        Not sure if this is a valid assumtion to make but I'm going to try it out and can change it later if I need to.. 
        other option is to use Roche lobe potiental at RL1 equation and solve numerically """

        '''distance from smallest mass star too 1st lagrange point'''
        self.RL1 = self.a*(self.mass_ratio * (1/3))**(1/3)
        self.Radius_2nd = self.RL1


        '''distance from 1st Lagarange point to primary white dwarf'''
        self.RL1_from_primary = self.a - self.RL1

        """Radius of the 1st star, probably won't need it but am making space just in case """
        #self.Radius_1st = None

        "Solid angle of each pixel"
        self.theta = (self.RL1/self.distance)*np.sqrt(np.cos(self.i))

        '''Making space for the grid values.. the grid is centered around the white dwarf 
        star and will have N^2 pixels'''
        self.pixel_grid = np.ndarray(shape = (2,2), dtype= np.int32)
        '''length of each pixel, space for later '''
        self.pixel_length = 1

        '''contains the 3d numpy grid of the starsystem's fractional visibility matrix produced by the imager object
            making space for later'''
        self.vis_matrix = np.ndarray(shape = (2,2,2))

    '''the shape of orbit in terms of x, y, z coordinates at phase phi rotated about the y axis at angle i...
        these will be functions that are called through the object, normalized by RL1'''

    def orb_x(self, phi):
        #dividing by RL1  puts the units interms of RL1

        return - (self.a * np.cos(phi) * np.sin(self.i))/self.RL1

    def orb_y(self, phi):
        return (self.b * np.sin(phi))/self.RL1

    def orb_z(self, phi):
        return (self.a * np.cos(phi) * np.cos(self.i))/self.RL1

    '''the following methods return the coordinates of the white dwarf in rotated orbit space, normalized by RL1'''
    def wd_x(self):
        return self.c * np.sin(self.i)/self.RL1

    def wd_y(self):
        return 0

    def wd_z(self):
        return - self.c * np.cos(self.i)/self.RL1


    ''' The following methods return the pixel_grid coordinates transformed into the space of the rotated binary orbit.
    so... pixel grid point [0,0] will be the spatial coordinates of the white dwarf star since the pixel grid is 
    centered around it'''


    def pixel_orb_x(self, xp, yp):

        '''
        :param xp: the x coordinate of the pixel grid
        :param yp: the y coordinate of the pixel grid, keeping it here for consitency
        :return: the orbit space coordinate for xp
        '''

        #factor which centers the coordinates of the pixel grid
        centered = ((self.pixel_grid.shape[0] - 1)/2)

        #centering xp:
        xg = (xp - centered)

        #a conversion factor ... is this right? units dont make sense to me .... im going with it..
        beta = self.RL1/self.pixel_grid.shape[0]

        #foci coordinates
        xf = -self.c / self.RL1

        xs = -(beta*xg + xf) * np.sin(self.i) #nevermind..units are right.. beta is multiplied by something in
                                              #pixel space, canceling the pixel units

        return xs


    def pixel_orb_y(self, xp, yp):
        '''

        :param xp: x coordinate form pixel grid
        :param yp: y coordinate from pixel grid
        :return: y coordinate of pixel grid in orbit space
        '''

        # factor which centers the coordinates of the pixel grid
        centered = ((self.pixel_grid.shape[0] - 1) / 2)

        # centering xp:
        yg = (centered - yp)

        # a conversion factor ... is this right? units dont make sense to me .... im going with it..
        beta = self.RL1 / self.pixel_grid.shape[0]

        # foci coordinates
        yf = 0 / self.RL1

        ys = beta * yg + yf

        return ys


    def pixel_orb_z(self, xp, yp):

        '''

        :param xp: x coordinate form pixel grid
        :param yp: y coordinate from pixel grid
        :return: z coordinate of pixel grid in orbit space
         '''

        # factor which centers the coordinates of the pixel grid
        centered = ((self.pixel_grid.shape[0] - 1) / 2)

        # centering xp:
        xg = (xp - centered)

        # a conversion factor ... is this right? units dont make sense to me .... im going with it..
        beta = self.RL1 / self.pixel_grid.shape[0]

        # foci coordinates
        xf = -self.c / self.RL1

        zs = (beta * xg + xf) * np.cos(self.i)
        return zs


    def check_period(self):
        #will calculate the system's period from the light curve if period == None, ill write it later
        pass




class Imager:

    def __init__(self, N):
        #dimensions used in the pixel grid
        self.N = N
        self.pixel_grid = np.zeros(shape= (N,N), dtype= int)

    def assign_grid(self, stars):
        '''
        :param stars: a star system object
        :return: assigns a grid to stars.grid with given parameters
        '''

        stars.pixel_grid = self.pixel_grid
        stars.pixel_length = stars.RL1_from_primary / stars.pixel_grid.shape[0]

    def build_vis_matrix(self, stars, num_frames):

        '''

        :param stars: star system object, see class for attribute definitions
        :param step: (float) step by which to move the orbit phase to the next, default is in degrees
        :param deg: boolean indicating whether the step is in degrees
        :param by_frames: boolean indicating whether the step is being manully input or determined from a number of
        desired frames
        :param num_frames: the number of frames in the eclipse map
        :return: a 3d np.ndarray giving the shadow of the secondary star over the white dwarf region for a defined set
        of orbit phase points
        '''
        phi = 0
        gshape = stars.pixel_grid.shape

        step = 2 * np.pi / num_frames
        vis_matrix = np.zeros(shape=(gshape[0], gshape[1], num_frames))


        while phi <= 2*np.pi:

            frame_number = 0
            frame = stars.pixel_grid
            orb_x, orb_y, orb_z = stars.orb_x(phi), stars.orb_y(phi), stars.orb_z(phi)

            for i in tqdm(range(len(stars.pixel_grid[:, 0] - 1)), "Building Visibility Matrix"):
                for j in range(len(stars.pixel_grid[0, :] - 1)):

                    pxl_x, pxl_y, pxl_z = stars.pixel_orb_x(i, j), stars.pixel_orb_y(i,j), stars.pixel_orb_z(i,j)

                    #getting distance between pixel and secondary star in the yz-plane
                    pxl_2ndry_dist = np.sqrt((orb_y - pxl_y)**2 + (orb_z - pxl_z)**2)

                    '''checking if the distance between the pixel and secondary star is less than or equal too 
                    the secondary star radius'''

                    if pxl_2ndry_dist <= stars.Radius_2nd:

                        frame[i,j] = 0

                    elif pxl_2ndry_dist >= stars.Radius_2nd:
                        frame[i,j] = 1

                    else:
                        print("Something went wrong checking if pixel, ", '[',i,' ',j,'], ',
                              'was visible at phase: ', phi)
                        sys.exit()

            vis_matrix[:, :, frame_number] = frame



            phi += step

        stars.vis_matrix = vis_matrix

    def build_default_map(self, stars):
        pass
        

class Weights:

    def __init__(self):
        pass

    def covolution(self, K, J, delta):
        '''
        :param K: (list) a reference pixels coordinates.. should these be in orbit space?
        :param J: (list) the iterable pixel's coordinates... in orbit space?
        :param delta: (float) some width value (mostly arbitrary)
        :return: an weight factor
        '''

        dist_KJ = np.sqrt((K[0] - J[0])**2 + (K[1] + J[1])**2)

        wjk = np.exp(-dist_KJ**2/2*delta**2)

        return wjk


    def radial(self, K, J, delta):

        '''
        :param K: (list) a reference pixels coordinates.. should these be in orbit space?
        :param J: (list) the iterable pixel's coordinates... in orbit space?
        :param delta: (float) some width value (mostly arbitrary)
        :return: an weight factor
        '''

        dist_center_k = np.sqrt(K[0]**2 + K[1]**2)

        dist_center_j = np.sqrt(J[0]**2 + J[1]**2)

        wjk = np.exp( -(dist_center_k - dist_center_j)**2/2 * delta**2)

        return wjk


system = Star_System(1,1,1,1,1,1)

imager = Imager(10)

imager.assign_grid(system)
imager.build_vis_matrix(system, 500)

print(system.vis_matrix.shape)

































def main():
    some_system = Star_System(mass_ratio= 0.5, period = 6,
                              semi_major_axis= 10, orbital_inclination= 0,
                              light_curve= [1,2,3,4,5,6], distance= 100,
                              eccentricity= 0.01 )

    imager = Imager(5)

    imager.assign_grid(some_system)

    for i in range(some_system.pixel_grid.shape[0]):
        print('x coordinate at', i , some_system.pixel_orb_x(i,0))




main()