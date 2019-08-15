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


def gkern(kernlen=21, nsig=3):
    """Returns a 2D Gaussian kernel array."""

    interval = (2*nsig+1.)/(kernlen)
    x = np.linspace(-nsig-interval/2., nsig+interval/2., kernlen+1)
    kern1d = np.diff(st.norm.cdf(x))
    kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
    kernel = kernel_raw/kernel_raw.sum()
    return kernel


def w_jk(j,k, dell, center):
    """
    :param j: a list or tuple containing the coordinates for the jth pixel
    :param k: a list or tuple containing the coordinates for the kth pixel
    :param dell: standard deviation (arbitrarily assigning it for now, will probably calculate it better later)
    :param center: coordinates for the center of the grid
    :return: a weight factor to create the default ecplise map
    """

    xj, yj, xk, yk, xc, yc = j[0], j[1], k[0], k[1], center[0], center[1]

    Rj = np.sqrt((xj - xc)**2 + (yj - yc)**2)
    Rk = np.sqrt((xk - xc)**2 + (yk - yc)**2)

    wjk = np.exp((-(Rj - Rk)**2)/(2*(dell**2)))
    return wjk


def get_orbit_shape(a, e, inclination):
    """
    :param a: semi major axis
    :param e: ecentricty
    :return: the elipse formula for the orbit and c, the distnce from the center to the focal point
    """
    b = np.sqrt((a**2)*(1-(e*e)))
    c = a*e
    def orbit_pos(phase):
        '''
        :param phase: phase in radians of the orbit
        :param inclination: inclination of orbital plane (in radians)
        :return: the list of x, y, z coordinates at any phase point
        '''
        x = a * np.cos(phase) * np.cos(inclination)
        y = b * np.sin(phase)
        z = - a * np.cos(phase) * np.sin(inclination)

        return [x,y,z]
    return orbit_pos, c




def build_fractional_visibility_matrix_test(a,R,N,e,inclination, num_frames, second_star_radius, gamma = 1):
    '''
    :param a: semi-major axis
    :param R: length of the pixel grid in orbit space
    :param N: number of pixels along the axis. Resolution will be N*N
    :param e: eccentricity of the orbit
    :param inclination: orbital inclination
    :param second_star_radius: radius of the second star
    :param gamma: a extra factor I multiply by when converting to orbit space. Arbitrarily/experimentally assigned, default is 1
    :return: a 3d matrix of ones and zeros, denoting whether a given pixel is visible from the yz plane at any point in the orbit
    '''

    if N%2 == 0:
        N -= 1
    matrix = np.zeros(shape = [N,N])
    inclination = np.deg2rad(inclination)
    #getting orbit position function
    orbit_pos, c = get_orbit_shape(a,e, inclination)
    #coordinates for the white dwarf at foci
    wdcoord = [-c * np.cos(inclination), 0 , -c * np.sin(inclination)]
    conv_factor = gamma * (R/N)
    phase = 0
    # this will probably be a for loop when I use actual lightcurve data, which have specific phase points
    # phase will be given in MJD so you will have to find the phase angle from the period and change in time
    print("Building fractional visibility matrix.....")
    while phase <= 2*np.pi:
        frame = np.zeros(shape = [N,N])
        # [x,y,z] orbit coordinates
        orbit_coords = orbit_pos(phase)
        for i in range(N):
            for j in range(N):
                # converting [x,y] grid coordinates so that center (0,0)
                # is in the center of the grid, not in the top left
                xp = (i - ((N-1)/2))
                yp = (((N-1)/2) - j)

                # converting from pixel space to orbit space, centering around white dwarf,
                #  and rotating about y axis based on inclination
                xp_orb = (conv_factor * xp - c) * np.cos(inclination)
                yp_orb = conv_factor * yp
                zp_orb = -(conv_factor * xp - c) * np.sin(inclination)
                p_orbcoords = [xp_orb, yp_orb, zp_orb]
                del xp_orb, yp_orb, zp_orb

                # checking if the second star is behind the white dwarf, if so, the whole grid is visible
                if np.pi / 2 < phase < 0.75 * np.pi:
                    one_or_zero = 1
                else:
                    one_or_zero = is_it_visible(p_orbcoords, orbit_coords, second_star_radius)
                frame[i,j] = one_or_zero
        matrix = np.vstack((matrix, frame))
        phase += 2*np.pi / num_frames

    matrix = np.reshape(matrix, (N, N, (num_frames+1)))
    matrix = np.delete(matrix, 0, 2)
    print("Done. Matrix shape is: ", matrix.shape)
    return matrix

def is_it_visible(p_orbcoords, orbit_coords, second_star_radius):
    '''
    :param p_orbcoords: [x,y,z] coordinates of each pixel in the orbit space
    :param orbit_coords: [x,y,z]  coordinates of the 2nd star orbit
    :return: 1 or 0 if the pixel is visible or invisible
    '''

    xp, yp, zp = p_orbcoords[0], p_orbcoords[1],p_orbcoords[2]
    xstar, ystar, zstar = orbit_coords[0], orbit_coords[1], orbit_coords[2]

    #checking distances in the yz plane
    R0 = np.sqrt((yp - ystar)**2 + (zp - zstar)**2)
    if R0 < second_star_radius:
        return 0
    elif R0 > second_star_radius:
        return 1
    else:
        print("Something went wrong checking if a pixel was visible...")
        sys.exit()



def animate_fac_vis_matrix(matrix, stop = None):
    num_frames = matrix.shape[2]
    for i in range(num_frames):
        if i == 0:
            frame = matrix[:,:,i]
            p = plt.imshow(frame)
            print("Step: 0")
        elif i == stop:
            break
        else:
            frame = matrix[:,:,i]
            p.set_data(frame)
            plt.pause(0.03)
            p = plt.imshow(frame)
            print("Step: ", np.rad2deg(2*np.pi/num_frames*i))


def main():
    matrix = build_fractional_visibility_matrix_test(a = 1, R = 1, N = 51, e = 0, num_frames= 50,
                                                     second_star_radius= 0.05, inclination= 50)

    animate_fac_vis_matrix(matrix)

main()



























































