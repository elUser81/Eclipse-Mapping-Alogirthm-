from astropy.table import Table
import numpy as np
from scipy.optimize import fsolve
import scipy.stats as st
from tqdm import tqdm
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

def relax(grid, stop_at):
    """
    :param grid: a grid with the intensity at phase phi distributed over it
    :param stop_at: the number of iterations
    :return: The same grid from 'grid' with averaged (blurred) pixels.
            Fancy talk for saying that this function makes the image blurry.
            Uses the 4 adjacent pixels for average (see below in the value: avg)
    """

    dim = len(grid[:,0])
    for k in tqdm(range(stop_at)):
        for i in range(len(grid[:,0]) - 1):
            for j in range(len(grid[0,:]) - 1):
                if i - 1 < 0:
                    continue
                if i + 1 > dim:
                    continue
                if j - 1 < 0:
                    continue
                if j + 1 >= dim:
                    continue

                avg = np.mean([grid[i,j], grid[i + 1, j], grid[i-1,j], grid[i,j+1], grid[i, j - 1]])
                grid[i,j] = avg


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


def build_weighted_grid(grid, dell):
    """
    :param grid: a seed grid  
    :param dell: the width of the gaussian curve (standard deviation)
    :return: a modified sead grid with an applied weight function wjk
                    modifiy this later so that wjk function is a parameter for this function
                    gunna debug and test this first
    """
    center = [(grid.shape[0])/2, (grid.shape[1])/2]
    D = np.zeros(shape = grid.shape)
    for i in tqdm(range(grid.shape[0])):
        for l in range(grid.shape[1]):
            J = [i,l]
            wjkIks = []
            wjks = []

            for m in range(grid.shape[0]):
                for n in range(grid.shape[1]):
                    K = [m,n]
                    wjk = w_jk(J,K, dell = dell, center = center)
                    wjks.append(wjk)
                    Ik = grid[K[0],K[1]]
                    wjkIk = wjk*Ik
                    wjkIks.append(wjkIk)

            Dj = (sum(wjkIks)/sum(wjks))
            D[J[0],J[1]] = Dj
    return D


def plot_things(a, inclination, e):
    theta = np.deg2rad(inclination)
    phase = np.linspace(0,(2*np.pi))
    b = np.sqrt((a ** 2) * (1 - (e * e)))
    c = a*e
    x = a * np.cos(phase) * np.cos(theta)
    y = b * np.sin(phase)
    z = a * np.cos(phase) * np.sin(theta)

    wdx, wdy, wdz = -c * np.cos(theta), 0 , -c* np.sin(theta)


    fig = plt.figure()
    ax = Axes3D(fig)
    ax.set_xlabel('X axis')
    ax.set_ylabel("Y axis")
    ax.set_zlabel('Z axis')
    ax.set_ylim3d(-1,1)
    ax.set_xlim3d(-1,1)
    ax.set_zlim3d(-1,1)
    ax.plot(xs = x, ys = y, zs = z)
    ax.plot([wdx],[wdy],[wdz], 'ro')
    plt.show()

def get_orbit_shape(a, e, inclination):
    """
    :param a: semi major axis
    :param e: ecentricty
    :return: the elipse formula for the orbit and c, the distnce from the center to the focal point
    """
    b = np.sqrt((a**2)*(1-(e*e)))
    c = a*e
    def orbit_pos(phase):
        y = np.cos(phase)*a
        x = np.sin(phase)*b
        z = - np.cos(phase) * np.sin(inclination)
        return x,y,z
    return orbit_pos, c

def scale(x,y,N,radius,c,mode):
    """
    :param x: x coordinate
    :param y: y coordinate
    :param N: pixel range from -N to N
    :param radius: radius of the secondary star
    :param c: distance from center to focii
    :param mode: a string designating which conversion to perform
    :return:
    """
    if mode == 'radiustoN':
        xn = (N/radius)*x
        yn = (N/radius)*y
        return xn, yn
    if mode == 'Ntoradius':
        xradius = (radius/N)*x
        yradius = (radius / N) * y
        return xradius, yradius
    if mode == "NtoP":
        xp = x + N
        yp = N - y
        return xp, yp
    if mode == "PtoN":
        xn = x - N
        yn = N - y
        return xn, yn
    if mode == "NtoA":
        len_p = radius/N
        xa = x*len_p + c
        ya = y*len_p
        return xa, ya
    if mode == "AtoN":
        len_p = radius / N
        xn = (x - c)/ len_p
        yn = y/len_p
        return xn, yn

    else:
        print("Please enter a valid mode")
        sys.exit()

def is_it_visible(orbcoord,pxcoord, second_star_radius, conv_factor, inclination, wdcoord, phase):
    if np.pi/2 < phase < 0.75 * np.pi:
        return 1
    else:
        xp, yp = pxcoord
        wdx, wdy = wdcoord[0], wdcoord[1]
        xL0 = conv_factor*xp + wdx
        zL0 = (conv_factor*yp + wdy)*np.sin(inclination)
        xL1, zL1 = orbcoord[0], orbcoord[2]
        #R0 = np.sqrt((xL1 - xL0)**2 + (zL1 - zL0)**2)
        R0 = np.sqrt((xL1 - xL0) ** 2)
        if R0 > second_star_radius:
            return 1
        if R0 < second_star_radius:
            return 0
        else:
            print("Something went wrong checking if a pixel was visible..")
            sys.exit()




def build_fractional_visibility_frame(grid, orbit_pos, phase, inclination, second_star_radius, conv_factor, c):
    frame = np.zeros(shape = grid.shape)
    inclination = np.deg2rad(inclination)
    #phase = np.deg2rad(phase)
    orbx , orby, orbz = orbit_pos(phase, inclination)
    orbcoord = [orbx, orby, orbz]
    wdcoord = [-c, 0]
    for i in range(frame.shape[0]):
        for j in range(grid.shape[1]):
            pxcoord = [(i - (frame.shape[0] - 1)/2),(((frame.shape[1] - 1)/2) - j)]
            one_or_zero = is_it_visible(orbcoord, pxcoord, second_star_radius, conv_factor, inclination, wdcoord, phase)
            frame[i,j] = one_or_zero
    return frame




def main():
    seed = gkern(50)
    #plt.imshow(seed, interpolation= 'none')
    weighted_seed = build_weighted_grid(seed,dell = 4)
    relax(weighted_seed, stop_at= 30)
    plt.imshow(weighted_seed, interpolation='none')
    plt.show()

def main2():
    a = 1
    e = 0.1
    inclination = 180
    fig = plt.figure()
    ax = Axes3D(fig)
    orbit_pos, c = get_orbit_shape(a, e)
    phi = np.linspace(0, 2*np.pi)
    xs, ys,zs = orbit_pos(phi, inclination)
    ax.plot(xs,ys,zs)
    plt.show()


def animate_fractional_visibility(grid_size,a,e,second_star_radius,inclination, gamma_factor, num_frames):
    grid = gkern(kernlen= grid_size)
    # phase = np.deg2rad(1)
    conv_factor = gamma_factor * a / grid.size
    orbit_pos, c = get_orbit_shape(a, e)
    phase = 0
    while phase <= 2*np.pi:
        if phase == 0:
            frame = build_fractional_visibility_frame\
                    (grid, orbit_pos, phase, inclination, second_star_radius, conv_factor,c)
            p = plt.imshow(frame)
            print("Step: ", phase)
            phase += 2*np.pi/num_frames
        else:
            frame = build_fractional_visibility_frame \
                (grid, orbit_pos, phase, inclination, second_star_radius, conv_factor,c)
            p.set_data(frame)
            plt.pause(0.03)
            p = plt.imshow(frame)
            print("Step: ", np.rad2deg(phase))
            phase += 2*np.pi/num_frames



        '''
        if phase == 0:
            frame = build_fractional_visibility_frame \
                (grid, orbit_pos, phase, inclination, second_star_radius, conv_factor)
            p = plt.imshow(frame)
        else:
            frame = build_fractional_visibility_frame \
                (grid, orbit_pos, phase, inclination, second_star_radius, conv_factor)
            p.set_data(frame)
            plt.pause(0.2)'''




def main3():
    animate_fractional_visibility(grid_size=85,a=1,e=0,second_star_radius= 0.0003,
                                  inclination =0, gamma_factor = 1, num_frames= 100)


def main4():
    plot_things(a = 1, inclination= 20, e = 0.8)



main4()
