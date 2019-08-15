from astropy.table import Table
import numpy as np
from scipy.optimize import fsolve
from tqdm import tqdm
import sys
import matplotlib.pyplot as plt
import astropy.units as u
import os
cwd = os.getcwd()

'''
this program is based of the paper, 
Eclipse Mapping:
Astrotomography of Accretion Discs
    by Raymond Baptista
    
 links to two versions of the paper:
https://arxiv.org/pdf/1508.03067.pdf
http://cds.cern.ch/record/466802/files/0009472.pdf
'''


def make_pixels(N):
    #makes a table of pixels with coordininates of a grid with resolution N^2
    #and range -N to N
    grid = np.mgrid[-N:N+1, -N:N+1]
    xs = grid[0].flatten()
    ys = grid[1].flatten()
    fill = [0 for fill in xs]

    return xs, ys, fill

def get_func(q ,R):
    '''
    :param q: mass ratio
    :param R: orbital radius
    :return: a function to solve for the distance to the secondary star to the first lagrange point,
            which we will use as the radius of the secondary star
    '''
    def func(r):
        x = (q/(r*r)) + 1/(R*R) - (r/(R*R*R))*(1+q) - 1/((R-r)*(R-r))
        return x
    return func


def wjk(x1, y1, x2, y2, delta = 1):
    '''
    :param x1: x coordinate for target pixel
    :param y1: y coordinate for target pixel
    :param x2: x coordinate for comparison pixel
    :param y2: y coordinate for comparison pixel
    :param delta: some weight parameter (not really sure how it effects the image yet)
    :return: the weight value wjk (see raymond baptistas paper on eclipse mapping)
    '''
    d = np.sqrt((x2-x1)**2 + (y2 - y1)**2)
    wjk = np.exp(-(d**2)/(2*(delta**2)))
    return wjk


def rand_dist_brightness(brightness, num_pixels):
    '''
    :param brightness: the intensity at some phase of the lightcurve
    :param num_pixels: the number of pixels in the grid
    :return: a grid where the brightness is randomly distributed throughout the map
    '''
    r = np.random.random_sample((num_pixels,))
    s = sum(r)
    r = [(i / s) * brightness for i in r]
    #scaler = MinMaxScaler( feature_range= (0,1))
    grid = np.array(r)
    #grid.reshape(1,-1)
    #grid = scaler.fit_transform(grid)
    shape = int(np.sqrt(num_pixels))
    #np.reshape(grid, newshape=(shape,shape))
    grid.resize((shape,shape))
    #grid = normalize(grid, axis = 0)

    return grid

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

def get_orbit_shape(a, e):
    """
    :param a: semi major axis
    :param e: ecentricty
    :return: the elipse formula for the orbit
    """
    b = np.sqrt((a**2)*(1-(e*e)))
    c = a*e
    def orbit_pos(phi):
        x = np.cos(phi)*a
        y = np.sin(phi)*b
        return x,y
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


def frac_vis(grid, phi, orbit_shape, radius,Inc):
    x1, y1 = orbit_shape(phi)
    for i in range(len(grid[:, 0]) - 1):
        for j in range(len(grid[0, :]) - 1):
            d = np.sqrt(((x1-i)*np.cos(Inc))**2 + (y1-j)**2)
            if d <= radius:
                grid[i, j] = 0

def build_table(LC_data):

    data = Table([LC_data[:,0], LC_data[:,1]], names = ("phase", "Ip"))
    return data


def get_frac_vis(phi, I, n_range, N, RL1,c, orbit_pos, i):
    n_visible = 0
    orb_x, orb_y = orbit_pos(phi)
    for x in n_range:
        for y in n_range:
            x, y = scale(x ,y, N, RL1, c, mode = "NtoA")
            d = np.sqrt((orb_x - x*np.sin(i))**2 + (orb_y - y)**2)




def dist_brightness(pixels, I, phi, orbit_pos, radius, N, c, theta ):
    orb_x, orb_y = orbit_pos(phi)
    for i in range(len(pixels['Ip'])):
        x, y = float(pixels['X'][i]), float(pixels['Y'][i])
        x, y = scale(x,y,N,radius,c, mode = 'NtoA')

        d = np.sqrt(((orb_x - x)*np.sin(theta))**2 + (orb_y - y)**2)
        N0 = pixels['num_avgs'][i]
        A0 = pixels["Ip"][i]

        if np.cos(phi) > 0:
            if d < radius:
                I = 0

        A = ((A0 * N0) + I) / (N0 + 1)
        pixels["Ip"][i] = A
        pixels['num_avgs'][i] += 1




def build_grid_from_pixels(pixels, n_dim):
    arr = np.array(pixels["Ip"])
    arr.resize((n_dim,n_dim))
    return arr


def get_phasemin(data):
    I_min = np.amin(data[:,1])
    indx = np.where(data[:,1] == I_min)[0][0]
    phasemin = data[:,0][indx]
    return phasemin



def e_map(N,  a, e, q , theta, light_curve_data, est = 1):

    #make_pixels creates a square grid from -N to N
    xs, ys, fill = make_pixels(N)
    index = [i for i in range(len(xs))]
    #pixels has dimensions (N+2)^2
    pixels = Table([fill, xs, ys, fill, index], names = ("Ip", 'X', 'Y', 'num_avgs', 'index'))
    num_pixels = len(pixels['X'])
    RL1_func = get_func(q, a)
    radius = fsolve(RL1_func, np.ndarray([est]))[0]
    rL1 = a - radius
    A = (2*rL1)**2
    Ap = A/num_pixels
    n_dim = np.sqrt(num_pixels)
    orbit_pos, c = get_orbit_shape(a,e)

    data = build_table(light_curve_data).group_by('phase')
    n_range = [i for i in range(-N,N+1)]
    phasemin = get_phasemin(light_curve_data)
    phase0 = data['phase'][0]
    for i in tqdm(range(len(data['phase'])), desc = 'Mapping...'):
        I = data['Ip'][i]
        phase = data['phase'][i]
        phi = 2*np.pi*(phase + phase0 - phasemin)
        if phi <= 0:
            continue
        if phi >= np.pi/4:
            break
        dist_brightness(pixels,I,phi,orbit_pos,radius,N,c,theta)
        #if i >= 20:
            #break
    pgrid = build_grid_from_pixels(pixels, int(n_dim))

    print(pgrid)
    print(pixels)
    fig, ax = plt.subplots()
    im = ax.imshow(pgrid)
    plt.show()








def main():
    path = os.path.join(cwd, 'Data', 'ktwo201903318-c01-phase.txt')
    light_curve = np.loadtxt(path)
    e_map(N = 10, a = 5, e = 0, q = 0.55, theta = 76, light_curve_data= light_curve)




main()