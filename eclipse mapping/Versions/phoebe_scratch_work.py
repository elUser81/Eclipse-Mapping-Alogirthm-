
import numpy as np







#the masses are constrained by the period and semi-major axis. Hence, you can't define them, they are always calculated.
#BUT...you can call the calculated value like you call everything else..
#print(b['mass@primary@star@component'])


#params = {'q': 0.21, 'period': 0.09112, 'sma': 0.85, 'ecc':0.0,
 #       'distance': 7.1428571428571, 'incl': 88}

import sys


'''
def __init__(self, system_name, incl, period, sma, q, ecc, N, distance, chi_AIM, delta,
                 light_curve_data = None, radius_secondary = None):

'''




class Binary:

    def __init__(self, orbital_params, resolution, chi_AIM, delta,
                 light_curve_data, method, radius_secondary_given = True, fluxes_not_mags = True):

        '''
        :param orbital_params: (dict) a dictionary of orbital parameters, see Processes.make_parameter_set

        :param method: (str) the method by which the algorithm will solve the multivariable entropy funciton
                             see scipy.optimize.minimize docs

        :param resolution: (int) the dimensions of the pixel grid, so reslution = 25 gives a 25 by 25 pixel grid.
                                 renamed to self.N in __init__ .. I know, I know... but I decided to name the parameter
                                 'resolution' instead of 'N' too late...

        :param chi_AIM: (float) a chi_squared value that determines how constrained the accretion image will be to the
                                observed data, typically a chi value thats a little bit more than the average uncertainty
                                of the observed data is best so the algorithm doesn't treat noise as real data

        :param delta: (float) the spread of the guassian initial guess, also effects the spread of the weight matrix.

        :param light_curve_data: (np.array) a 2d array of the observed light curve data.. organize such that
                                  light_curve_data[0] = observation times (days), light_curve_data[1] = flux,
                                  light_curve_data[2] = uncertainties

        :param radius_secondary_given: (bool) tells Binary object whether secondary radius is given.. if not given,
                                        binary object will estimate the radius with the roche lobe equation
                                        ***Not yet implemeted****

        :param fluxes_not_mags: (bool) tells Binary object that input brightness data is flux, not magnitude
                                Processes.build_compute_save will normalize the data if it is True.
        '''


        import phoebe
        from scipy.optimize import minimize
        #building phoebe model inside binary class

        #building a dictionary of all parameters passed in the object, will use later for logging...
        self.all_params = orbital_params
        self.all_params['resolution'] = resolution
        self.all_params['chi_AIM'] = chi_AIM
        self.all_params['delta'] = delta
        self.all_params['radius_secondary_given'] = radius_secondary_given
        self.all_params['method'] = method

        b = phoebe.default_binary()
        b.set_value('q', orbital_params['q'])
        b.set_value('incl@binary@orbit@component', orbital_params['incl'])
        b.set_value('period@binary@orbit@component', orbital_params['period'])
        b.set_value('sma@binary@orbit@component', orbital_params['sma'])
        b.set_value('ecc@binary@orbit@component', orbital_params['ecc'])
        b.set_value('distance@system', orbital_params['distance'])
        self.mass_primary = b.get_value('mass@primary@star@component')
        self.mass_secondary = b.get_value('mass@secondary@star@component')
        self.sma_primary = b.get_value('sma@primary@star@component')
        self.sma_secondary = b.get_value('sma@secondary@star@component')
        b.add_dataset('orb')


        #b.run_compute()
        #defining values.. so many values
        self.system_name = orbital_params['system_name']
        self.radius_primary = orbit_params['radius_primary']
        self.fluxes_not_mags = fluxes_not_mags
        self.method = method
        self.delta = delta
        self.chi_AIM = chi_AIM
        self.distance = orbital_params['distance']
        self.period = orbital_params['period']
        self.binary_sma = orbital_params['sma']
        self.incl = np.deg2rad((90 - orbital_params['incl'])) #orbital inclination, 90 - incl becase 90 degrees is considered perpendicular
        self.ecc = orbital_params['ecc']
        self.q = orbital_params['q']
        self.N = resolution #pixel grid dimension
        self.accretion_grid = np.ndarray(shape = (self.N,self.N))
        self.light_curve_data = light_curve_data #should be an array of shape (3, number_of_datapoints)
                                                 # axis [0] is times, [1] is fluxes [2] is uncertainties
        self.JDtimes = light_curve_data[:, 0]
        self.observed_fluxes = light_curve_data[:, 1]
        self.observed_uncertainties = light_curve_data[:, 2]

        '''CALL self.construct_phasepoints() BEFORE use'''
        self.constructed_phasepoints = None

        '''CALL self.get_synthetic_fluxes BEFORE use'''
        self.synthetic_fluxes = None


        self.delta = 3 ###going to make this a variable parameter later

        #creating space to compute orbital models

        #self.primary_orbit, and self.secondary orbit return a list of [x, y, z] coordinates
        #easier for handling when building the eclipse map

        '''***MUST CALL self.build_orbital_model BEFORE using these methods***'''
        self.primary_orbit = None
        self.secondary_orbit = None

        #these functions return the individual points on each axis, easier for parametric graphing
        self.x_primary = None
        self.y_primary = None
        self.z_primary = None

        self.x_secondary = None
        self.y_secondary = None
        self.z_secondary = None

        '''***MUST CALL self.build_fractional_visibility_matrix BEFORE using these attributes***'''
        self.image_entropy = None
        self.image_solution = None

        #gotta solve for RL1 soooooo.... use fsolve fromm scipy? (not working)  ---> play the system, make the equation
        #an absolute value function and minimize for 0 ---- ayyyyyooooo

        #python wont let me call self.attribute in a method definition inside __init__ so i'm renaming values


        def to_kg(sol_mass):
            #converts_solar mass to kg
            #1.9889200011445836e+30

            return sol_mass * 1.9889200011445836 * 10**30


        M1, M2 = map(to_kg, [self.mass_primary, self.mass_secondary])
        sma = self.binary_sma * 1.496 * 10**11

        #1 AU is 1.496e+11 m

        def get_RL1_func(M1, M2, sma):
            #defining args of the 1st lagrange point function

            def RL1_func(r):
                #this equation gives the roots of the RL1 function measured from the 2nd star. it is set equal to zero
                #takning the absolute value and minimizing for 0 and taking the r value that gives me 0 because fsolve i
                #is being dumb

                return abs((M2 / (r ** 2)) + M1 / (sma ** 2) - ((r * (M1 + M2)) / (sma ** 3)) - (M1 / ((sma - r) ** 2)))  # == 0

            return RL1_func

        #getting estimate value to start fsolve and building lagrange point function
        guesstimate = np.asarray([0])
        RL1_func = get_RL1_func(M1, M2, sma)

        #solving for RL1..

        RL1_from_secondary = minimize(RL1_func, x0= guesstimate, method= 'Nelder-Mead').x
        RL1_from_primary = sma - RL1_from_secondary

        self.RL1_from_primary = RL1_from_primary[0]/(1.496 * 10**11) #back to solar radii
        self.RL1_from_secondary = RL1_from_secondary[0]/(1.496 * 10**11)

        self.pixel_to_RL1 = self.RL1_from_primary/self.N #conversion factor, converts pixel units to RL1 units

        if radius_secondary_given:
            self.radius_secondary  = orbital_params['radius_secondary']
        else:
            #self.radius_secondary = None
            print 'Roche Lobe Equation estimate not yet implemented, please add a secondary radius to your orbital ' \
                  'parameter set'
            sys.exit()

        #making room for visibility matrix
        self.visibility_matrix = None

        #making room for temperature solutions
        self.temperature_solutions = None



    def to_RL1(self, a_distance):
        #converts a distance in solar-radius to units of RL1 from primary star

        return a_distance/self.RL1_from_primary


    def build_orbital_model(self, phi0 = 0):
        #creates two position functions for the primary and secondary star orbits
        #at some phase phi, rotates about the y-axis at some inclination 'incl'
        #and assigns them to self.primary_orbit, and self.secondary orbit

        # I would like this to be in __init__ but it isnt cooperating ....

        #also creates indivudual axis funtctions and assigns them to self.x1 ,x2, y1, y2, z1, z2

        from numpy import sin, cos

        a1 = self.sma_primary #semi-major axis of primary orbit
        a2 = self.sma_secondary #semi-major axis of secondary orbit
        incl = self.incl #orbital inclination
        RL1_from_primary = self.RL1_from_primary


        def primary_orbit(phi):

            #returns the x-y-z coordinates as a list

            x = -a1 * cos(phi + phi0) * cos(incl)/RL1_from_primary
            y = -a1 * sin(phi + phi0)/RL1_from_primary
            z = -a1 * cos(phi + phi0) * sin(incl)/RL1_from_primary

            return [x, y, z]

        def secondary_orbit(phi):

            x = a2 * cos(phi + phi0) * cos(incl)/RL1_from_primary
            y = a2 * sin(phi + phi0)/RL1_from_primary
            z = a2 * cos(phi + phi0) * sin(incl)/RL1_from_primary

            return [x, y, z]


        self.primary_orbit = primary_orbit
        self.secondary_orbit = secondary_orbit

    #theres probably a better way to do this but oh well... I gotta get this working
    #I know I could use decorators, but id rather just copy and paste /RL1_from_primary at this point

        def x1(phi):
            return -a1 * cos(phi + phi0) * cos(incl)/RL1_from_primary

        def y1(phi):
            return -a1 * sin(phi + phi0)/RL1_from_primary

        def z1(phi):
            return -a1 * cos(phi + phi0) * sin(incl)/RL1_from_primary

        def x2(phi):
            return a2 * cos(phi + phi0) * cos(incl)/RL1_from_primary

        def y2(phi):
            return a2 * sin(phi + phi0)/RL1_from_primary

        def z2(phi):
            return a2 * cos(phi + phi0) * sin(incl)/RL1_from_primary

        self.x_primary = x1
        self.y_primary = y1
        self.z_primary = z1

        self.x_secondary = x2
        self.y_secondary = y2
        self.z_secondary = z2


    def pixel_orb_coords(self, xp, yp, phi):
        '''
        ***MUST CALL self.build_orbital_model() BEFORE use****
        :param xp: (int) x coordinate in pixel space
        :param yp: (int) y coordinate in pixel space
        :param phi: (float) phase of the orbit
        :return: the x,y,z coordinates of the pixel in orbit space and phase phi as a list in units of RL1

        '''

        from numpy import sin, cos

        #factor that centers the coordinates of the pixel grid around 0,0
        centered = (self.N -1)/2

        #centering xp and yp:
        xg = (xp - centered) ##changed from xp - centered##
        yg = (centered - yp) ###changed from centered - yp###

        #converting pixel to orbit space

        x_space = xg * self.pixel_to_RL1
        y_space = yg * self.pixel_to_RL1

        #now to rotate the pixels about the y axis

        x_space_abt_y = x_space * cos(self.incl)
        y_space_abt_y = y_space
        z_space_abt_y = x_space * sin(self.incl)

        #now to rotate about relative z axis:

        x_space_abt_y_z = x_space_abt_y * cos(phi) - y_space_abt_y * sin(phi)
        y_space_abt_y_z = x_space_abt_y * sin(phi) + y_space_abt_y * cos(phi)
        z_space_abt_y_z = z_space_abt_y

        xwd, ywd, zwd = self.primary_orbit(phi) #pulling coordinates for primary star at current phase
        # --> xwd = x coordinate of the white dwarf (primary star)

        x_final, y_final, z_final = (x_space_abt_y_z + xwd), (y_space_abt_y_z + ywd), (z_space_abt_y_z + zwd)

        return [x_final, y_final, z_final]

    @staticmethod
    def invert(frame):
        for x in range(frame.shape[0]):
            for y in range(frame.shape[1]):

                if frame[x, y] == 0:

                    frame[x, y] = 1

                elif frame[x, y] == 1:
                    frame[x, y] = 0

                else:
                    print('A frame value was neither 1 or 0 when inverting... what.. ')
                    sys.exit()

        return frame

    def gkern(self):
        import scipy.stats as st
        """Returns a 2D Gaussian kernel array."""
        nsig = self.delta
        kernlen = self.N
        interval = (2 * nsig + 1.) / (kernlen)
        x = np.linspace(-nsig - interval / 2., nsig + interval / 2., kernlen + 1)
        kern1d = np.diff(st.norm.cdf(x))
        kernel_raw = np.sqrt(np.outer(kern1d, kern1d))
        kernel = kernel_raw / kernel_raw.sum()
        return kernel

    def build_visibility_matrix(self):
        from numpy import amin, where, pi, ndarray, sqrt
        '''***MUST CALL self.build_orbital_model() BEFORE use'''


        # also need to implement roche lobe equation.... right now, assuming radius will be given

        #assuming the data structure will be an array, so I'm treating it that way for now, will change syntax when
        #I start integrating actual data.


        num_frames = len(self.constructed_phasepoints)

        visibility_matrix = ndarray(shape = (self.N, self.N, num_frames))


        for i in range(num_frames):

            matrix_frame = ndarray(shape=(self.N, self.N))

            #pulling phase point from list of phase points and  converting to radians
            phi = self.constructed_phasepoints[i] * 2 * pi

            #getting secondary star coordinates and putting them in units of RL1 from primary
            x_2nd, y_2nd, z_2nd = map(self.to_RL1, self.secondary_orbit(phi))

            for x in range(self.N): #range or range(len())??
                for y in range(self.N):  # range or range(len())??

                    x_px, y_px, z_px = self.pixel_orb_coords(x, y, phi)

                    #getting the y-z projection of the distance between the orbit
                    y_z_distance = sqrt((y_2nd - y_px)**2 + (z_2nd - z_px)**2)

                    secondary_radius_in_RL1 = self.radius_secondary/self.RL1_from_primary

                    if y_z_distance <= secondary_radius_in_RL1:

                       if phi >= 0.5 * pi and phi <= 1.5 * pi:

                        matrix_frame[x, y] = 1

                       if phi >= 1.5 * pi or phi <= 0.5 * pi:
                        matrix_frame[x, y] = 0


                    else:
                        matrix_frame[x, y] = 1

            visibility_matrix[:, :, i] = matrix_frame

        self.visibility_matrix = visibility_matrix


    def animate_visibility_matrix(self):
        import os
        import matplotlib.pyplot as plt
        import matplotlib.animation as animation
        from numpy import rad2deg

        visibility_matrix = self.visibility_matrix
        num_frames = visibility_matrix.shape[2]


        fig = plt.figure()

        ims = []
        for i in range(num_frames):
            frame = visibility_matrix[:, :, i]
            im = plt.imshow(frame, animated = True)

            ims.append([im])

        ani = animation.ArtistAnimation(fig, ims, interval = 50, blit = True,
                                            repeat_delay= 1000)

        target = ['animation_','incl_', str(rad2deg(90 - self.incl)),'.html']
        cwd = os.getcwd()
        path = os.path.join(cwd, 'Outputs', 'eclipse_map_animations', ''.join(target))
        ani.save('movie.html')

    @staticmethod
    def dist_from_center(coords):
        ##pretty sure these should be spatial coordinates, not pixel coordinates.. debug this later
        #its only called in self.get_wjks, so you should only have to add the conversion on at the end below
        x, y = coords[0], coords[1]
        dist = np.sqrt(x ** 2 + y ** 2)
        return dist


    def build_default_map(self, wjk):
        '''

        :param I: list of all the brightnesses. shape = (n**2)
        :param wjk: weight matrix flattened into 2d with shape = (n**2, n**2)
        :return: list of the weighted default map of the brightness distribution
        '''
        I = self.image_solution.reshape(int(self.N**2))
        D = np.asarray([sum(a * b for a, b in zip(I, wjk[j, :])) / sum(wjk[j, :]) for j in range(len(I))])

        return D

    def get_wjks(self):
        from tqdm import tqdm

        N = self.N
        delta = self.delta

        wjk_matrix = np.ndarray(shape=(N, N, int(N ** 2)))

        frame_number = 0

        # these i and j loops define the jth pixel for the weight function
        for i in tqdm(range(N), ''.join(["Building weight matrix for ", self.system_name])):
            for j in range(N):
                # centering grid around coordinates [0,0]
                xj = i - ((N - 1) / 2)
                yj = ((N - 1) / 2) - j

                jframe = np.ndarray(shape=(N, N))

                Rj = self.dist_from_center([xj, yj])


                # these k and l loops define the kth pixel for the weight function
                # because these loops are inside the i and j loops, k and l will complete a full cycle for every i and j

                for k in range(N):
                    for l in range(N):
                        xk = k - ((N - 1) / 2)
                        yk = ((N - 1) / 2) - l

                        Rk = self.dist_from_center([xk, yk])

                        wjk = np.exp(-(Rj - Rk) ** 2 / (2 * delta ** 2))

                        jframe[k, l] = wjk

                wjk_matrix[:, :, frame_number] = jframe

                frame_number += 1

        return wjk_matrix


    @staticmethod
    def get_jacobian(f, dx = 10 ** -8):

        def jacobian(I):
            n = len(I)
            func = f(I)
            jac = np.zeros((n, n))
            for j in range(n):  # through columns to allow for vector addition
                Dxj = (abs(I[j]) * dx if I[j] != 0 else dx)
                I_plus = [(Ii if k != j else Ii + Dxj) for k, Ii in enumerate(I)]
                jac[:, j] = (f(I_plus) - func) / Dxj
            return jac
        return jacobian

    @staticmethod
    def normalize(array, over = None):
        from numpy import asarray, amax
        if not over:
            return asarray(map(lambda x: x/amax(array), array))
        else:
            return asarray(map(lambda x: x/over, array))

    def get_entropy_function(self):
        N = self.N
        wjks = self.get_wjks()
        wjks = wjks.reshape((int(N ** 2), int(N ** 2)))

        def entropy_function(I):
            D = np.asarray([sum(a * b for a, b in zip(I, wjks[j, :])) / sum(wjks[j, :]) for j in range(len(I))])

            q = np.asarray([i / sum(D) for i in D])
            p = np.asarray([i / sum(I) for i in I])
            # actual function is - sum(....  but im minimizing the positive version to find the maximum of its opposite
            S = sum([pj * np.log(pj / qj) for pj, qj in zip(p, q)])
            return S


        return entropy_function



    def construct_phasepoints(self):
        '''
        :return: a list of phase points assigned to self.phase_points
        '''
        from numpy import amin, where, asarray

        #fluxes = asarray(self.observed_fluxes)
        jdtimes = asarray(self.JDtimes)

        #assuming lowest flux wil give the time of super conjunction
        #time_of_superconj = jdtimes[where(fluxes == amin(fluxes))]
        #ephemeris = jdtimes[0]

        #taking the difference between time of super conjunction and ephemeris to find the phase angle at ephemeris.
        #roll_back = time_of_superconj - ephemeris
        old_phase_point = 0
        constructed_phases = []

        for i in range(len(jdtimes)):

            #was old_phase_point - roll_back
            new_phase_point = (old_phase_point)
            constructed_phases.append(new_phase_point)

            if i + 1 >= len(jdtimes):
                break

            delta_phase = (jdtimes[i + 1] - jdtimes[i])/self.period
            old_phase_point += delta_phase

        self.constructed_phasepoints = asarray(constructed_phases)

    def estimate_chi_sqrd(self):
        from numpy import mean, amax, amin
        M = len(self.observed_uncertainties)
        sig_max_sqrd = amax(self.observed_uncertainties)
        chi_sqrd_estimate = sum(map(lambda x: sig_max_sqrd/(x**2), self.observed_uncertainties))/M
        return chi_sqrd_estimate


    def get_synthetic_fluxes2(self):
        '''
        :param I: the solution set of brightneses
        :return: finds synthetic fluxes with the given maximum entropy solution and assigns them to
                    self.synthetic fluxes
        '''

        '''***there should ALREADY be an image solution BEFORE this method is called, call self.get_image_solution FIRST***'''
        #phi might have to be a frame slice instead of a phase value.. but I can find said slice with the phase value

        from numpy import cos, asarray
        theta_sqrd = ((self.RL1_from_primary/self.distance) ** 2) * cos(self.incl)

        phase_points = self.constructed_phasepoints
        Is = self.image_solution.reshape((int(self.N**2)))

        #reshaping eclipse map (visibility matrix) to 2d with dimensions N**2, length of list of light_curve_data
        flattened_matrix = self.visibility_matrix.reshape((self.N**2, len(phase_points)))

        fluxes = []

        for phi in range(len(phase_points)):
            #getting eclipse map for phase phi
            V_phi = flattened_matrix[:, phase_points[phi]]

            I_times_V = [I * V for I, V in zip(Is, V_phi)]
            flux = (theta_sqrd/4 * self.N**2) * sum(I_times_V)

            fluxes.append(flux)

        self.synthetic_fluxes = asarray(fluxes)


    def get_synthetic_fluxes(self):
        from numpy import cos, where, asarray

        phase_points = self.constructed_phasepoints
        N = self.N
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)
        flattened_matrix = self.visibility_matrix.reshape((N ** 2, len(phase_points)))
        chi_AIM = self.chi_AIM
        I = self.image_solution.reshape((N**2))

        def synth_flux(I, phi):
            index = where(phase_points == phi)
            V_phi = flattened_matrix[:, index]
            I_times_V = [i * v for i, v in zip(I, V_phi)]

            flux = (theta_sqrd / 4 * N ** 2) * sum(I_times_V)
            return flux

        synth_fluxes = [synth_flux(I, phi) for phi in phase_points]

        self.synthetic_fluxes = asarray(synth_fluxes).reshape((len(self.observed_fluxes)))


    def build_chi_constraint(self):
        from numpy import cos, where, mean, amax

        N = self.N
        phase_points = self.constructed_phasepoints
        obs_fluxes = self.observed_fluxes
        sigmas = self.observed_uncertainties
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)
        flattened_matrix = self.visibility_matrix.reshape((N ** 2, len(phase_points)))
        chi_AIM = self.chi_AIM


        def chi_constraint(I):

            M = len(phase_points)

            def synth_flux(I, phi):
                index = where(phase_points == phi)
                V_phi = flattened_matrix[:,index]
                I_times_V = [i * v for i, v in zip(I, V_phi)]

                flux = (theta_sqrd / 4 * N ** 2) * sum(I_times_V)
                return flux


            synth_fluxes = [synth_flux(I, phi) for phi in phase_points]


            synth_minus_observed_sqrd = [(synth - obs)**2 for synth, obs in zip(synth_fluxes, obs_fluxes)]
            sigma_sqrd = map(lambda x: x**2, sigmas)
            fractions = [num/denom for num, denom in zip(synth_minus_observed_sqrd, sigma_sqrd)]
            chi_sqrd = sum(fractions)/M

            return float(chi_sqrd - chi_AIM * (N**2))


        return chi_constraint


    def get_image_solution(self):
        from scipy.optimize import minimize
        from numpy import amax, asarray
        import time
        from numdifftools import Jacobian

        N = self.N
        S = self.get_entropy_function()
        guess = self.gkern().reshape((int(N)**2))
        chi_sqrd = self.build_chi_constraint()

        #, {'type':'ineq', 'fun': ineq}
        method = self.method
        cons = [{'type': 'eq', 'fun': chi_sqrd}]
        print"Solving entropy equation with resolution: ", self.N
        start = time.time()
        res = minimize(S, x0 = guess, method = method, constraints= cons)
        stop = time.time()

        solve_time = start - stop
        self.all_params['solve time'] = solve_time

        image_entropy = S(res.x)
        image_solution = res.x
        sol_max = amax(image_solution)
        normalized_solution = map(lambda x: x/sol_max, image_solution)


        self.image_entropy = image_entropy
        self.image_solution = asarray(normalized_solution).reshape((N,N))


    def normalize_fluxes(self):
        '''
        :return: a normalized data set of fluxes and uncertainties...only use on fluxes, not magnitudes
        '''

        from numpy import amax, asarray

        fluxes = self.observed_fluxes
        sigmas = self.observed_uncertainties

        flux_max = amax(fluxes)
        sig_max = amax(sigmas)

        norm_fluxes = asarray(map(lambda flux: flux/flux_max, fluxes))
        norm_sigmas = asarray(map(lambda sig: sig/sig_max, sigmas))

        self.observed_fluxes = norm_fluxes
        self.observed_uncertainties = norm_sigmas

    def flux_to_temp(self):
        from numpy import asarray, amax

        image_solutions = self.image_solution.reshape(self.N**2)
        #temp1 = self.primary_temperature # need to make this an attribute

        def temp_func(I):
            steph_boltz = 5.67 * 10**-8
            T = (I / steph_boltz) ** (1/4)
            return T

        temp_solutions = map(temp_func, image_solutions)

        #normalizing
        temp_solutions = asarray(map(lambda x: x/amax(temp_solutions), temp_solutions))
        #temp_solutions = asarray(map(lambda x: x * temp1, temp_solutions))

        self.temperature_solutions = temp_solutions.reshape((self.N, self.N))


    def save_bbody_curve(self, path):
        '''returns a list of radii and average azimuthal temperature to plot'''
        from numpy import mean, asarray, floor, log
        import matplotlib.pyplot as plt

        bb_curve = []
        N = self.N
        conv = self.pixel_to_RL1 #a conversion factor
        temp_solutions = self.temperature_solutions
        radius_primary = self.radius_primary
        radii = asarray([i for i in reversed(map(lambda x: x * conv, range(int(floor(N/2)))))])
        steph_boltz = 5.67 * 10 ** -8
        viscosity = 1 # will make an attribute later

        def dist(x,y):
            from numpy import sqrt
            #takes pixel coordinates, centers them around 0,0 and calculates the distance from the center
            centered = (N-1)/2
            x = x - centered
            y = centered - y
            d = sqrt(x**2 + y**2)
            return d

        for n in range(int(floor(N/2))):
            rng = range(N)
            pixels = asarray([temp_solutions[x, y] for x in rng for y in rng if dist(x, y) == n])
            bb_curve.append(mean(pixels))

        reved_bb_curve = reversed(bb_curve)
        bb_curve = self.normalize(asarray([ i for i in reved_bb_curve])) ## need to reverse data order because of the way the loop iterates
                                                # through the data, reverse orders data from 0 to max radius

        def acc_temp_theory(R):

            #consts = (3 * G * M1 * M_trans) / (8 * pi * steph_boltz)
            # T_to_fourth = consts * (R**-3) * (1 - viscosity * (Rwd/R) ** (1/2))
            # just want the shape, so only using the stephon_boltzman constant

            T = (((steph_boltz * R**3) ** (-1)) * (1 - viscosity * (radius_primary/R) ** (1/2))) ** (1/4)
            return T

        theory_temps = self.normalize(asarray(map(acc_temp_theory, radii)))
        radii = self.normalize(radii, over = self.RL1_from_primary)

        fig, ax = plt.subplots()

        ln_T1 = map(lambda x: log(x), theory_temps)
        ln_T2 = map(lambda x: log(x), bb_curve)

        ax.plot(radii, ln_T1, 'rx', label = 'Theoretical')
        ax.plot(radii, ln_T2, 'b^', label = 'Reconstructed')
        ax.legend()
        ax.set_xlabel('[R/RL1]')
        ax.set_ylabel('Temperature [K] (Normalized)')
        plt.title(''.join(['Temperature profile of ', self.system_name]))
        plt.savefig(path)

    def save_synth_vs_obs(self, path):
        import matplotlib.pyplot as plt

        synth_fluxes = self.synthetic_fluxes
        observed_fluxes = self.observed_fluxes
        times = self.JDtimes

        fig, ax = plt.subplots()
        ax.set_ylim(0, 1.25)
        ax.plot(times, synth_fluxes, 'rx', label='Synthetic')
        ax.plot(times, observed_fluxes, 'b^', label='Observed')
        ax.set_xlabel('Time [d]')
        plt.title(''.join(['Observed and Reconstructed Light Curves: ', self.system_name]))

        if self.fluxes_not_mags:
            ax.set_ylabel('Flux (Normalized)')
        else:
            ax.set_ylabel('Magnitude')

        ax.legend()

        # chi_text = r"\begin{eqnarray}" + r"\chi^2 &=& ", str(self.chi_AIM)
        # ax.text(1, 0.9, chi_text, {'color': 'C2', 'fontsize': 18}, va="top", ha="left")

        plt.savefig(path)




class Processes:

    def __init__(self):
        pass

    @staticmethod
    def save_synth_vs_obs(obj, path):
        import matplotlib.pyplot as plt

        synth_fluxes = obj.synthetic_fluxes
        observed_fluxes = obj.observed_fluxes
        times = obj.JDtimes

        fig, ax = plt.subplots()
        ax.set_ylim(0,1.25)
        ax.plot(times, synth_fluxes, 'rx', label='Synthetic')
        ax.plot(times, observed_fluxes, 'b^', label='Observed')
        ax.set_xlabel('Time [d]')
        plt.title(''.join(['Observed and Reconstructed Light Curves: ', obj.system_name]))

        if obj.fluxes_not_mags:
            ax.set_ylabel('Flux (Normalized)')
        else:
            ax.set_ylabel('Magnitude')

        ax.legend()

        #chi_text = r"\begin{eqnarray}" + r"\chi^2 &=& ", str(obj.chi_AIM)
        #ax.text(1, 0.9, chi_text, {'color': 'C2', 'fontsize': 18}, va="top", ha="left")

        plt.savefig(path)


    @staticmethod
    def build_compute_save(binary_object):
        from os.path import join, isdir
        from os import getcwd, makedirs
        import matplotlib.pyplot as plt
        '''
        :param binary: a binary object
        :return: calls everything needed to compute a maximum entropy image and saves it in a file corresponding
                    to the system name
        '''



        def save_image(obj, path, is_temp = False):
            from numpy import asanyarray

            fig, ax = plt.subplots()
            if is_temp:
                image = obj.temperature_solutions
                title = ''.join(['Accretion image of ', obj.system_name, ' ', '(Teff)'])
            else:
                image = obj.image_solution
                title = ''.join(['Accretion image of ', obj.system_name, ' ', '(Flux)'])

            im = ax.imshow(image)
            ax.set_title(title)
            plt.savefig(path)
            plt.show()


        def save_parameters(obj, path):
            import csv
            parameters = obj.all_params
            w = csv.writer(open(path, "w"))
            for key, val in parameters.items():
                w.writerow([key, val])

        def save_light_curve_data(obj, path):
            import csv
            from numpy import zeros, savetxt
            import pandas as pd

            w = csv.writer(open(path, "w"))
            w.writerow(['JDtimes', 'observed fluxes', 'observed uncertainties'])

            for i in range(len(obj.JDtimes)):
                w.writerow([obj.JDtimes[i], obj.observed_fluxes[i],
                            obj.observed_uncertainties[i]])






            #savetxt(path, X = data, delimiter= ',', header = 'times (JD), flux, uncertainty')

        def save_solutions(obj, path, is_temp = False):
            from numpy import savetxt
            if is_temp:
                savetxt(path, obj.temperature_solutions.reshape((int(obj.N ** 2))), delimiter=',')
            else:
                savetxt(path, obj.image_solution.reshape((int(obj.N**2))), delimiter= ',')

        print "Reconstructing Accretion disk for system: ", str(binary_object.system_name)

        if binary_object.fluxes_not_mags:
            print "\nLight curve data was in units of flux. Normalizing... ",
            binary_object.normalize_fluxes()
            print "Done."

        print "\n Building orbital model...",
        binary_object.build_orbital_model()
        print "Done."

        print '\n Constructing phase points...',
        binary_object.construct_phasepoints()
        print 'Done.'

        print '\n Building Eclipse Map...',
        binary_object.build_visibility_matrix()
        print "Done."

        print '\n Solving Maximum Entropy Equation...',
        binary_object.get_image_solution()
        print 'Done.'

        binary_object.get_synthetic_fluxes()

        binary_object.flux_to_temp()




        cwd = getcwd()
        binary_dir = join(cwd, 'Outputs', 'Results', binary_object.system_name)


        if not isdir(binary_dir):
            makedirs(binary_dir)

        print "\nSaving data in directory: ", binary_dir

        image_path = join(binary_dir, 'accretion_image.png')
        params_path = join(binary_dir, 'parameters.csv')
        data_path = join(binary_dir, 'light_curve_data.csv')
        solutions_path = join(binary_dir, 'solution_set.csv')
        temp_solutions_path = join(binary_dir, 'temperature_solutions.csv')
        temp_image_path = join(binary_dir, 'temperature_image.png')
        LCplot_path = join(binary_dir, 'LC_plot.png')
        bbody_path = join(binary_dir, 'temp_profile.png')

        save_image(binary_object, image_path)
        save_image(binary_object, temp_image_path, is_temp= True)
        save_solutions(binary_object, temp_solutions_path, is_temp= True)
        save_parameters(binary_object, params_path)
        save_light_curve_data(binary_object, data_path)
        save_solutions(binary_object, solutions_path)
        binary_object.save_synth_vs_obs(LCplot_path)
        binary_object.save_bbody_curve(bbody_path)


        print "\nAnalysis complete for system: ", binary_object.system_name

    @staticmethod
    def make_parameter_set(system_name, incl, period, sma, q, ecc, distance, radius_secondary, radius_primary):

        '''
        :param system_name: (str) name of the system
        :param incl: (float) orbital inclination (in degrees)
        :param period: (float) orbital period (in days)
        :param sma: (float) semi major axis of the system (in solar radii)
        :param q: (float) mass ratio (M2/M1)
        :param ecc: (float) eccentricity
        :param distance: (float) distance to system (in parsecs)
        :param radius_secondary: (float) radius of secondary star (in solar radii)
        :param radius_primary: (float) radius of primary star (in solar radii)
        :return: a dictionary of orbital parameters to pass into a Binary object
        '''

        params = {'system_name': system_name, 'incl': incl, 'period': period, 'sma':sma, 'q': q,
                  'distance': distance, 'radius_secondary': radius_secondary, 'radius_primary':radius_primary,
                  'ecc': ecc}

        return params

    @staticmethod
    def get_params_from_file(params_path):

        '''
        :param params_path: (file directory) the path to the parameters.csv file
        :return: a dictionary of the parameter set
        '''

        from numpy import loadtxt

        param_file = loadtxt(params_path)
        param_dict = {}

        for key, val in zip(param_file[: 0], param_file[:, 1]):
            param_dict[key] = val

        return param_dict

    @staticmethod
    def get_LCdata_from_file(data_path):

        '''
        :param data_path: (file directory) path to the light curve data file
        :return: an np.array that has the times, brightness values, and uncertainties organized for the algorithm
        '''

        from numpy import loadtxt

        LC_data = loadtxt(data_path, delimiter= ',')

        return LC_data



from os import getcwd
from os.path import join

cwd = getcwd()
LC_path = join(cwd, 'data', 'KIC_201325107_flux.csv')

orbit_params = Processes.make_parameter_set('KIC_201325107(22)', 76.19118, 0.7443, 5.00352, 0.5489, 0.0,
                                            12.5, 1.643034, 1.670691)
light_curve_data = Processes.get_LCdata_from_file(LC_path)

#print light_curve_data[:, 2]

b = Binary(orbit_params, resolution= 20, chi_AIM=0.125,
           delta = 3, light_curve_data = light_curve_data, method = 'SLSQP')

Processes.build_compute_save(b)

#print b.estimate_chi_sqrd()



















'''
#incl = 76.19118

params = {'q': 0.5489, 'period': 0.7443, 'sma': 5.00352, 'ecc':0.0,
        'distance': 12.5, 'incl': 76.19118, 'rad_second': 1.643034}

print('Initializing...')
b = Binary('test_system', incl = params['incl'], period= params['period'],
           distance=params['distance'], sma= params['sma'], q = params['q'],
           ecc = params['ecc'], N = 25, radius_secondary= params['rad_second'])

print('Done...')

print('Building orbital model...')
b.build_orbital_model()

print('Building visibility matrix..')
b.build_visibility_matrix()
print('Done... ')


b.animate_visibility_matrix()
'''



























































'''
def build_and_compute_binary(incl, period, sma, q, ecc):
    import numpy as np
    b = phoebe.default_binary()
    b.set_value('q', q)
    b.set_value('incl@binary@orbit@component', incl)
    b.set_value('period@binary@orbit@component', period)
    b.set_value('sma@binary@orbit@component', sma)
    b.set_value('ecc@binary@orbit@component', ecc)
    times = np.linspace(0,10,10)
    b.add_dataset(phoebe.dataset.orb, times = times, dataset= 'orb1', component= ['primary', 'secondary'])
    b.add_compute(phoebe.compute.phoebe, compute= 'preview', irrad_method = None)
    b.run_compute()
    return b
'''





'''
test_b.build_orbital_model()

prime_orb = test_b.primary_orbit
second_orb = test_b.secondary_orbit

phi = np.linspace(0,2*np.pi)

x1, y1, z1 = prime_orb(phi)

print(test_b.RL1_from_primary)
print(x1)
'''








#binary = build_and_compute_binary(params['incl'], params['period'], params['sma'], params['q'], params['ecc'])
#binary.plot(show = True)

'''
b = phoebe.default_binary()
b.add_dataset('orb')
b.set_value_all('times', np.linspace(0,3,201))
b.run_compute()

times = b['times@primary@orb01@orb@model'].get_value()
us = b['us@primary@orb01@orb@model'].get_value()
vs = b['vs@primary@orb01@orb@model'].get_value()
ws = b['ws@primary@orb01@orb@model'].get_value()
'''

'''
fig = plt.figure()
ax = Axes3D(fig)


ax.plot(xs = us, ys = vs, zs = ws)


from os import getcwd
from os.path import join
cwd = getcwd()
fname = join(cwd, 'TEST.png')

plt.savefig(fname)
'''





















