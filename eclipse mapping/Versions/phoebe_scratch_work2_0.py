


import sys





class Binary:


    def __init__(self, orbital_params,
                 light_curve_data, radius_secondary_given = True, fluxes_not_mags = True):
        import numpy as np

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
        # building phoebe model inside binary class


        b = phoebe.default_binary()
        b.set_value('q', orbital_params['q'])
        b.set_value('incl@binary@orbit@component', orbital_params['incl'])
        b.set_value('period@binary@orbit@component', orbital_params['period'])
        b.set_value('sma@binary@orbit@component', orbital_params['sma'])
        b.set_value('ecc@binary@orbit@component', orbital_params['ecc'])
        b.set_value('distance@system', orbital_params['distance'])

        if 'mass_primary' not in orbital_params.keys():
            self.mass_primary = b.get_value('mass@primary@star@component')
        else:
            self.mass_primary = orbital_params['mass_primary']

        if 'mass_secondary' not in orbital_params.keys():
            self.mass_secondary = b.get_value('mass@secondary@star@component')
        else:
            self.mass_secondary = orbital_params['mass_secondary']

        self.sma_primary = b.get_value('sma@primary@star@component')
        self.sma_secondary = b.get_value('sma@secondary@star@component')

        # building a dictionary of all parameters passed in the object, will use later for logging...
        self.all_params = orbital_params
        self.all_params['resolution'] = orbital_params['resolution']
        self.all_params['chi_AIM'] = orbital_params['chi_AIM']
        self.all_params['delta'] = orbital_params['delta']
        self.all_params['radius_secondary_given'] = radius_secondary_given
        self.all_params['method'] = orbital_params['method']
        self.all_params['mass_primary'] = self.mass_primary
        self.all_params['mass_secondary'] = self.mass_secondary

        if not self.all_params['rho']:
            print 'Rho value not found, setting rho to 1, ignore if constrained solution is being used.'
        else:
            self.rho = self.all_params['rho']

        if 'phi0' not in self.all_params.keys():

            print 'phi0 not in orbital_params, setting phi0 equal to 0'
            self.phi0 = 0

        else:
            self.phi0 = self.all_params['phi0']

        # b.run_compute()
        # defining values.. so many values
        self.system_name = orbital_params['system_name']
        self.temperature_primary = orbital_params['temperature_primary']
        self.temperature_secondary = orbital_params['temperature_secondary']
        self.radius_primary = orbital_params['radius_primary']
        self.fluxes_not_mags = fluxes_not_mags
        self.method = orbital_params['method']
        self.delta = orbital_params['delta']
        self.chi_AIM = orbital_params['chi_AIM']
        self.distance = orbital_params['distance']
        self.period = orbital_params['period']
        self.binary_sma = orbital_params['sma']  # sma already in sol rad
        self.incl = np.deg2rad((90 - orbital_params
        ['incl']))  # orbital inclination, 90 - incl becase 90 degrees is considered perpendicular
        self.ecc = orbital_params['ecc']
        self.q = orbital_params['q']
        self.N = orbital_params['resolution']  # pixel grid dimension
        self.accretion_grid = np.ndarray(shape = (self.N, self.N))
        self.light_curve_data = light_curve_data  # should be an array of shape (3, number_of_datapoints)
        # axis [0] is times, [1] is fluxes [2] is uncertainties
        """MUST CALL Processes.chop_LC BEFORE use."""
        self.epochs = None  # containts a list of arrays containing light curve data epochs
        self.JDtimes = light_curve_data[:, 0]
        self.observed_fluxes = light_curve_data[:, 1]
        self.observed_uncertainties = light_curve_data[:, 2]
        self.wjks = None

        '''CALL self.construct_phasepoints() BEFORE use'''
        self.constructed_phasepoints = None

        '''CALL self.get_synthetic_fluxes BEFORE use'''
        self.synthetic_fluxes = None


        self.delta = self.all_params['delta'] ###going to make this a variable parameter later

        # creating space to compute orbital models

        # self.primary_orbit, and self.secondary orbit return a list of [x, y, z] coordinates
        # easier for handling when building the eclipse map

        '''***MUST CALL self.build_orbital_model BEFORE using these methods***'''
        self.primary_orbit = None
        self.secondary_orbit = None

        # these functions return the individual points on each axis, easier for parametric graphing
        self.x_primary = None
        self.y_primary = None
        self.z_primary = None

        self.x_secondary = None
        self.y_secondary = None
        self.z_secondary = None

        '''***MUST CALL self.build_fractional_visibility_matrix BEFORE using these attributes***'''
        self.image_entropy = None
        self.image_solution = None
        self.secondary_contributions = []

        # gotta solve for RL1 soooooo.... use fsolve fromm scipy? (not working)  ---> play the system, make the equation
        # an absolute value function and minimize for 0 ---- ayyyyyooooo

        # python wont let me call self.attribute in a method definition inside __init__ so i'm renaming values


        def to_kg(sol_mass):
            # converts_solar mass to kg
            # 1.9889200011445836e+30

            return sol_mass * 1.9889200011445836 * 10**30


        M1, M2 = map(to_kg, [self.mass_primary, self.mass_secondary])
        sma = self.binary_sma * 6.957 * (10**8)  # converting from sol rad to meters


        def get_RL1_func(M1, M2, sma):
            # defining args of the 1st lagrange point function

            def RL1_func(r):
                # this equation gives the roots of the RL1 function measured from the 2nd star. it is set equal to zero
                # takning the absolute value and minimizing for 0 and taking the r value that gives me 0 because fsolve i
                # is being dumb

                return abs((M2 / (r ** 2)) + M1 / (sma ** 2) - ((r * (M1 + M2)) / (sma ** 3)) -
                           (M1 / ((sma - r) ** 2)))  # == 0

            return RL1_func

        # getting estimate value to start fsolve and building lagrange point function
        guesstimate = np.asarray([0])
        RL1_func = get_RL1_func(M1, M2, sma)

        # solving for RL1..

        RL1_from_secondary = minimize(RL1_func, x0=guesstimate, method='Nelder-Mead').x
        RL1_from_primary = sma - RL1_from_secondary

        self.RL1_from_primary = RL1_from_primary[0] / (695.51 * 10 ** 6)  # meters to solar radii
        self.RL1_from_secondary = RL1_from_secondary[0] / (695.51 * 10 ** 6)
        # should I multiply this by 2???
        self.pixel_to_RL1 = 2 * self.RL1_from_primary / self.N  # conversion factor, converts pixel units to RL1 units
        # all distances to convert should be in solar radii!

        if radius_secondary_given:
            self.radius_secondary = orbital_params['radius_secondary']
        else:
            # self.radius_secondary = None
            print 'Roche Lobe Equation estimate not yet implemented, please add a secondary radius to your orbital ' \
                  'parameter set (in solar radii), and make sure secondary_radius_given = True'
            sys.exit()

        # making room for visibility matrix
        self.visibility_matrix = None

        # making room for temperature solutions
        self.temperature_solutions = None


    def get_secondary_contributionv2(self, y_secondary, pixel_len = 50):
        '''
        :param y_secondary: y coordinate of secondary star
        :param z_secondary:
        :param pixel_len:
        :return: lays a grid on top of a circle with the same radius of the secondary star.
                the method loops through each pixel coordinate to check if it is eclipsing the accretion grid.
                the amount of flux to contribute to the lightcurve will be based on the ratio of eclipsing pixels to
                uneclipsing pixels
        '''
        from numpy import zeros, mean
        from pandas import Series

        secondary_grid = zeros(shape = (pixel_len))
        secnd_rad_RL1 = self.radius_secondary/self.RL1_from_primary
        RL1_to_px = pixel_len/secnd_rad_RL1

        for i in range(pixel_len):
            ypx = i
            ypx_space = abs(ypx/RL1_to_px + y_secondary)

            if ypx_space <= 1:
                secondary_grid[i] = 1

        flux_ratio = mean(secondary_grid)

        return float(flux_ratio)







    def randomized_array(self):
        from numpy.random import rand
        N = self.N
        rand_array = rand(N, N)
        return rand_array

    def to_RL1(self, a_distance):
        # converts a distance in solar-radius to units of RL1 from primary star

        return a_distance / self.RL1_from_primary

    def build_orbital_model(self):
        # creates two position functions for the primary and secondary star orbits
        # at some phase phi, rotates about the y-axis at some inclination 'incl'
        # and assigns them to self.primary_orbit, and self.secondary orbit

        # I would like this to be in __init__ but it isnt cooperating ....

        # also creates indivudual axis funtctions and assigns them to self.x1 ,x2, y1, y2, z1, z2

        from numpy import sin, cos, deg2rad

        a1 = self.sma_primary  # semi-major axis of primary orbit
        a2 = self.sma_secondary  # semi-major axis of secondary orbit
        incl = self.incl  # orbital inclination
        RL1_from_primary = self.RL1_from_primary

        phi0 = 0  # usually, self.phi0, but I think I should only account
        # for phi0 in self.construct_phasepoints()
        phi0 = deg2rad(phi0)

        def primary_orbit(phi):
            # returns the x-y-z coordinates as a list
            # units in RL1

            x = -a1 * cos(phi + phi0) * cos(incl) / RL1_from_primary
            y = -a1 * sin(phi + phi0) / RL1_from_primary
            z = -a1 * cos(phi + phi0) * sin(incl) / RL1_from_primary

            return [x, y, z]

        def secondary_orbit(phi):
            x = a2 * cos(phi + phi0) * cos(incl) / RL1_from_primary
            y = a2 * sin(phi + phi0) / RL1_from_primary
            z = a2 * cos(phi + phi0) * sin(incl) / RL1_from_primary

            return [x, y, z]

        self.primary_orbit = primary_orbit
        self.secondary_orbit = secondary_orbit

        # theres probably a better way to do this but oh well... I gotta get this working
        # I know I could use decorators, but id rather just copy and paste /RL1_from_primary at this point

        def x1(phi):
            return -a1 * cos(phi + phi0) * cos(incl) / RL1_from_primary

        def y1(phi):
            return -a1 * sin(phi + phi0) / RL1_from_primary

        def z1(phi):
            return -a1 * cos(phi + phi0) * sin(incl) / RL1_from_primary

        def x2(phi):
            return a2 * cos(phi + phi0) * cos(incl) / RL1_from_primary

        def y2(phi):
            return a2 * sin(phi + phi0) / RL1_from_primary

        def z2(phi):
            return a2 * cos(phi + phi0) * sin(incl) / RL1_from_primary

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

        # factor that centers the coordinates of the pixel grid around 0,0
        centered = (self.N - 1) / 2
        to_RL1 = self.to_RL1

        # centering xp and yp:
        xg = (xp - centered)  ##changed from xp - centered##
        yg = (centered - yp)  ###changed from centered - yp###

        # converting pixel to orbit space

        x_space = xg * self.pixel_to_RL1
        y_space = yg * self.pixel_to_RL1

        # now to rotate the pixels about the y axis

        x_space_abt_y = x_space * cos(self.incl)
        y_space_abt_y = y_space
        z_space_abt_y = x_space * sin(self.incl)

        # now to rotate about relative z axis:

        x_space_abt_y_z = to_RL1(x_space_abt_y * cos(phi) - y_space_abt_y * sin(phi))
        y_space_abt_y_z = to_RL1(x_space_abt_y * sin(phi) + y_space_abt_y * cos(phi))
        z_space_abt_y_z = to_RL1(z_space_abt_y)

        xwd, ywd, zwd = self.primary_orbit(phi)  # pulling coordinates for primary star at current phase
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
        import numpy as np
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
        from numpy import pi, ndarray, sqrt
        '''***MUST CALL self.build_orbital_model() BEFORE use'''

        # also need to implement roche lobe equation.... right now, assuming radius will be given

        # assuming the data structure will be an array, so I'm treating it that way for now, will change syntax when
        # I start integrating actual data.

        visibility_matrix = ndarray(shape=(self.N, self.N, len(self.constructed_phasepoints)))

        to_RL1 = self.to_RL1

        for i in range(len(self.constructed_phasepoints)):

            matrix_frame = ndarray(shape=(self.N, self.N))

            # pulling phase point from list of phase points and converting to radians
            phi = self.constructed_phasepoints[i] * 2 * pi

            # getting secondary star coordinates and putting them in units of RL1 from primary
            # **** was map(self.to_RL1, self.secondary_orbit(phi)), but moved conversion to RL1 to
            #      build orbital model
            '''units ok, in RL1'''
            x_2nd, y_2nd, z_2nd = self.secondary_orbit(phi)
            #self.secondary_contributions.append(0)
            #contrib = self.get_secondary_contributionv2(y_2nd)
            self.secondary_contributions.append(0)

            for x in range(self.N):  # range or range(len())??
                for y in range(self.N):  # range or range(len())??
                    '''units ok, in RL1'''
                    x_px, y_px, z_px = self.pixel_orb_coords(x, y, phi)

                    # getting the y-z projection of the distance between the orbit
                    y_z_distance = sqrt((y_2nd - y_px) ** 2 + (z_2nd - z_px) ** 2)

                    '''units ok, in RL1'''
                    secondary_radius_in_RL1 = to_RL1(self.radius_secondary)

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
            im = plt.imshow(frame, animated=True)

            ims.append([im])

        ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
                                        repeat_delay=1000)

        target = ['animation_', 'incl_', str(rad2deg(90 - self.incl)), '.html']
        cwd = os.getcwd()
        path = os.path.join(cwd, 'Outputs', 'eclipse_map_animations', ''.join(target))
        ani.save('movie.html')

    @staticmethod
    def dist_from_center(coords):
        from numpy import sqrt
        ##pretty sure these should be spatial coordinates, not pixel coordinates.. debug this later
        # its only called in self.get_wjks, so you should only have to add the conversion on at the end below
        x, y = coords[0], coords[1]
        dist = sqrt(x ** 2 + y ** 2)
        return dist

    def build_default_map(self, wjk):
        '''

        :param I: list of all the brightnesses. shape = (n**2)
        :param wjk: weight matrix flattened into 2d with shape = (n**2, n**2)
        :return: list of the weighted default map of the brightness distribution
        '''
        from numpy import asarray
        I = self.image_solution.reshape(int(self.N ** 2))
        D = asarray([sum(a * b for a, b in zip(I, wjk[j, :])) / sum(wjk[j, :]) for j in range(len(I))])

        return D



    def get_wjks(self):
        from tqdm import tqdm
        from numpy import ndarray, exp

        N = self.N
        delta = self.delta

        wjk_matrix = ndarray(shape=(N, N, int(N ** 2)))

        frame_number = 0

        # these i and j loops define the jth pixel for the weight function
        for i in tqdm(range(N), ''.join(["Building weight matrix for ", self.system_name])):
            for j in range(N):
                # centering grid around coordinates [0,0]
                xj = i - ((N - 1) / 2)
                yj = ((N - 1) / 2) - j

                jframe = ndarray(shape=(N, N))

                Rj = self.dist_from_center([xj, yj])

                # these k and l loops define the kth pixel for the weight function
                # because these loops are inside the i and j loops, k and l will complete a full cycle for every i and j

                for k in range(N):
                    for l in range(N):
                        xk = k - ((N - 1) / 2)
                        yk = ((N - 1) / 2) - l

                        Rk = self.dist_from_center([xk, yk])

                        wjk = exp(-(Rj - Rk) ** 2 / (2 * delta ** 2))

                        jframe[k, l] = wjk

                wjk_matrix[:, :, frame_number] = jframe

                frame_number += 1
        self.wjks = wjk_matrix
        return wjk_matrix

    def normalize(self, array, over=None):
        from numpy import asarray, amax
        len_arr = len(array)
        if not over:
            return asarray(map(lambda x: x / amax(array), array)).reshape((len_arr))
        else:
            return asarray(map(lambda x: x / over, array)).reshape((len_arr))

    def get_entropy_function(self):
        from numpy import asarray
        N = self.N
        wjks = self.get_wjks()
        wjks = wjks.reshape((int(N ** 2), int(N ** 2)))

        def entropy_function(I):
            from numpy import log, nan, inf, nan_to_num
            D = asarray([sum(a * b for a, b in zip(I, wjks[j, :])) / sum(wjks[j, :]) for j in range(len(I))])

            q = asarray([i / sum(D) for i in D])
            p = asarray([i / sum(I) for i in I])
            # actual function is - sum(....  but im minimizing the positive version to find the maximum of its opposite
            S = sum([pj * log(pj / qj) for pj, qj in zip(p, q)])

            if S == inf:
                S = nan_to_num(inf)
            if S == nan:
                S = nan_to_num(nan)

            return S

        return entropy_function

    def construct_phasepoints(self):
        '''
        :return: a list of phase points assigned to self.phase_points
        '''
        from numpy import amin, where, asarray

        # fluxes = asarray(self.observed_fluxes)
        jdtimes = asarray(self.JDtimes)

        # assuming lowest flux wil give the time of super conjunction
        # --> UPDATE --> now im going to account for any phi0, so data curve can start at any phase
        # time_of_superconj = jdtimes[where(fluxes == amin(fluxes))]
        # ephemeris = jdtimes[0]

        # taking the difference between time of super conjunction and ephemeris to find the phase angle at ephemeris.
        # roll_back = time_of_superconj - ephemeris

        # self.phi0 is in degrees, this will be the starting point of
        phase0 = self.phi0 / 360
        constructed_phases = []

        for i in range(len(jdtimes)):

            # was old_phase_point - roll_back
            new_phase_point = phase0
            constructed_phases.append(new_phase_point)

            if i + 1 >= len(jdtimes):
                break

            delta_phase = (jdtimes[i + 1] - jdtimes[i]) / self.period
            phase0 += delta_phase

        self.constructed_phasepoints = asarray(constructed_phases)

    def estimate_chi_sqrd(self):
        from numpy import mean, amax, amin
        M = len(self.observed_uncertainties)
        sig_max_sqrd = amax(self.observed_uncertainties)
        chi_sqrd_estimate = sum(map(lambda x: sig_max_sqrd / (x ** 2), self.observed_uncertainties)) / M
        return chi_sqrd_estimate

    def get_chi_actual(self):
        from numpy import cos, where

        N = self.N
        phase_points = self.constructed_phasepoints
        obs_fluxes = self.observed_fluxes
        synth = self.synthetic_fluxes
        sigmas = self.observed_uncertainties
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)
        flattened_matrix = self.visibility_matrix.reshape((N ** 2, len(phase_points)))
        chi_AIM = self.chi_AIM
        I = self.image_solution

        synth_minus_obs_sqrd = [((obs - syn) / sig) ** 2 for obs, syn, sig in zip(obs_fluxes, synth, sigmas)]

        chi = sum(synth_minus_obs_sqrd) / len(synth)

        self.all_params['chi_actual'] = float(chi / (self.N ** 2))

    def get_synthetic_fluxes2(self):
        '''
        :param I: the solution set of brightneses
        :return: finds synthetic fluxes with the given maximum entropy solution and assigns them to
                    self.synthetic fluxes
        '''

        '''***there should ALREADY be an image solution BEFORE this method is called, call self.get_image_solution FIRST***'''
        # phi might have to be a frame slice instead of a phase value.. but I can find said slice with the phase value

        from numpy import cos, asarray
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)

        phase_points = self.constructed_phasepoints
        Is = self.image_solution.reshape((int(self.N ** 2)))

        # reshaping eclipse map (visibility matrix) to 2d with dimensions N**2, length of list of light_curve_data
        flattened_matrix = self.visibility_matrix.reshape((self.N ** 2, len(phase_points)))

        fluxes = []

        for phi in range(len(phase_points)):
            # getting eclipse map for phase phi
            V_phi = flattened_matrix[:, phase_points[phi]]

            I_times_V = [I * V for I, V in zip(Is, V_phi)]
            flux = (theta_sqrd / 4 * self.N ** 2) * sum(I_times_V)

            fluxes.append(flux)

        self.synthetic_fluxes = asarray(fluxes)
        self.all_params['num_datapoints'] = len(self.synthetic_fluxes)

    def get_synthetic_fluxes(self):
        from numpy import cos, where, asarray

        phase_points = self.constructed_phasepoints
        N = self.N
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)
        flattened_matrix = self.visibility_matrix.reshape((N ** 2, len(phase_points)))
        #scnd_contribs = self.secondary_contributions
        chi_AIM = self.chi_AIM
        I = self.image_solution.reshape((N ** 2))

        def synth_flux(I, phi):
            index = where(phase_points == phi)[0][0]
            #scnd_contrib = scnd_contribs[index]
            V_phi = flattened_matrix[:, index]
            I_times_V = [i * v for i, v in zip(I, V_phi)]

            flux = (theta_sqrd / 4 * N ** 2) * sum(I_times_V) #+scnd_contrib
            return flux #Why is that happening I have no idea what is going on????? an embedded array?  why python

        synth_fluxes = asarray([synth_flux(I, phi) for phi in phase_points])

        self.synthetic_fluxes = self.normalize(synth_fluxes)
        self.all_params['num_datapoints'] = len(self.synthetic_fluxes)

    def get_uncosntrained_chi(self):
        from numpy import cos, where

        N = self.N
        phase_points = self.constructed_phasepoints
        obs_fluxes = self.observed_fluxes
        sigmas = self.observed_uncertainties
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)
        flattened_matrix = self.visibility_matrix.reshape((N ** 2, len(phase_points)))
        chi_AIM = self.chi_AIM
        #scnd_contribs = self.secondary_contributions

        def chi(I):
            M = len(phase_points)

            def synth_flux(I, phi):
                index = int(where(phase_points == phi)[0][0])
                #scnd_contrib = scnd_contribs[index]
                V_phi = flattened_matrix[:, index]
                I_times_V = [i * v for i, v in zip(I, V_phi)]

                flux = (theta_sqrd / 4 * N ** 2) * sum(I_times_V) #+ scnd_contrib
                return flux

            synth_fluxes = [synth_flux(I, phi) for phi in phase_points]

            synth_minus_observed_sqrd = [((synth - obs) / sig) ** 2 for synth, obs, sig in
                                         zip(synth_fluxes, obs_fluxes, sigmas)]

            chi_sqrd = sum(synth_minus_observed_sqrd) / M

            return float(chi_sqrd * (N**2))
        return chi


    def build_chi_constraint(self):
        from numpy import cos, where, mean, amax, inf, nan, nan_to_num

        N = self.N
        phase_points = self.constructed_phasepoints
        obs_fluxes = self.observed_fluxes
        sigmas = self.observed_uncertainties
        theta_sqrd = ((self.RL1_from_primary / self.distance) ** 2) * cos(self.incl)
        flattened_matrix = self.visibility_matrix.reshape((N ** 2, len(phase_points)))
        chi_AIM = self.chi_AIM
       # scnd_contribs = self.secondary_contributions



        def chi_constraint(I):
            M = len(phase_points)

            def synth_flux(I, phi):
                index = int(where(phase_points == phi)[0][0])
                #scnd_contrib = scnd_contribs[index]
                V_phi = flattened_matrix[:, index]
                I_times_V = [i * v for i, v in zip(I, V_phi)]


                flux = (theta_sqrd / 4 * N ** 2) * sum(I_times_V) #+ scnd_contrib
                return flux

            synth_fluxes = [synth_flux(I, phi) for phi in phase_points]

            synth_minus_observed_sqrd = [((synth - obs) / sig) ** 2 for synth, obs, sig in
                                         zip(synth_fluxes, obs_fluxes, sigmas)]

            chi_sqrd = sum(synth_minus_observed_sqrd) / M

            return float(chi_sqrd - chi_AIM * (N ** 2))

        return chi_constraint


    def get_jacobian(self, f):
        '''
        :param f: a function to find the jacobian of
        :return: the jacobian function of f
        '''
        '''
        from numpy import zeros, asarray
        def J(x, dx = 10 ** -8):
         '''   '''
            :param x: an array of values
            :param dx: dx increment
            :return: the jacobian matrix as a flattened array
            '''
        '''
            n = len(x)
            func = f(x)
            jac = zeros((n, n))
            for j in range(n):  # through columns to allow for vector addition
                Dxj = (abs(x[j]) * dx if x[j] != 0 else dx)
                x_plus = [(xi if k != j else xi + Dxj) for k, xi in enumerate(x)]
                jac[:, j] = (f(x_plus) - func) / Dxj
            return asarray(jac).reshape((n**2))
        '''
        from numdifftools import Jacobian
        def J(x):
            return Jacobian(f)(x).ravel()
        return J




    def get_image_solution(self):
        from scipy.optimize import minimize, NonlinearConstraint
        from numpy import amax, asarray, inf
        import time
        from numdifftools import Jacobian


        N = self.N
        S = self.get_entropy_function()
        guess = self.gkern().reshape((int(N) ** 2))
        chi_sqrd = self.build_chi_constraint()

        # , {'type':'ineq', 'fun': ineq}
        method = self.method
        cons = [{'type': 'eq', 'fun': chi_sqrd}]

        lb = (self.chi_AIM - 0.01) * (
                self.N ** 2)  # defining bounds, moving chi AIM out of self.build chi constraint in this version
        ub = (self.chi_AIM + 0.01) * (self.N ** 2)
        # changed cons function to the line below in this version
        # nonlincons = NonlinearConstraint(chi_sqrd, lb, ub)
        print"Solving entropy equation with resolution: ", self.N, ' by ', self.N
        start = time.time()
        # method usually = method, but im changing to 'trust-constr' for testing, since i think it might use better coonstraints

        if self.method == 'Newton-CG':
            jac = self.get_jacobian(S)
            res = minimize(S, x0=guess, method=self.method, constraints=cons, jac = jac)
        else:
            res = minimize(S, x0=guess, method= self.method, constraints = cons)
        stop = time.time()

        solve_time = (stop - start)/(60**2)
        time_text = ''.join([str(solve_time), ' hrs'])
        self.all_params['solve time'] = time_text

        image_entropy = S(res.x)
        image_solution = res.x
        sol_max = amax(image_solution)
        normalized_solution = map(lambda x: x / sol_max, image_solution)

        self.image_entropy = image_entropy
        self.image_solution = asarray(normalized_solution).reshape((N, N))

    def normalize_fluxes(self):
        '''
        :return: a normalized data set of fluxes and uncertainties...only use on fluxes, not magnitudes
        '''

        from numpy import amax, asarray

        fluxes = self.observed_fluxes
        sigmas = self.observed_uncertainties

        flux_max = amax(fluxes)

        def norm_sig(sig, flux):
            return sig / flux

        norm_fluxes = asarray(map(lambda flux: flux / flux_max, fluxes))
        norm_sigmas = asarray(map(lambda sig: sig / flux_max, sigmas))
        # norm_sigmas = asarray([norm_sig(sig, flux) for sig, flux in zip(sigmas, fluxes)])

        self.observed_fluxes = norm_fluxes
        self.observed_uncertainties = norm_sigmas

    def flux_to_temp(self):
        from numpy import asarray, amax

        image_solutions = self.image_solution.reshape(self.N ** 2)
        temp_primary = self.temperature_primary

        # temp1 = self.primary_temperature # need to make this an attribute

        def temp_func(I):
            steph_boltz = 5.67 * 10 ** -8  # removing setphan boltzman constant for now, was (I / steph_boltz) ** (1/4)
            flux_primary = steph_boltz * temp_primary ** 4
            T = ((I * flux_primary) / steph_boltz) ** (0.25)
            return T / temp_primary

        temp_solutions = map(temp_func, image_solutions)

        # normalizing
        temp_solutions = asarray(map(lambda x: x / amax(temp_solutions), temp_solutions))
        # temp_solutions = asarray(map(lambda x: x * temp1, temp_solutions))

        self.temperature_solutions = temp_solutions.reshape((self.N, self.N))

    def save_bbody_curve(self, path):
        '''returns a list of radii and average azimuthal temperature to plot'''
        from numpy import mean, asarray, floor, log
        import matplotlib.pyplot as plt

        bb_curve = []
        N = self.N
        conv = self.pixel_to_RL1  # a conversion factor
        temp_solutions = self.temperature_solutions
        radius_primary = self.radius_primary
        radii = asarray([i for i in reversed(map(lambda x: x * conv, range(int(floor(N / 2)))))])
        steph_boltz = 5.67 * 10 ** -8
        viscosity = 0.25  # will make an attribute later
        temp_primary = self.temperature_primary

        def dist(x, y):
            from numpy import sqrt
            # takes pixel coordinates, centers them around 0,0 and calculates the distance from the center
            centered = (N - 1) / 2
            x = x - centered
            y = centered - y
            d = sqrt(x ** 2 + y ** 2)
            return d

        for n in range(int(floor(N / 2))):
            rng = range(N)
            pixel_ring = asarray([temp_solutions[x, y] for x in rng for y in rng if dist(x, y) == n])
            bb_curve.append(mean(pixel_ring))

        reved_bb_curve = reversed(bb_curve)
        bb_curve = asarray(
            [i for i in reved_bb_curve])  ## need to reverse data order because of the way the loop iterates

        # through the data, reverse orders data from 0 to max radius

        def acc_temp_theory(R):
            # consts = (3 * G * M1 * M_trans) / (8 * pi * steph_boltz)
            # T_to_fourth = consts * (R**-3) * (1 - viscosity * (Rwd/R) ** (1/2))
            # just want the shape, so only using the stephon_boltzman constant

            # (((steph_boltz * R**3) ** (-1)), removing stephan boltzman for now
            # T = (((R_in_m**3) ** (-1)) * (1 - viscosity * ((radius_primary*(695.51 * 10**6))/R) ** (0.5))) ** (0.25)

            term1 = (radius_primary / R) ** 0.75
            numer = 1 - viscosity * ((radius_primary / R) ** 0.5)
            denom = 1 - viscosity

            T_ratio = term1 * ((numer / denom) ** 0.25)

            return T_ratio * temp_primary

        theory_temps = asarray(map(acc_temp_theory, radii))
        # theory_temps = self.normalize(asarray(map(acc_temp_theory, radii)))
        radii = self.normalize(radii, over=self.RL1_from_primary)

        fig, ax = plt.subplots()

        # ln_T1 = map(lambda x: log(x), theory_temps)
        # ln_T2 = map(lambda x: log(x), bb_curve)

        T1 = self.normalize(theory_temps)
        T2 = self.normalize(bb_curve)

        ax.plot(radii, T1, 'rx', label='Theoretical')
        # ax.plot(radii, T2, 'b^', label = 'Reconstructed')
        # ax.legend()
        ax.set_xlabel('[R/RL1]')
        ax.set_ylabel('T/Twd')
        plt.title(''.join(['Temperature profile of ', self.system_name]))
        plt.savefig(path)

    def get_bbody_curve(self, path):
        from numpy import asarray, floor, log, pi, amax, amin, linspace, sort, arange
        import matplotlib.pyplot as plt
        N = self.N
        conv = self.pixel_to_RL1  # a conversion factor
        temps = self.temperature_solutions.reshape(N**2)
        radius_primary = self.radius_primary
        temp_primary = self.temperature_primary
        #radii = asarray([i for i in reversed(map(lambda x: x * conv, range(int(floor(N / 2)))))])
        steph_boltz = 5.67 * 10 ** -8
        viscosity = 0.25  # will make an attribute later

        def dist(x, y):
            from numpy import sqrt
            # takes pixel coordinates, centers them around 0,0 and calculates the distance from the center
            centered = (N - 1) / 2
            x = x - centered
            y = centered - y
            d = sqrt(x ** 2 + y ** 2)
            return d



        radii = []
        #getting distance for each temperature solution
        for j in range(N):
            for i in range(N):

                radius = dist(i, j) * conv
                radii.append(radius)

        def acc_temp_theory(R):
            # consts = (3 * G * M1 * M_trans) / (8 * pi * steph_boltz)
            # T_to_fourth = consts * (R**-3) * (1 - viscosity * (Rwd/R) ** (1/2))
            # just want the shape, so only using the stephon_boltzman constant

            # (((steph_boltz * R**3) ** (-1)), removing stephan boltzman for now
            # T = (((R_in_m**3) ** (-1)) * (1 - viscosity * ((radius_primary*(695.51 * 10**6))/R) ** (0.5))) ** (0.25)

            term1 = (radius_primary / R) ** 0.75
            numer = 1 - viscosity * ((radius_primary / R) ** 0.5)
            denom = 1 - viscosity

            T_ratio = term1 * ((numer / denom) ** 0.25)

            temp =  T_ratio * temp_primary
            return temp

        def temp_theory2(R):
            Tmax = amax(temps) * temp_primary
            Tmin = amin(temps) * temp_primary
            return (abs(R - 1)**(0.75))*(Tmax - Tmin) + Tmin

        '''
        theory_temps = self.normalize(asarray(map(temp_theory2, radii)))
        theory_temps = asarray(map(lambda x: abs(x - amax(theory_temps)) * temp_primary, theory_temps))
        theory_temps = self.normalize(theory_temps)
        temps = asarray(map(lambda x: x * self.temperature_primary, temps))
        radii = self.normalize(radii, over = self.RL1_from_primary)
        
        
        Tmax = amax(temps)
        T0 = amin(temps) + 110

        theory_temps = asarray([a * (Tmax - T0) + T0 for a  in theory_temps])
        '''


        #radii_theory = sort(asarray(list(set(radii))))
        #radii_theory = asarray(map(lambda x: -1*x, radii_theory))
        #theory_temps = asarray(map(temp_theory2, radii_theory))
        #radii_theory = asarray(map(lambda x: -1 * x, radii_theory))

        #radii_theory = arange(0,1,0.1)
        #theory_temps = temp_theory2(radii_theory)



        plt.figure()
        plt.scatter(radii, temps, s=pi * 3, c=(0, 0, 0), label = 'Calculated')
        #plt.scatter(radii_theory, theory_temps, s = pi*3, c = (1, 0, 0), alpha = 0.7, marker = 'x',
        #           label = 'Theoretical')

        #plt.title('Temperature vs. Radius')
        plt.ylabel('Temperature [K]')
        plt.xlabel('Radius [R/RL1]')
        plt.xlim(right = 1)
        plt.legend()
        plt.savefig(path)




    def get_bbody_curve2(self, path):
        from numpy import asarray, floor, log, pi, amax, amin, linspace, nan, isnan
        import matplotlib.pyplot as plt
        N = self.N
        conv = self.pixel_to_RL1  # a conversion factor
        temps = self.temperature_solutions
        radius_primary = self.radius_primary
        temp_primary = self.temperature_primary
        # radii = asarray([i for i in reversed(map(lambda x: x * conv, range(int(floor(N / 2)))))])
        steph_boltz = 5.67 * 10 ** -8
        viscosity = 0.25  # will make an attribute later

        def dist(x, y):
            from numpy import sqrt
            # takes pixel coordinates, centers them around 0,0 and calculates the distance from the center
            centered = (N - 1) / 2
            x = x - centered
            y = centered - y
            d = sqrt(x ** 2 + y ** 2)
            return d


        radii = []
        #getting distance for each temperature solution
        for j in range(N):
            for i in range(N):
                radius = dist(i, j) * conv
                if radius > self.RL1_from_primary:
                    temps[i, j] = nan
                    continue
                radii.append(radius)

        temps = temps.reshape(N**2)
        temps = temps[~isnan(temps)]
        temps_calc = [temp * temp_primary for temp in temps]
        radii_calc = self.normalize(radii, over = self.RL1_from_primary)


        def temp_theory2(R):
            Tmax = amax(temps) * temp_primary
            Tmin = amin(temps) * temp_primary
            const = amax(radii)/self.RL1_from_primary
            return (abs(R - 1)**(0.75))*(Tmax - Tmin) + Tmin

        radii_theory = linspace(0,1)
        temps_theory = temp_theory2(radii_theory)


        plt.figure()
        plt.scatter(radii_calc, temps_calc, s=pi * 3, c=(0, 0, 0), label = 'Calculated')
        plt.plot(radii_theory, temps_theory, 'r--', label = 'Theoretical')
        #plt.scatter(radii_theory, temps_theory, s = pi*3, c = (1, 0, 0), alpha = 0.7, marker = 'x',
        #            label = 'Theoretical')
        plt.title('Temperature vs. Radius')
        plt.ylabel('Temperature [K]')
        plt.xlabel('Radius [R/RL1]')
        plt.ylim(bottom = amin(temps) * temp_primary, top = amax(temps)* temp_primary)
        #plt.xlim(right = 1)
        plt.legend()
        plt.savefig(path)


    def save_synth_vs_obs(self, path):
        import matplotlib.pyplot as plt

        synth_fluxes = self.synthetic_fluxes
        observed_fluxes = self.observed_fluxes
        times = self.JDtimes

        fig, ax = plt.subplots()
        ax.set_ylim(0, 1.25)
        ax.plot(times, observed_fluxes, 'b^', label='Observed')
        ax.plot(times, synth_fluxes, 'rx', label='Synthetic')
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


    def get_unconstrained(self, rho):

        N = self.N
        chi_aim = self.chi_AIM
        chi_func = self.get_uncosntrained_chi()
        S = self.get_entropy_function()

        def unconstrained(params):
            '''
            :param params: the brightness values as a 1d array,
                           with rho and chi parameter taced on to the end.
            :return: returns the function: F = S - C**2/rho to be maximized
                        where C = (chi**2 - chi_aim**2)/chi_aim**2
            '''

            I = params[:int(N**2 - 1)]
            chi_param = params[int(N**2)]

            def C(chi):
                return (chi - chi_aim * (N**2))/(chi_aim ** N)

            F = (S(I) - ((C(chi_func(I))**2))/rho)

            return F
        return unconstrained


    def get_unconstrained_solution(self):
        from scipy.optimize import minimize
        from numdifftools import Jacobian, Hessian
        from numpy import append, amin, asarray
        import time

        print 'Solving unconstrained entropy equation with resolution '\
            , self.N, ' by ', self.N

        F = self.get_unconstrained(rho = self.rho)
        params_guess = append(self.gkern().reshape(self.N**2), self.chi_AIM)
        start = time.time()
        if self.method == 'Newton-CG':
            jac = self.get_jacobian(F)
            hess = Hessian(F)
            solutions = minimize(F, x0=params_guess, method=self.method, jac = jac).x
        else:
            solutions = minimize(F, x0 = params_guess, method = self.method).x
        stop = time.time()

        solve_time = (stop - start) / (60 ** 2)
        time_text = ''.join([str(solve_time), ' hrs'])
        self.all_params['solve time'] = time_text

        chi_actual = solutions[-1]
        image_solutions = solutions[:self.N**2]
        image_solutions = asarray(map(lambda x: x - amin(image_solutions), image_solutions))
        image_solutions = self.normalize(image_solutions)
        self.image_entropy = F(solutions)

        self.image_solution = image_solutions.reshape((self.N, self.N))
        self.all_params['chi_actual'] = chi_actual

    #@staticmethod
    #def dsdk(k):









class Processes:

    def __init__(self):
        pass

    @staticmethod
    def get_params_from_file(params_path):
        '''
        :param params_path: (file directory) the path to the parameters.csv file
        :return: a dictionary of the parameter set
        '''

        from numpy import loadtxt
        import pandas as pd

        param_file = pd.read_csv(params_path, sep=',', skiprows=0).values

        param_dict = {}

        for i in range(len(param_file[:, 0])):
            key, val = param_file[i, 0], param_file[i, 1]
            param_dict[key] = val

        to_float = ['ecc', 'chi_AIM', 'period', 'distance', 'incl',
                    'radius_primary', 'q', 'radius_secondary', 'sma', 'delta',
                    'mass_primary', 'mass_secondary', 'temperature_primary', 'temperature_secondary', 'phi0', 'rho']

        to_int = ['resolution']

        to_bool = ['radius_secondary_given']

        # setting values to correct data types
        for key in param_dict.keys():
            if key in to_float:
                param_dict[key] = float(param_dict[key])
            elif key in to_int:
                param_dict[key] = int(param_dict[key])
            elif key in to_bool:
                param_dict[key] = bool(param_dict[key])

        return param_dict

    @staticmethod
    def get_LC_data_from_file(LC_data_path):
        from pandas import read_csv
        from numpy import asarray
        data = read_csv(LC_data_path, skiprows=0).values

        return asarray(data)

    @staticmethod
    def get_closest(val, array):
        from numpy import argmin
        index = argmin(map(lambda x: abs(x - val), array))
        return index

    @staticmethod
    def save_synth_vs_obs(obj, path):
        import matplotlib.pyplot as plt

        synth_fluxes = obj.synthetic_fluxes
        observed_fluxes = obj.observed_fluxes
        times = obj.JDtimes

        fig, ax = plt.subplots()
        ax.set_ylim(0, 1.25)
        ax.plot(times, synth_fluxes, 'rx', label='Synthetic')
        ax.plot(times, observed_fluxes, 'b^', label='Observed')
        ax.set_xlabel('Time [d]')
        plt.title(''.join(['Observed and Reconstructed Light Curves: ', obj.system_name]))

        if obj.fluxes_not_mags:
            ax.set_ylabel('Flux (Normalized)')
        else:
            ax.set_ylabel('Magnitude')

        ax.legend()

        # chi_text = r"\begin{eqnarray}" + r"\chi^2 &=& ", str(obj.chi_AIM)
        # ax.text(1, 0.9, chi_text, {'color': 'C2', 'fontsize': 18}, va="top", ha="left")

        plt.savefig(path)

    @staticmethod
    def plot_save_LC(fluxes, times, path, title='Flux vs. Time', xlabel='Time', ylabel='Flux'):
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(times, fluxes)
        plt.title = title
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.savefig(path)

    @staticmethod
    def save_image(obj, path, is_temp=False):
        from numpy import asanyarray
        import matplotlib.pyplot as plt
        from matplotlib import cm

        fig, ax = plt.subplots()
        if is_temp:
            image = obj.temperature_solutions
            title = ''.join(['Accretion image of ', obj.system_name, ' ', '(Teff)'])
        else:
            image = obj.image_solution
            title = ''.join(['Accretion image of ', obj.system_name, ' ', '(Flux)'])
        #gist_heat is a good one
        cmap = plt.get_cmap('magma')

        im = ax.imshow(image, cmap= cmap)
        ax.set_title(title)
        plt.savefig(path)
        plt.show()

    @staticmethod
    def save_parameters(obj, path):
        import csv
        parameters = obj.all_params
        w = csv.writer(open(path, "w"))
        w.writerow(['', ''])
        for key, val in parameters.items():
            w.writerow([key, val])

    @staticmethod
    def save_light_curve_data(obj, path, synthetic=False):
        import csv
        from numpy import zeros, savetxt
        import pandas as pd

        w = csv.writer(open(path, "w"))

        if synthetic:
            w.writerow(['JDtimes', 'synthetic fluxes'])
            for i in range(len(obj.JDtimes)):
                w.writerow([obj.JDtimes[i], obj.synthetic_fluxes[i]])
        else:
            w.writerow(['JDtimes', 'observed fluxes', 'observed uncertainties'])
            for i in range(len(obj.JDtimes)):
                w.writerow([obj.JDtimes[i], obj.observed_fluxes[i], obj.observed_uncertainties[i]])

        # savetxt(path, X = data, delimiter= ',', header = 'times (JD), flux, uncertainty')

    @staticmethod
    def save_LC_data_from_array(data, path, syntheric = False):
        import csv

        w = csv.writer(open(path, "w"))
        w.writerow(['JDtimes', 'observed fluxes', 'observed uncertainties'])
        for i in range(len(data[:, 0])):
            w.writerow([data[i, 0], data[i, 1], data[i, 2]])

    @staticmethod
    def save_solutions(obj, path, is_temp=False):
        from numpy import savetxt
        if is_temp:
            savetxt(path, obj.temperature_solutions.reshape((int(obj.N ** 2))), delimiter=',')
        else:
            savetxt(path, obj.image_solution.reshape((int(obj.N ** 2))), delimiter=',')

    @staticmethod
    def save_synthetic_flux_data(binary_object, path):
        from numpy import savetxt
        savetxt(path, binary_object.synthetic_fluxes, delimiter=',')

    @classmethod
    def build_compute_save(self, binary_object):
        from os.path import join, isdir
        from os import getcwd, makedirs
        import matplotlib.pyplot as plt
        '''
        :param binary: a binary object
        :return: calls everything needed to compute a maximum entropy image and saves it in a file corresponding
                    to the system name

                    saves parameters, light curve, image solution, temperature solution, temperature profile,
                    along with all corresponding csv files of data for later use
        '''

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

        binary_object.get_chi_actual()

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
        LC_all_plot_path = join(binary_dir, 'LC_all_plot.png')
        bbody_path = join(binary_dir, 'temp_profile.png')
        LC_synth_data_path = join(binary_dir, 'LC_synthetic_data.csv')
        LC_synth_plot_path = join(binary_dir, 'LC_synth_plot.png')
        LC_obs_plot_path = join(binary_dir, 'LC_obs_plot.png')

        self.save_image(binary_object, image_path)

        self.save_image(binary_object, temp_image_path, is_temp=True)
        self.save_solutions(binary_object, temp_solutions_path, is_temp=True)

        self.save_parameters(binary_object, params_path)
        self.save_light_curve_data(binary_object, data_path)
        self.save_light_curve_data(binary_object, LC_synth_data_path, synthetic=True)
        self.save_solutions(binary_object, solutions_path)
        binary_object.save_synth_vs_obs(LC_all_plot_path)
        synth_title = 'Synthetic Flux vs. Time [JD]'
        self.plot_save_LC(binary_object.synthetic_fluxes, binary_object.JDtimes, LC_synth_plot_path, title=synth_title,
                          xlabel='Time [JD]', ylabel='Flux (Normalized)')

        obs_title = 'Observed Flux vs. Time [JD]'
        self.plot_save_LC(binary_object.synthetic_fluxes, binary_object.JDtimes, LC_synth_plot_path, title=obs_title,
                          xlabel='Time [JD]', ylabel='Flux (Normalized)')
        binary_object.get_bbody_curve(bbody_path)
        # binary_object.save_bbody_curve(bbody_path)

        print "\nAnalysis complete for system: ", binary_object.system_name



    @classmethod
    def build_compute_save_unconstrained(self, binary_object):
        from os.path import join, isdir
        from os import getcwd, makedirs
        import matplotlib.pyplot as plt
        '''
        :param binary: a binary object
        :return: calls everything needed to compute a maximum entropy image and saves it in a file corresponding
                    to the system name

                    saves parameters, light curve, image solution, temperature solution, temperature profile,
                    along with all corresponding csv files of data for later use
        '''

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
        binary_object.get_unconstrained_solution()
        print 'Done.'

        binary_object.get_synthetic_fluxes()

        binary_object.flux_to_temp()

        # binary_object.get_chi_actual()

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
        LC_all_plot_path = join(binary_dir, 'LC_all_plot.png')
        bbody_path = join(binary_dir, 'temp_profile.png')
        LC_synth_data_path = join(binary_dir, 'LC_synthetic_data.csv')
        LC_synth_plot_path = join(binary_dir, 'LC_synth_plot.png')
        LC_obs_plot_path = join(binary_dir, 'LC_obs_plot.png')

        self.save_image(binary_object, image_path)

        self.save_image(binary_object, temp_image_path, is_temp=True)
        self.save_solutions(binary_object, temp_solutions_path, is_temp=True)

        self.save_parameters(binary_object, params_path)
        self.save_light_curve_data(binary_object, data_path)
        self.save_light_curve_data(binary_object, LC_synth_data_path, synthetic=True)
        self.save_solutions(binary_object, solutions_path)
        binary_object.save_synth_vs_obs(LC_all_plot_path)
        synth_title = 'Synthetic Flux vs. Time [JD]'
        self.plot_save_LC(binary_object.synthetic_fluxes, binary_object.JDtimes, LC_synth_plot_path,
                          title=synth_title,
                          xlabel='Time [JD]', ylabel='Flux (Normalized)')

        obs_title = 'Observed Flux vs. Time [JD]'
        self.plot_save_LC(binary_object.synthetic_fluxes, binary_object.JDtimes, LC_synth_plot_path,
                          title=obs_title,
                          xlabel='Time [JD]', ylabel='Flux (Normalized)')
        binary_object.get_bbody_curve(bbody_path)
        # binary_object.save_bbody_curve(bbody_path)

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

        params = {'system_name': system_name, 'incl': incl, 'period': period, 'sma': sma, 'q': q,
                  'distance': distance, 'radius_secondary': radius_secondary, 'radius_primary': radius_primary,
                  'ecc': ecc}

        return params

    @staticmethod
    def get_LCdata_from_file(data_path):

        '''
        :param data_path: (file directory) path to the light curve data file
        :return: an np.array that has the times, brightness values, and uncertainties organized for the algorithm
        '''

        from numpy import loadtxt

        LC_data = loadtxt(data_path, delimiter=',')

        return LC_data

    @classmethod
    def reload_binary(self, binary_dir):

        from os.path import join, exists
        from numpy import asarray, genfromtxt

        solutions_path = join(binary_dir, 'solution_set.csv')
        LC_path = join(binary_dir, 'light_curve_data.csv')
        params_path = join(binary_dir, 'parameters.csv')

        light_curve_data = self.get_LC_data_from_file(LC_path)
        params = self.get_params_from_file(params_path)

        # image_solutions = asarray(pd.read_csv(solutions_path))
        # image_solutions = asarray(genfromtxt(solutions_path))
        b = Binary(params, light_curve_data)
        print '\n Reloaded Binary System: ', b.system_name

        # b.image_solution = asarray(image_solutions).reshape((b.N, b.N))

        if exists(solutions_path):
            image_solution = asarray(genfromtxt(solutions_path)).reshape((b.N, b.N))
            b.image_solution = image_solution
        else:
            print("Image solutions didnt exist... using temperature solutions instead")
            temp_path = join(binary_dir, 'temperature_solutions.csv')
            image_solution = asarray(genfromtxt(temp_path)).reshape((b.N, b.N))
            image_solution = b.normalize(image_solution)

            b.image_solution = image_solution


        print 'Loading phasepoints, orbital model, visibility matrix, and synthetic fluxes...',

        b.construct_phasepoints()
        b.build_orbital_model()
        b.build_visibility_matrix()
        b.get_synthetic_fluxes()

        print 'Done'

        return b

    @classmethod
    def analyze_solved_binary(self, binary_dir):
        '''

        :param binary_dir: a directory to a star system file, containing an parameters, a lightcurve,
                            and a previously solved image solution
        :return: saves temperature solution data, image, temperature profile, and Lightcurve plot
        '''

        from os.path import join

        binary_object = self.reload_binary(binary_dir)

        image_path = join(binary_dir, 'accretion_image.png')
        params_path = join(binary_dir, 'parameters.csv')
        data_path = join(binary_dir, 'light_curve_data.csv')
        solutions_path = join(binary_dir, 'solution_set.csv')
        LCplot_path = join(binary_dir, 'LC_plot.png')
        bbody_path = join(binary_dir, 'temp_profile.png')
        temp_solutions_path = join(binary_dir, 'temperature_solutions.csv')
        temp_image_path = join(binary_dir, 'temperature_image.png')

        print"Analyzing previously solved image for system: ", binary_object.system_name
        print '\n Getting synthetic fluxes, temperature profile and light curve plot...',

        binary_object.flux_to_temp()

        self.save_image(binary_object, temp_image_path, is_temp=True)
        self.save_solutions(binary_object, temp_solutions_path, is_temp=True)
        binary_object.save_synth_vs_obs(LCplot_path)
        binary_object.get_bbody_curve2(bbody_path)
        binary_object.get_chi_actual()
        self.save_parameters(binary_object, params_path)
        print 'Done.'
        print '\n Analysis complete for system: ', binary_object.system_name
        print 'Saved in ', binary_dir

    @classmethod
    def load_binary(self, binary_dir):
        from os.path import join

        LC_path = join(binary_dir, 'light_curve_data.csv')
        params_path = join(binary_dir, 'parameters.csv')

        light_curve_data = self.get_LC_data_from_file(LC_path)
        params = self.get_params_from_file(params_path)

        b = Binary(params, light_curve_data)
        print '\n Loaded Binary System: ', b.system_name
        return b

    @staticmethod
    def chop_LC(b, stop=None, chop_at='all', time_error=0.1):
        # chops LC data into epochs, assumes time incriment of observation is consistant throughout dataset
        from scipy.signal import find_peaks
        from numpy import floor, mean, diff, asarray, amin

        epochs = []
        data = b.light_curve_data

        # subtracting smallest time to start time at 0
        data[:, 0] = asarray(map(lambda x: x - amin(data[:, 0]), data[:, 0]))
        # flipping flux data set to find peaks as primary minima
        data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

        # period time in units of index
        period = floor(b.period / data[:, 0][1] - data[:, 0][0])
        # correcttime between observations
        t_diff = data[1:0] - data[0:0]

        # peaks is a list of indicies where each index is a peak
        peaks, _ = find_peaks(data[:, 1], distance=period)
        data = data[peaks[0]:peaks[-1] + 1, :]

        # flipping flux data back over
        data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

        for i in range(len(peaks - 2)):
            if chop_at == 'all':
                if i == stop:
                    break
                elif data[0: i + 1] - data[0: i] >= t_diff * time_error:
                    pass
                else:
                    start_chop, stop_chop = peaks[i], peaks[i + 1]

                    epoch = asarray(data[start_chop, stop_chop + 1])
                    epochs.append(epoch)

            elif chop_at != 'all':
                if type(chop_at) != isinstance(chop_at, int):
                    print('Enter an integer value for "chop_at, and restart, chop_at can either be an integer or string'
                          ' "all"')
                else:
                    if chop_at != i:
                        pass
                    else:
                        start_chop, stop_chop = peaks[i], peaks[i + 1]

                        epoch = asarray(data[start_chop, stop_chop + 1])
                        epochs.append(epoch)
        b.epochs = epochs

    @staticmethod
    def chop_LC_from_file(path, period, stop=None, chop_at=None, time_error=0.1):
        # chops LC data into epochs, assumes time incriment of observation is consistant throughout dataset
        from scipy.signal import find_peaks
        from numpy import floor, mean, diff, asarray, amin, loadtxt

        epochs = []
        data = loadtxt(path, delimiter=',', skiprows=1)

        # subtracting smallest time to start time at 0
        data[:, 0] = asarray(map(lambda x: x - amin(data[:, 0]), data[:, 0]))
        # flipping flux data set to find peaks as primary minima
        data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

        # period time in units of index
        period_in_points = floor(period / (data[:, 0][1] - data[:, 0][0]))
        # correcttime between observations
        t_diff = data[1:0] - data[0:0]

        # peaks is a list of indicies where each index is a peak
        peaks, _ = find_peaks(data[:, 1], distance=period_in_points)
        # data = data[peaks[0]:peaks[-1] + 1, :]

        # flipping flux data back over
        data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

        for i in range(len(peaks) - 2):
            start_chop, stop_chop = peaks[i - 1], peaks[i + 1]

            if chop_at == None:
                if i == stop:
                    break
                elif data[i + 1, 0] - data[i, 0] >= t_diff * time_error:
                    pass
                else:

                    epoch = asarray(data[start_chop:stop_chop + 1, :])
                    epochs.append(epoch)

            elif chop_at != None:
                if type(chop_at) is not int:
                    print('Enter an integer value for "chop_at, and restart, chop_at can either be an integer or string'
                          ' "all"')
                else:
                    if chop_at != i:
                        pass
                    else:

                        epoch = asarray(data[start_chop:stop_chop + 1, :])
                        epochs.append(epoch)
        return epochs

    @staticmethod
    def save_epochs(dir, epochs):
        from pandas import DataFrame as df
        import matplotlib.pyplot as plt
        from os.path import join

        def save_data(data, path):
            import csv

            w = csv.writer(open(path, "w"))
            w.writerow(['JDtimes', 'observed fluxes', 'observed uncertainties'])

            for i in range(len(data[:, 0])):
                w.writerow([data[i, 0], data[i, 1],
                            data[i, 2]])

        for i in range(len(epochs)):
            if i == 0:
                pass
            elif i % 2 == 0:
                pass
            else:
                data = epochs[i]
                plt.figure()
                plt.plot(data[:, 0], data[:, 1], 'bx')
                title = ''.join(['KIC_201325107 ', 'Epoch at: ', str(data[0, 0]), ' days'])
                plt.title(title)

                fname = ''.join(['epoch_at_', str(data[0, 0]), '.png'])
                img_path = join(dir, fname)
                plt.savefig(img_path)

                fname = ''.join(['epoch_at_', str(data[0, 0]), '.csv'])
                data_path = join(dir, fname)
                save_data(data, data_path)
                plt.close()

    @classmethod
    def plot_peaks(self, LC_path, period, start=None, stop=None, color='bx', file_name='plotted_peaks.png'):
        from numpy import asarray, floor, loadtxt
        from scipy.signal import find_peaks
        import matplotlib.pyplot as plt
        from os.path import join

        data = loadtxt(LC_path, delimiter=',', skiprows=1)
        period_in_points = floor(period / (data[:, 0][1] - data[:, 0][0]))

        if start != None:
            start = self.get_closest(start, data[:, 0])

        if stop != None:
            stop = self.get_closest(stop, data[:, 0])

        if start == None and stop == None:
            pass

        elif start != None and stop == None:
            data = data[start:, :]

        elif start == None and stop != None:
            data = data[:stop, :]

        elif start != None and stop != None:
            data = data[start:stop, :]
        else:
            print 'Something strange happened with the start and stop indicies, plotting all instead'
            data = data[:, :]

        # flipping flux data set to find peaks as primary minima
        data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

        # period time in units of index
        period = floor(period / data[:, 0][1] - data[:, 0][0])
        # correcttime between observations
        t_diff = data[1:0] - data[0:0]

        # peaks is a list of indicies where each index is a peak
        peaks, _ = find_peaks(data[:, 1], distance=period_in_points)
        # data = data[peaks[0]:peaks[-1] + 1, :]

        # flipping flux data back over
        data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

        plt.plot(data[:, 1], color)
        plt.plot(peaks, data[:, 1][peaks], 'yx')

        # fname = join(binary_dir, file_name)
        plt.savefig(file_name)
        print "Plotted peaks, saved in: ", file_name

    @staticmethod
    def get_data_paths(data_dir, start_file, range, ext = '.csv'):
        from os import listdir
        from os.path import join
        '''
        :param data_dir: directory to a list of csv files, other things can be in the directory
        :param start_file: the file at which to start grabbing paths
        :param range: the number of files to grab after start_file
        :param ext: a string denoting the file extentions to grab
        :return: a list of file paths starting at start_file, and all subsequent paths in range listed after start file
                    designes specifically for chopped epoch data saved by self.save_epochs
        '''

        csv_files = list(filter(lambda x: ext in x, listdir(data_dir)))

        start_index = csv_files.index(start_file)
        end_index = start_index + range

        if end_index > len(csv_files):
            print 'end index was larger than the list of files, getting all files to the very end instead'
            end_index = len(csv_files) - 1

        csv_files = csv_files[start_index:end_index]

        csv_paths = [join(data_dir, csv_file) for csv_file in csv_files]

        return csv_paths


    @classmethod
    def merge_save_epochs(self, paths, save_path):
        from numpy import asarray, vstack, zeros
        '''
        :param paths: a list of paths to data files to merge
        :param save_path: a save path for the merged data
        :return: saves a merged set of epoch data
        '''

        merged_data = zeros(shape= (1,3))
        for i in range(len(paths)):
            data = self.get_LC_data_from_file(paths[i])
            if i == 0:
                merged_data = data
            else:
                merged_data = vstack((merged_data, data))

        self.save_LC_data_from_array(merged_data, save_path)
        print 'Saved merged data in ', save_path

    @classmethod
    def error_handler(self):
        import logging
        from logging.handlers import RotatingFileHandler
        from os.path import join
        from os import getcwd
        log_filename = join(getcwd(), 'log.txt')
        logging.basicConfig(filename= log_filename, level= logging.DEBUG)


    @classmethod
    def build_compute_save_multiple(self, binary_dirs):
        '''
        :param dirs: a list of system direcories
        :return: analyzes and saves the data for multiple star sytems
        '''

        import logging
        from datetime import datetime
        from os.path import join
        from os import getcwd
        log_filename = join(getcwd(), 'log.txt')
        logging.basicConfig(filename= log_filename, level= logging.DEBUG)
        logging.debug('This logger file everytime the logging protocall is called. Save this text to another file'
                      ' if you need it for later. ')

        for binary_dir in binary_dirs:

            b = self.load_binary(binary_dir)
            try:
                self.build_compute_save(b)
            except Exception:
                print 'Encountered error for system: ',b.system_name, ' at ', datetime.now()
                logging.exception('Encountered error for system: ',b.system_name, ' at ', datetime.now())
                continue

    @staticmethod
    def average_lightcurve_files(data_dir, file_names, save_path):
        '''
        :param dir: a directory to the list of data files to be averaged together
        :param file_names: a list of filenames to average together
        :param save_path: a save path for the averaged lightcurve
        :return: saves an averaged lightcurve
        '''

        from os import listdir
        from os.path import join

        #putting arrays in a list
        data_sets = [Processes.get_LC_data_from_file(join(data_dir, f)) for f in listdir(data_dir) if f in file_names]

        #making sure everything is the same shape as the first dataset in the list
        first_shape = data_sets[0].shape
        data_sets = [data_set[:first_shape[0], :] for data_set in data_sets]

        averaged_data = []
        averaged_time = []

        #for i in range(len())
    @staticmethod
    def list_subdir_paths(a_dir):
        '''
        :param a_dir: a directory
        :return: looks in a directory and returns a list of complete files paths withing that directory
        '''
        from os.path import join
        from os import walk

        subdirs =  [join(a_dir, subdir) for subdir in [x[0] for x in walk(a_dir)]]

        if a_dir in subdirs:
            del subdirs[subdirs.index(a_dir)]

        return subdirs



















#21.371629 - 33.283291


from os import getcwd
from os.path import join

cwd = getcwd()
#binary_dir = join(cwd, 'Outputs', 'Results', 'EPIC_201325107_anomaly_phi0_1')
# Processes.plot_peaks(LC_path, period = 0.1, start = 2795, stop = 2799, color = 'bx')
# LC_path = join(binary_dir, 'light_curve_data.csv')
# dir = join(cwd, 'data', 'KIC_201325107_epochs_secondary_minima')
# LC_path = join(cwd, 'data', 'Systems', 'KIC_201325107_flux.csv')
# epochs = Processes.chop_LC_from_file(LC_path, period = 0.09)
'''
binary_dir = join(cwd, 'data', 'Systems', 'EPIC_201325107_epoch')

b = Processes.load_binary(binary_dir)
Processes.build_compute_save(b)

binary_dir = join(cwd, 'data', 'Systems', 'EPIC_201325107_epoch(2)')

b = Processes.load_binary(binary_dir)
Processes.build_compute_save(b)
'''


bin_dir = join(cwd, 'Outputs', 'Results', 'constrained_SLSQP', 'EPIC_201325107_anomaly_phi0')
b = Processes.load_binary(bin_dir)
Processes.build_compute_save(b)

#a_dir = join(cwd, 'data', 'Systems', 'Batch2')
#binary_dirs = Processes.list_subdir_paths(a_dir)
#Processes.build_compute_save_multiple(binary_dirs)

#binary_dir = join(cwd, 'Outputs', 'Results', 'EPIC_201325107_epoch_SLSQP_del_1')
#Processes.analyze_solved_binary(binary_dir)

#data_files = ['epoch_at_21.371629000000212.csv', 'epoch_at_22.12760100000014.csv', 'epoch_at_22.863142000000153.csv']
#dir = join(cwd, 'data', 'KIC_201325107_epochs_primary_minima')
#Processes.average_lightcurve_files(dir, data_files, save_path= None)



'''
dir = join(cwd, 'data', 'KIC_201325107_epochs_primary_minima')
save_path = join(cwd, 'data', 'KIC_201325107_anomaly.csv' )
csv_paths = Processes.get_data_paths(dir, 'epoch_at_42.5592700000002.csv', range = 7)

Processes.merge_save_epochs(csv_paths, save_path)

Processes.plot_peaks(save_path, period= 10000, file_name= 'anomaly.png')
'''



