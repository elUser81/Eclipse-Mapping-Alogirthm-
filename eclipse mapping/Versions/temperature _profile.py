def normalize(self, array, over=None):

    '''normalizes an array over the max value if over = None, or over whatever over equals if over != None'''
    from numpy import asarray, amax
    len_arr = len(array)
    if not over:
        return asarray(map(lambda x: x / amax(array), array)).reshape((len_arr))
    else:
        return asarray(map(lambda x: x / over, array)).reshape((len_arr))


def flux_to_temp(self):
    from numpy import asarray, amax
    '''Takes flux image solutions from binary object 'self' and creates a corresponding set of temperature solutions
        and assigns them to self.temperature_solutions, which are later called in self.save_bbody_curve
        
        see line 783 of phoebe_scratch_work2_0.py for context
        primary use of this method is in Processephoebe_scratch_work2_0.pys.build_compute_save and Processes.analyze_solved_binary, and
        Processes.build_compute_save, called in lines 1104 and 1255 respectively in phoebe_scratch_work2_0.py
        '''

    image_solutions = self.image_solution.reshape(self.N**2)
    temp_primary = self.temperature_primary


    def temp_func(I):

        '''takes a pixel brightness solution and calculates its corresponding temperature solution'''
        steph_boltz = 5.67 * 10 ** -8
        flux_primary = steph_boltz * temp_primary ** 4
        T = (( I *flux_primary ) /steph_boltz) ** (0.25)
        return T/ temp_primary

    temp_solutions = map(temp_func, image_solutions)

    # normalizing.... this normalization might be an issue since temp_func already returns a ratio
    temp_solutions = asarray(map(lambda x: x / amax(temp_solutions), temp_solutions))
    # temp_solutions = asarray(map(lambda x: x * temp1, temp_solutions))

    #reshaping to same dims as original solution set
    self.temperature_solutions = temp_solutions.reshape((self.N, self.N))




def save_bbody_curve(self, path):
    '''
    :param self: binary object
    :param path: path to save temp profile
    :return: computes and saves temperature profile to path

    this method is defined at line 805 in phoebe_scratch_work2_0.py, called in Processes.build_compute_save, and
    Processes.analyze_solved_binary, lines 1150 and 1260 respectively

    '''
    from numpy import mean, asarray, floor, log
    import matplotlib.pyplot as plt

    bb_curve = [] #empty list to store mean temperature values
    N = self.N
    conv = self.pixel_to_RL1  # factor to convert pixel to RL1 units
    temp_solutions = self.temperature_solutions #temperature solutions from self.flux_to_temp()
    radius_primary = self.radius_primary

    '''taking a range from 0 to N/2 and rounding down with floor, in case N is odd.
        this is a set radii for the azimuthal temperature in pixel space. Then I multiply by conv (see line 42)
        and reverse the list. radii should be a list of radii in RL1 units going from innermost to outermost'''

    radii = asarray([i for i in reversed(map(lambda x: x * conv, range(int(floor(N / 2)))))])
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

    for n in range(int(floor(N / 2))):
        rng = range(N)
        pixel_ring = asarray([temp_solutions[x, y] for x in rng for y in rng if dist(x, y) == n])
        bb_curve.append(mean(pixel_ring))

    reved_bb_curve = reversed(bb_curve)
    bb_curve = asarray([i for i in reved_bb_curve])  ## need to reverse data order because of the way the loop iterates

    # through the data, reverse orders data from 0 to max radius

    def acc_temp_theory(R):
        # consts = (3 * G * M1 * M_trans) / (8 * pi * steph_boltz)
        # T_to_fourth = consts * (R**-3) * (1 - viscosity * (Rwd/R) ** (1/2))
        # just want the shape, so only using the stephon_boltzman constant

        # (((steph_boltz * R**3) ** (-1)), removing stephan boltzman for now
        # T = (((R_in_m**3) ** (-1)) * (1 - viscosity * ((radius_primary*(695.51 * 10**6))/R) ** (0.5))) ** (0.25)
        #comments just above are different things I tried

        #this is the ratio that we talked about in our meeting (assuming my math is correct...)
        term1 = (radius_primary / R) ** 0.75
        numer = 1 - viscosity * ((radius_primary / R) ** 0.5)
        denom = 1 - viscosity

        T_ratio = term1 * ((numer / denom) ** 0.25)

        return T_ratio

    theory_temps = asarray(map(acc_temp_theory, radii))
    # theory_temps = self.normalize(asarray(map(acc_temp_theory, radii)))
    radii = self.normalize(radii, over=self.RL1_from_primary) #normalizing radii from 0 to 1

    fig, ax = plt.subplots()

    # ln_T1 = map(lambda x: log(x), theory_temps) ---> this is what I got rid of, I couldnt find any other place where
                                                    # I called log, yet the profile still looks logarithmic
    # ln_T2 = map(lambda x: log(x), bb_curve)

    T1 = self.normalize(theory_temps) #maybe normalizing is the issue?
    T2 = self.normalize(bb_curve)

    ax.plot(radii, T1, 'rx', label='Theoretical')
    # ax.plot(radii, T2, 'b^', label = 'Reconstructed')
    # ax.legend()
    ax.set_xlabel('[R/RL1]')
    ax.set_ylabel('T/Twd')
    plt.title(''.join(['Temperature profile of ', self.system_name]))
    plt.savefig(path)
