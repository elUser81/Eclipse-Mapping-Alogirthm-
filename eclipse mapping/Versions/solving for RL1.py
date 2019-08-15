from scipy.optimize import minimize, fsolve
import numpy as np
'''
('M1: ', 0.8202204949371832)
('M2: ', 0.17224630393680848)
('sma: ', 0.85)
'''



def to_kg(sol_mass):
    # converts_solar mass to kg
    # 1.9889200011445836e+30

    return sol_mass * 1.9889200011445836 * 10 ** 30



M1, M2 = map(to_kg, [0.8202204949371832, 0.17224630393680848]) #masses of primary and secondary stars in kg
sma = 0.85 * 1.496 * 10**11 #semi major axis in m

estimate = sma * (M2/(3*M1))**(1/3) #a good estimate from wikipedia
estimate = np.asarray([estimate])
x0 = np.asarray([0])
def get_RL1_func(M1, M2, sma):
    # defining args of the 1st lagrange point function

    def func(r):
        # this equation gives the roots of the RL1 function

        return abs((M2 / (r ** 2)) + M1 / (sma ** 2) - ((r * (M1 + M2)) / (sma ** 3)) - (M1 / ((sma - r) ** 2)))  # == 0

    return func



RL1_func = get_RL1_func(M1, M2, sma)



#res = minimize(fun = RL1_func, x0 = np.asarray[0], method= 'SLSQP')

res = minimize(RL1_func, x0 = x0, method= 'Nelder-Mead')

#print(RL1_func(res.x))
#print(res.x/( 1.496 * 10**11))

#IT WORKS IM CLEVER KANSAS IS A GREAT BAND YEET THIS WHEAT HA. 44444444444444444


