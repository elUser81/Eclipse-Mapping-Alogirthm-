from phoebe_scratch_work import Binary
import matplotlib.pyplot as plt
from  mpl_toolkits.mplot3d import Axes3D
from os import getcwd
from os.path import join
import numpy as np
from phoebe_scratch_work import Binary
from phoebe_scratch_work import Processes as p
from phoebe_scratch_work import Binary

'''
params = {'q': 0.21, 'period': 0.09112, 'sma': 0.85, 'ecc':0.0,
        'distance': 7.1428571428571, 'incl': 88}
'''

'''
#found the parallax on SIMBAD and computed distance in pcs

test_b = Binary(params['incl'], params['period'], params['sma'], params['q'], params['ecc'], 25, params['distance'])
print(test_b.RL1_from_primary)
print(test_b.RL1_from_secondary)
test_b.build_orbital_model(phi0 = np.pi)


x1, y1, z1 = test_b.x_primary, test_b.y_primary, test_b.z_primary
x2, y2, z2 = test_b.x_secondary, test_b.y_secondary, test_b.z_secondary

phi = np.arange(0, 2*np.pi, 0.2)


fig = plt.figure()
ax = Axes3D(fig)

ax.plot(xs = x1(phi), ys = y1(phi), zs = z1(phi))
ax.plot(xs = x2(phi), ys = y2(phi), zs = z2(phi))


plt.plot(y1(phi), z1(phi))
plt.plot(y2(phi), z2(phi))

cwd = getcwd()
fname = join(cwd, 'TEST.png')

plt.savefig(fname)
'''

cwd = getcwd()
LC_path = join(cwd, 'data', 'KIC_201325107_flux.csv')

orbit_params = p.make_parameter_set('KIC_201325107', 76.19118, 0.7443, 5.00352, 0.5489, 0.0, 12.5, 1.643034)
light_curve_data = p.get_LCdata_from_file(LC_path)

#print light_curve_data[:, 2]

b = Binary(orbit_params, resolution= 25, chi_AIM=1,
           delta = 3, light_curve_data = light_curve_data, method = 'SLSQP')

p.build_compute_save(b)


