
from numpy import genfromtxt
from os import getcwd
from os.path import join
import phoebe
from numpy import nan_to_num, inf, nan
from os.path import join
from os import listdir

from os import getcwd
import os
import datetime
from tomography import Processes, Binary
from os.path import join
import numpy as np



pdir = os.path.abspath(os.path.join(getcwd(), os.pardir))
base_dir = os.path.join(pdir, 'Outputs', 'Results', 'constrained_SLSQP','EPIC_201325107_anomaly_phi0')
base = Processes.reload_binary(base_dir)

signal_dir = join(pdir, 'Outputs', 'Results', 'EPIC_201325107_epoch_anomoly_del_1')
signal = Processes.reload_binary(signal_dir)

save_dir = join(pdir, 'Outputs', 'Convolutions')
Processes.convolve(base, signal, save_dir)




















































'''
system_name = 'EPIC_201325107(1)'
binary_dir = join(cwd, 'Outputs', 'Results', system_name)
params_path = join(binary_dir, 'parameters.csv')
LC_path = join(binary_dir, 'light_curve_data.csv')

data_path = join(cwd, 'data', 'Systems', 'KIC_201325107_flux.csv')
from numpy import amin, asarray, where, ceil, floor, diff, mean
import matplotlib.pyplot as plt
from pandas import read_csv
from scipy.signal import find_peaks



def get_data_from_txt(path):
    from numpy import loadtxt

    data = loadtxt(path)

    return data

def mag_to_flux(mag):
    m1 = 5
    flux_ratio = 10 ** ((m1 - mag)/2.5)
    return flux_ratio




data = asarray(read_csv(data_path).values)

data[:, 0] = asarray(map(lambda x: x - amin(data[:, 0]), data[:, 0]))
data[:, 1] = asarray(map(lambda x: (-1)*x, data[:, 1]))

time_diff = data[:, 0][1] - data[:, 0][0]

index_start = 1250
index_end = 1500

#flux, time = data[:, 1][index_start:index_end], data[:, 0][index_start:index_end]

data = data[index_start:index_end, :]
flux, time = data[:, 1], data[:, 0]

print diff(time)

period = 0.09113
distance = ceil(floor(period/time_diff))

#print 'time difference: ', distance


peaks, _ = find_peaks(flux, distance = distance)
avg =  mean(diff(peaks))
print 'average time: ', avg
#print 'avg differences btw peaks: ',diff(peaks)

#peaks, _ = find_peaks(flux, distance = float(avg))

flux = asarray(map(lambda x: (-1)*x, flux))


print peaks

plt.plot(flux, 'bx')
plt.plot(peaks, flux[peaks], 'yx')
plt.savefig(join(cwd, 'SYNTHETIC1.png'))
print 'Plot Saved'


def chop_LC(b, stop = None, chop_at = 'all', time_error = 0.1):
    #chops LC data into epochs, assumes time incriment of observation is consistant throughout dataset
    from scipy.signal import find_peaks
    from numpy import floor, mean, diff, asarray

    epochs = []
    data = b.light_curve_data

    #subtracting smallest time to start time at 0
    data[:, 0] = asarray(map(lambda x: x - amin(data[:, 0]), data[:, 0]))
    #flipping flux data set to find peaks as primary minima
    data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))

    #period time in units of index
    period = floor(b.period / data[:, 0][1] - data[:, 0][0])
    #correcttime between observations
    t_diff = data[1:0] - data[0:0]

    #peaks is a list of indicies where each index is a peak
    peaks, _ = find_peaks(data[:, 1], distance = period)
    data = data[peaks[0]:peaks[-1] + 1, :]

    #flipping flux data back over
    data[:, 1] = asarray(map(lambda x: (-1) * x, data[:, 1]))



    for i in range(len(peaks - 2)):
        if chop_at == 'all':
            if i == stop:
                break
            elif data[0: i + 1] - data[0: i] >= time_diff * time_error:
                pass
            else:
                start_chop, stop_chop = peaks[i], peaks[i+1]

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

'''

'''
fig, ax = plt.subplots()

data1 = [i for i in range(10)]
data2 = [i for i in reversed(data1)]
times = [i for i in range(10)]

print data1
print data2
'''




