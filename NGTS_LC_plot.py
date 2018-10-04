'''Python programme for quick plotting of NGTS lightcurves, both binned and unbinned, from the data text files'''

#Initial Imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

#Load Data file and extract JD, mag, flux
file_name = sys.argv[1]
DATA_WHOLE = np.loadtxt(file_name)
HJD = DATA_WHOLE[:, 0]                        #Heliocentric Julian Day
JD_MOD = HJD - 2450000                        #Modified Julian Day
MAG = DATA_WHOLE[:, 1]                        #magnitude
FLUX_REL = DATA_WHOLE[:, 3]                   #Relative Flux (=1 out of transit)

#Phase Fold the data
epoch = float(sys.argv[2])                    #epoch of first transit and orbital period
period = float(sys.argv[3])                          

#print(sys.argv)
#print(epoch)
#print(period)

timebase = JD_MOD - epoch                     #sets time relative to first transit
phase = np.zeros_like(timebase)
for i in range(len(phase)):
    phase[i] = timebase[i] / period - np.int(timebase[i] / period)    #converts time to phase

#clean up phases
for i in range(len(phase)):
    if phase[i] < 0:
        phase[i] = 1 + phase[i]
       
    elif phase[i] > 0.8:
        phase[i] = phase[i] - 1

#Bin the data
bin_size = int(float(sys.argv[4]))
binned_dataset_length = np.int(len(phase) / bin_size)

binned_phase = np.zeros(binned_dataset_length)
binned_flux = np.zeros(binned_dataset_length)

for j in range(binned_dataset_length):
    binned_phase[j] = np.mean(phase[j * bin_size : (j+1) * bin_size])
    binned_flux[j] = np.mean(FLUX_REL[j * bin_size : (j+1) * bin_size])

#Produce the plots
title = sys.argv[5]
axis_font = {'fontname':'Times New Roman', 'size':'20'}

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.plot(phase, FLUX_REL, 'ko')
ax1.set_ylabel('Relative Flux', **axis_font)
ax1.set_title(title, **axis_font)

ax2 = fig.add_subplot(212)
ax2.plot(binned_phase, binned_flux, 'bo')
ax2.set_xlabel('Phase', **axis_font)
ax2.set_ylabel('Relative Flux', **axis_font)

plt.show()
