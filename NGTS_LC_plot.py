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
print(sys.argv)
#Phase Fold the data
epoch = float(sys.argv[2])                    #epoch of first transit and orbital period
period = float(sys.argv[3])                          

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

#Produce the plots
title = sys.argv[4]
axis_font = {'fontname':'Times New Roman', 'size':'20'}

fig = plt.figure()

ax1 = fig.add_subplot(211)
ax1.plot(JD_MOD, FLUX_REL, 'ko', markersize=1)
ax1.set_ylabel('Relative Flux', **axis_font)
ax1.set_xlabel('Time [HJD - 2450000]')

ax2 = fig.add_subplot(212)
ax2.plot(phase * period, FLUX_REL, 'ko', markersize=1)
ax2.set_ylabel('Relative Flux', **axis_font)
ax2.set_xlabel('Phase [days]')
ax1.set_title(title, **axis_font)

plt.show()
