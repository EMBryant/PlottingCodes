'''Python programme for quick plotting of NGTS lightcurves, both binned and unbinned, from the data text files'''

#Initial Imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

#Load Data file and extract JD, mag, flux
file_name = sys.argv[1]
DATA_WHOLE = np.loadtxt(file_name)
HJD = DATA_WHOLE[:, 0]                         #Heliocentric Julian Day
JD_MOD = HJD - 2450000                         #Modified Julian Day
MAG = DATA_WHOLE[:, 1]                         #magnitude
FLUX_REL = DATA_WHOLE[:, 3]                    #Relative Flux (=1 out of transit)
print sys.argv
#Phase Fold the data
epoch = float(sys.argv[2])                     #epoch of first transit and orbital period
period = float(sys.argv[3])                          

timebase = JD_MOD - epoch                      #sets time relative to first transit
phase = np.zeros_like(timebase)
for i in range(len(phase)):
    phase[i] = timebase[i] / period - np.int(timebase[i] / period)    #converts time to phase

#clean up phases - ie. set all phases in range (-0.25 to 0.75) 
for i in range(len(phase)):
    if phase[i] < 0:
        phase[i] = 1 + phase[i]
       
    if phase[i] > 0.75:
        phase[i] = phase[i] - 1

#Bin Data into bins of width 0.01 in phase
binning_data = np.zeros((len(phase), 3))		#Empty 2D array used to hold data for binning
binning_data[:, 0] = phase
binning_data[:, 1] = FLUX_REL

#This loop assigns each data point an integer value from -250 to 749, each integer corresponding to a phase bin
for i in range(len(phase)):
	if binning_data[i, 0] < 0:
		binning_data[i, 2] = np.int(binning_data[i, 0] * 1000) - 1		#Ensures correct binning for phases < 0
	
	else:
		binning_data[i, 2] = np.int(binning_data[i, 0] * 1000)

bin_means = np.zeros((1000, 2))				#Empty array to hold mean values of phase and flux for each bin
for i in range(len(bin_means[:, 0])):
	bin_means[i, 0] = np.mean(binning_data[np.where(binning_data[:, 2] == i - 250), 0])
	bin_means[i, 1] = np.mean(binning_data[np.where(binning_data[:, 2] == i - 250), 1])

#Relabel binned data sets
phase_binned = bin_means[:, 0]
flux_binned = bin_means[:, 1]


#Produce the plots
title = sys.argv[4]
axis_font = {'fontname':'DejaVu Sans', 'size':'20'}

plt.figure()

plt.plot(phase_binned, flux_binned, 'ko', markersize=2)
plt.ylabel('Relative Flux', **axis_font)
plt.xlabel('Phase', **axis_font)
plt.gca().set_xticks([-0.25, 0.0, 0.25, 0.50, 0.75])
plt.title(title, **axis_font)

plt.show()

fig = plt.figure()

ax1 = fig.add_subplot(211)

ax1.plot(phase, FLUX_REL, 'ko', markersize=1)
ax1.set_ylabel('Relative Flux', **axis_font)
#ax1.set_xlabel('Phase', **axis_font)
ax1.set_xticks([-0.25, 0.0, 0.25, 0.50, 0.75])
ax1.set_title(title, **axis_font)

ax2 = fig.add_subplot(212)

ax2.plot(phase_binned, flux_binned, 'ko', markersize=2)
ax2.set_ylabel('Relative Flux', **axis_font)
ax2.set_xlabel('Phase', **axis_font)
ax2.set_xticks([-0.25, 0.0, 0.25, 0.50, 0.75])

plt.show()
