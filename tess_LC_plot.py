'''Code to read in and plot TESS LCs from FITS files

Command line arguments needed: (1) File name; (2) First transit epoch; 
(3) Period; (4) Plot title'''

#Preliminary imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys

#Load command line arguments
file_name = sys.argv[1]               		#name of FITS file containing LC data
epoch = float(sys.argv[2])		      		#Time of first transit centre [BJD - 2457000]
period = float(sys.argv[3])			  		#Orbital Period [days]
TITLE = sys.argv[4]					  		#Title for plot

#Load FITS file
hdul = fits.open(file_name)           		#Read in FITS file
hdr = hdul[0].header                  		#Primary FITS header
DATA_WHOLE = hdul[1].data             		#Extracts whole data

#Extract desired columns
time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
#time_corr = DATA_WHOLE['TIMECORR']    		#Time correction: time - time_corr gives light arrival time at spacecraft
PDC_FLUX = DATA_WHOLE['PDCSAP_FLUX']  		#PDC corrected flux from target star

#remove zero entries
zero_entries = np.where(PDC_FLUX == 0)						#Locate entries in the FLUX column which have a value of 0
PDC_FLUX_zeroremoved = np.delete(PDC_FLUX, zero_entries)	#Remove corresponding entries from both FLUX and TIME columns
time_zeroremoved = np.delete(time, zero_entries)

#remove null entries
null_entries = np.where(np.isnan(time_zeroremoved))						#Locate entries in the TIME column which have a value of 'nan'
PDC_FLUX_nullremoved = np.delete(PDC_FLUX_zeroremoved, null_entries)	#Remove corresponding entries from both FLUX and TIME columns
time_nullremoved = np.delete(time_zeroremoved, null_entries)

phase = np.zeros_like(time_nullremoved)					#Empty array to hold phase values

#Perform 'Phase Fold'
for i in range(len(phase)):
	phase[i] = (time_nullremoved[i] - epoch) / period  -  np.int((time_nullremoved[i] - epoch) / period) 
	
	if phase[i] < 0:
		phase[i] = phase[i] + 1
		
	if phase[i] > 0.75:
		phase[i] = phase[i] - 1

median_flux_value = np.median(PDC_FLUX_nullremoved)

#Plot data
fig = plt.figure()

#Top subplot: Unfolded LC
ax1 = fig.add_subplot(211)

ax1.plot(time_zeroremoved, PDC_FLUX_zeroremoved / median_flux_value, 'ko', markersize='1')
ax1.set_title(TITLE)
ax1.set_ylabel('Relative Flux')
ax1.set_xlabel('Time [BJD - 2457000]')

#Bottom subplot: Folded LC
ax2 = fig.add_subplot(212)

ax2.plot(phase * period, PDC_FLUX_nullremoved / median_flux_value, 'bo', markersize='1')
ax2.set_ylabel('Relative Flux')
ax2.set_xlabel('Phase [days]')

plt.show()
