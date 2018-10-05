'''Code to read in and plot TESS LCs from FITS files

Command line arguments needed: (1) File name; (2) First transit epoch; 
(3) Period; (4) Plot title; (5) Pipeline used'''

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
pipeline = sys.argv[5]						#Pipeline used - allows appropriate flux header label to be used

#Load FITS file
hdul = fits.open(file_name)           		#Read in FITS file
hdr = hdul[0].header                  		#Primary FITS header
DATA_WHOLE = hdul[1].data             		#Extracts whole data

#Extract desired columns
time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
#time_corr = DATA_WHOLE['TIMECORR']    		#Time correction: time - time_corr gives light arrival time at spacecraft

if pipeline == 'spoc':
	FLUX = DATA_WHOLE['PDCSAP_FLUX']  		#PDC corrected flux from target star
if pipeline == 'qlp':
	FLUX = DATA_WHOLE['SAP_FLUX']

#remove zero entries
zero_entries = np.where(FLUX == 0)						#Locate entries in the FLUX column which have a value of 0
FLUX_zeroremoved = np.delete(FLUX, zero_entries)	#Remove corresponding entries from both FLUX and TIME columns
time_zeroremoved = np.delete(time, zero_entries)

#remove null entries
null_entries = np.where(np.isnan(time_zeroremoved))						#Locate entries in the TIME column which have a value of 'nan'
FLUX_nullremoved = np.delete(FLUX_zeroremoved, null_entries)	#Remove corresponding entries from both FLUX and TIME columns
time_nullremoved = np.delete(time_zeroremoved, null_entries)

phase = np.zeros_like(time_nullremoved)					#Empty array to hold phase values

#Perform 'Phase Fold'
for i in range(len(phase)):
	phase[i] = (time_nullremoved[i] - epoch) / period  -  np.int((time_nullremoved[i] - epoch) / period) 
	
	if phase[i] < 0:
		phase[i] = phase[i] + 1
		
	if phase[i] > 0.75:
		phase[i] = phase[i] - 1

#Plot data
fig = plt.figure()

#Top subplot: Unfolded LC
ax1 = fig.add_subplot(211)
ax1.plot(time_zeroremoved, FLUX_zeroremoved, 'ko', markersize='1')
if pipeline == 'spoc':
	ax1.set_ylabel('Flux [e$^-$ / s]')
elif pipeline == 'qlp':
	ax1.set_ylabel('Relative Flux')

ax1.set_title(TITLE)
ax1.set_xlabel('Time [BJD - 2457000]')

#Bottom subplot: Folded LC
ax2 = fig.add_subplot(212)

ax2.plot(phase * period, PDC_FLUX_nullremoved, 'bo', markersize='1')
if pipeline = 'spoc':
	ax2.set_ylabel('Flux [e$^-$ / s]')
elif pipeline = 'qlp':
	ax2.set_ylabel('Relative Flux')
ax2.set_xlabel('Phase [days]')

plt.show()
