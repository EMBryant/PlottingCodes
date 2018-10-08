'''Code to use quality flags to remove "bad" points from TESS spoc LCs'''

from __future__ import division
import numpy as np
from astropy.io import fits
import pandas
import sys
import matplotlib.pyplot as plt

file_name = sys.argv[1]
TOI = float(sys.argv[2])

df = pandas.read_csv('../TESS/TOIs_Sec1_20180905.csv', index_col='toi_id')


epoch = df.loc[TOI, 'Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
period = df.loc[TOI, 'Period']			  	#Orbital Period [days]
T_dur = df.loc[TOI, 'Duration']				#Transit duration [hours]
pipeline = df.loc[TOI, 'src']				#Pipeline used to reduce data - so that can call correct columns
TIC_ID = df.loc[TOI, 'tic_id']              #TIC ID for the object - used for plot title


#Load FITS file
hdul = fits.open(file_name)           		#Read in FITS file
hdr = hdul[0].header                  		#Primary FITS header
DATA_WHOLE = hdul[1].data             		#Extracts whole data

#Extract desired columns
time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
#time_corr = DATA_WHOLE['TIMECORR']    		#Time correction: time - time_corr gives light arrival time at spacecraft


FLUX = DATA_WHOLE['PDCSAP_FLUX']  		#PDC corrected flux from target star

#Clean up "bad" points:

#remove zero entries
zero_entries = np.where(FLUX == 0)				#Locate entries in the FLUX column which have a value of 0
FLUX_zeroremoved = np.delete(FLUX, zero_entries)	#Remove corresponding entries from both FLUX and TIME columns
time_zeroremoved = np.delete(time, zero_entries)

#remove null entries - time
null_entries_time = np.where(np.isnan(time_zeroremoved))						#Locate entries in the TIME column which have a value of 'nan'
FLUX_nullremoved_intermediate = np.delete(FLUX_zeroremoved, null_entries_time)	#Remove corresponding entries from both FLUX and TIME columns
time_nullremoved_intermediate = np.delete(time_zeroremoved, null_entries_time)

#remove null entries - flux
null_entries_flux = np.where(np.isnan(FLUX_nullremoved_intermediate))						#Locate entries in the TIME column which have a value of 'nan'
FLUX_nullremoved = np.delete(FLUX_nullremoved_intermediate, null_entries_flux)	#Remove corresponding entries from both FLUX and TIME columns
time_nullremoved = np.delete(time_nullremoved_intermediate, null_entries_flux)	
	
flags = DATA_WHOLE['QUALITY']

time_good = np.array([])
time_med = np.array([])
time_bad = np.array([])
flux_good = np.array([])
flux_med = np.array([])
flux_bad = np.array([])

for i in range(len(flags)):
	if flags[i] == 0:
		time_good = np.append(time_good, time[i])
		flux_good = np.append(flux_good, FLUX[i])
	
	elif flags[i] >= 128:

		time_bad = np.append(time_bad, time[i])
		flux_bad = np.append(flux_bad, FLUX[i])
	
	else:
		time_med = np.append(time_med, time[i])
		flux_med = np.append(flux_med, FLUX[i])

print np.where(np.isnan(time_good))

print time_good[9200], time_good[10040]

plt.figure()

plt.plot(time_good, flux_good, 'ko', markersize=1.5)
plt.plot(time_med, flux_med, 'ro', markersize=1.5)
plt.plot(time_bad, flux_bad, 'go', markersize=1.5)

plt.show()




