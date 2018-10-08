'''Code to read in and plot TESS LCs from FITS files

Command line arguments needed: (1) File name; (2) TOI ID'''

#Preliminary imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import pandas

#Load command line arguments
file_name = sys.argv[1]               		#name of FITS file containing LC data
TOI = float(sys.argv[2])

df = pandas.read_csv('../TESS/TOIs_Sec1_20180905.csv', index_col='toi_id')


epoch = df.loc[TOI, 'Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
period = df.loc[TOI, 'Period']			  	#Orbital Period [days]
T_dur = df.loc[TOI, 'Duration']				#Transit duration [hours]
pipeline = df.loc[TOI, 'src']				#Pipeline used to reduce data - so that can call correct columns
TIC_ID = df.loc[TOI, 'tic_id']              #TIC ID for the object - used for plot title
TITLE = "TOI ID: TOI; "				  		#Title for plot

print epoch, period, T_dur

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

#Clean up "bad" points:

if pipeline == 'spoc':
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
	
	

if pipeline == 'qlp':
	flag = DATA_WHOLE['QUALITY']
	bad_indices = np.where(flag == 1)
	FLUX_nullremoved = np.delete(FLUX, bad_indices)
	time_nullremoved = np.delete(time, bad_indices)
	
	FLUX_bad = FLUX[bad_indices]
	time_bad = time[bad_indices]

phase = np.zeros_like(time_nullremoved)					#Empty array to hold phase values

#Perform 'Phase Fold'
for i in range(len(phase)):
	phase[i] = (time_nullremoved[i] - epoch) / period  -  np.int((time_nullremoved[i] - epoch) / period) 
	
	if phase[i] < 0:
		phase[i] = phase[i] + 1
		
	if phase[i] > 0.75:
		phase[i] = phase[i] - 1	

phase_days = phase * period								#two additional phase arrays, one in units of days, the other hours
phase_hours = phase_days * 24

#Now we want to mask out the transit points, to do statistices on the rest of the LC
transit_indices = np.where(np.abs(phase_hours) <= T_dur / 2)	#Array indices of all phase/flux values during the transit
FLUX_OOT = np.delete(FLUX_nullremoved, transit_indices)			#"Out Of Transit" flux values
phase_OOT = np.delete(phase_days, transit_indices)				#"Out Of Transit" phase values [units of days]

sigma = np.std(FLUX_OOT)										#Standard deviation of out-of-transit flux values
median = np.median(FLUX_OOT)									#median of all out-of-transit flux values

check_indices = np.where(np.abs(FLUX_nullremoved - median) > 3*sigma)     #Indices of all flux values > +5sigma from median - includes points during transits!
print(len(transit_indices[0]), len(check_indices[0]))
outlier_indices = np.array([])
for i in range(len(check_indices[0])):
	if len(np.where(transit_indices[0] == check_indices[0][i])[0]) == 0:		  #ie. if the index corresponds to a point not during a transit
		outlier_indices = np.append(outlier_indices, check_indices[0][i])

print(len(outlier_indices))
phase_cleaned = np.delete(phase, outlier_indices)
FLUX_cleaned = np.delete(FLUX_nullremoved, outlier_indices)
time_cleaned = np.delete(time_nullremoved, outlier_indices)

axis_font = {'fontname':'Times New Roman', 'size':'20'}

#Test plot
plt.figure()

plt.plot(phase, FLUX_nullremoved / median, 'ro', markersize=1)
plt.plot(phase_cleaned, FLUX_cleaned / median, 'ko', markersize=1)

plt.xlabel('Phase', **axis_font)
plt.gca().set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
plt.ylabel('Relative Flux', **axis_font)


plt.show()



#Plot data
fig = plt.figure()

#Top subplot: Unfolded LC
ax1 = fig.add_subplot(211)
ax1.plot(time_cleaned, FLUX_cleaned / median, 'ko', markersize='1.5')
if pipeline == 'qlp':
	ax1.plot(time_bad, FLUX_bad / median, 'ro', markersize='1.5')
ax1.set_ylabel('Relative Flux', **axis_font)

ax1.set_title('TOI: {} ;  TIC ID: {} ;  Period: {} days'.format(TOI, TIC_ID, period), **axis_font)
ax1.set_xlabel('Time [BJD - 2457000]', **axis_font)

#Bottom subplot: Folded LC
ax2 = fig.add_subplot(212)

ax2.plot(phase_cleaned, FLUX_cleaned / median, 'bo', markersize='1.5')
ax2.set_ylabel('Relative Flux', **axis_font)
ax2.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
ax2.set_xlabel('Phase', **axis_font)

plt.show()

