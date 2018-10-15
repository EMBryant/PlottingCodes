'''Code to read in and plot TESS LCs from FITS files

Command line arguments needed: (1) File name; (2) TOI ID'''

#Preliminary imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import pandas
import rogues

#Load command line arguments
file_name = sys.argv[1]               		#name of FITS file containing LC data
TOI = float(sys.argv[2])					#TOI ID of the object
file_name2 = sys.argv[3]

df = pandas.read_csv('../TESS/TOIs_Sec1_20180905.csv', index_col='toi_id')		#.csv file containing info on parameters (period, epoch, ID, etc.) of all TOIs


epoch = df.loc[TOI, 'Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
period = df.loc[TOI, 'Period']			  	#Orbital Period [days]
T_dur = df.loc[TOI, 'Duration']				#Transit duration [hours]
pipeline = df.loc[TOI, 'src']				#Pipeline used to reduce data - so that can call correct columns
TIC = df.loc[TOI, 'tic_id']              #TIC ID for the object - used for plot title
comments = df.loc[TOI, 'Comment']			#Any existing comments on the object

print "Epoch of first transit is {} [BJD - 2457000]".format(epoch)
print "Orbital period is {} days".format(period)
print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
print "Pipeline used to process data is {}".format(pipeline)
print "Existing comments on this object are: {}".format(comments)

#Load FITS file
hdul = fits.open(file_name)           		#Read in FITS file
hdr = hdul[0].header                  		#Primary FITS header
DATA_WHOLE = hdul[1].data             		#Extracts whole data

hdul2 = fits.open(file_name2)
data2 = hdul2[1].data
time2 = data2['TIME']
flux_raw2 = data2['SAP_FLUX']
flux2 = data2['PDCSAP_FLUX']
flux_err2 = data2['PDCSAP_FLUX_ERR']
flags2 = data2['QUALITY']

#Extract desired columns
time = np.append(DATA_WHOLE['TIME'], time2)			  		#Time [BJD - 2457000]
#time_corr = DATA_WHOLE['TIMECORR']    		#Time correction: time - time_corr gives light arrival time at spacecraft

if pipeline == 'spoc':
	FLUX = np.append(DATA_WHOLE['PDCSAP_FLUX'], flux2) 		#PDC corrected flux from target star
	FLUX_ERR = np.append(DATA_WHOLE['PDCSAP_FLUX_ERR'], flux_err2) #Error in PDC corrected flux
	raw_flux = np.append(DATA_WHOLE['SAP_FLUX'], flux_raw2)		#Simple Aperture Photometry flux from target star
if pipeline == 'qlp':
	FLUX = DATA_WHOLE['SAP_FLUX']

#Load Quality flags in to remove flagged data points
flags = np.append(DATA_WHOLE['QUALITY'], flags2)
flag_indices = np.where(flags > 0)

flux_flagremoved = np.delete(FLUX, flag_indices)
#fluxerr_flagremoved = np.delete(FLUX_ERR, flag_indices)
time_flagremoved = np.delete(time, flag_indices)


#Remove time points during central gap
null_indices = np.where(np.isnan(time_flagremoved))
time_nullremoved = np.delete(time_flagremoved, null_indices)
flux_nullremoved = np.delete(flux_flagremoved, null_indices)
#fluxerr_nullremoved = np.delete(fluxerr_flagremoved, null_indices)

#Perform a phase fold
phase, phase_days = rogues.phase_fold(time_nullremoved, epoch, period, 0.75)

#Now we want to mask out the transit points, to do statistices on the rest of the LC
transit_indices = np.where(np.abs(phase_days) <= T_dur / (2 * 24))	#Array indices of all phase/flux values during the transit
FLUX_OOT = np.delete(flux_nullremoved, transit_indices)				#"Out Of Transit" flux values
phase_OOT = np.delete(phase_days, transit_indices)					#"Out Of Transit" phase values [units of days]

sigma = np.std(FLUX_OOT)											#Standard deviation of out-of-transit flux values
median = np.median(FLUX_OOT)										#median of all out-of-transit flux values

check_indices = np.where(np.abs(flux_nullremoved - median) > 5*sigma)     #Indices of all flux values > +5sigma from median - includes points during transits!
outlier_indices = np.array([])
for i in range(len(check_indices[0])):
	if len(np.where(transit_indices[0] == check_indices[0][i])[0]) == 0:		  	#ie. if the index corresponds to a point not during a transit
		outlier_indices = np.append(outlier_indices, check_indices[0][i])		  	#All points > +5sigma from median - NOT including transit points

phase_cleaned = np.delete(phase, outlier_indices)									#Remove all 5sigma points from phase, flux, and time
FLUX_cleaned = np.delete(flux_nullremoved, outlier_indices)
#fluxerr_cleaned = np.delete(fluxerr_nullremoved, outlier_indices)
time_cleaned = np.delete(time_nullremoved, outlier_indices)


#Calculate the model Light Curve
#phase_ordered = np.sort(phase_cleaned)
#rp = 0.1105
#a = 23.01
#flux_model = epm.light_curve_model(phase_ordered, rp, a)

#cat_table = rogues.TIC_byID(TIC)
#Tmag = cat_table['Tmag']
#print Tmag
#Plot data

axis_font = {'fontname':'DejaVu Sans', 'size':'20'}

fig = plt.figure()
if pipeline == 'spoc':
	#Top Subplot: Raw SAP Flux
	ax1 = fig.add_subplot(311)
	ax1.plot(time, raw_flux, 'ko', markersize=1.5)
	ax1.set_ylabel('Raw SAP Flux [e$^-$ / s]', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font) 
	ax1.set_title('TOI: {} ;  TIC ID: {} ;  Period: {} days'.format(TOI, TIC, period), **axis_font)


	#Middle subplot: Unfolded LC
	ax2 = fig.add_subplot(312)
	ax2.plot(time_cleaned, FLUX_cleaned / median, 'ko', markersize='1.5')
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.set_xlabel('Time [BJD - 2457000]', **axis_font)

	#Bottom subplot: Folded LC
	ax3 = fig.add_subplot(313)

	ax3.plot(phase_cleaned, FLUX_cleaned / median, 'bo', markersize='1.5')
	#ax3.plot(phase_ordered, flux_model, 'r-')
	ax3.set_ylabel('Relative Flux', **axis_font)
	ax3.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax3.set_xlabel('Phase', **axis_font)

if pipeline == 'qlp':
	
	#Top subplot: Unfolded LC
	ax1 = fig.add_subplot(211)
	ax1.plot(time_cleaned, FLUX_cleaned / median, 'ko', markersize='1.5')
	ax1.set_ylabel('Relative Flux', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax1.set_title('TOI: {} ;  TIC ID: {} ;  Period: {} days'.format(TOI, TIC, period), **axis_font)

	#Bottom subplot: Folded LC
	ax2 = fig.add_subplot(212)

	ax2.plot(phase_cleaned, FLUX_cleaned / median, 'bo', markersize='1.5')
#	ax2.plot(t, flux_model, 'r-')
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax2.set_xlabel('Phase', **axis_font)

plt.show()

