'''Code to read in and plot TESS Sector 2 SAP and PDC LCs from FITS files

Command line arguments needed: (1) File name; (2) TOI ID'''

#Preliminary imports
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import pandas
import eds_phd_module as epm
import itertools
import argparse

def tess_LC_dataload_spoc(file_name):
	'''Loads TESS LC data for a given object and uses quality flags etc. to remove bad points 
		'''
	#Load FITS file
	hdul = fits.open(file_name)           		#Read in FITS file
	hdr = hdul[0].header                  		#Primary FITS header
	DATA_WHOLE = hdul[1].data             		#Extracts whole data

	#Extract desired columns
	time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
	#time_corr = DATA_WHOLE['TIMECORR']    		#Time correction: time - time_corr gives light arrival time at spacecraft

	raw_flux = DATA_WHOLE['SAP_FLUX']		#Simple Aperture Photometry flux from target star
	
	#Load Quality flags in to remove flagged data points
	flags = DATA_WHOLE['QUALITY']
	#flag_indices = np.where(flags > 0)
	flag_indices = np.array([])
	flag_vals = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]

	flux_RW_desat_vals = np.array([])
	time_RW_desat_vals = np.array([])

	for i in range(len(flags)):
		if flags[i] != 0:
			result = [seq for k in range(len(flag_vals), 0, -1) for seq in itertools.combinations(flag_vals, k) if sum(seq) == flags[i]]
			if 32 not in result[0]:
				flag_indices = np.append(flag_indices, i)

	time_flagremoved = np.delete(time, flag_indices)
	raw_flux_flagremoved = np.delete(raw_flux, flag_indices)

	#Remove time points during central gap
	null_indices = np.where(np.isnan(time_flagremoved))[0]
	if len(null_indices) == 0:
		return time_flagremoved, raw_flux_flagremoved, flux_RW_desat_vals, time_RW_desat_vals
	else:
		time_nullremoved = np.delete(time_flagremoved, null_indices)
		raw_flux_nullremoved = np.delete(raw_flux_flagremoved, null_indices)
	
		return time_nullremoved, raw_flux_nullremoved, flux_RW_desat_vals, time_RW_desat_vals


#Perform a phase fold
def phase_fold(time, epoch, period, max_phase=0.75):
	'''Function to convert a given set of times to phases, based off of a zero phase time and a period.
	
		INPUTS:
			time:					numpy array containing the time values to be phase folded [days]
			epoch:					the time at which to define phase = 0.0 [days]
			period:					time period over which to perform the phase fold [days]
			max_phase:				max value of phase to have in the output array
			
		OUTPUTS:
			phase:					numpy array of same dimensions as 'time' containing the calculated phase values
			phase_days:				same as 'phase' but in units of days'''
		
	phase = np.zeros_like(time)					#Empty array to hold phase values

	#Perform 'Phase Fold'
	for i in range(len(phase)):
		phase[i] = (time[i] - epoch) / period  -  np.int((time[i] - epoch) / period)	#populates 'phase' array with time values converted to phases
			
		if phase[i] < 0:																#takes any phases initially < 0 - from time points < first epoch - and makes positive	
			phase[i] = phase[i] + 1
		
		if phase[i] > max_phase:														#puts all phases in range (max_phase - 1) <= phase <= max_phase. This makes plots look slightly nicer
			phase[i] = phase[i] - 1	

	phase_days = phase * period						#additional output phase array in units of days

	return phase, phase_days

def plot_RW_desat_period_spoc(time, flux, phase, TOI, TIC, period, Tmag, Rs, save=False):
	
	fig = plt.figure(figsize=(30, 20))
	
	#Top Subplot: Raw SAP Flux
	ax1 = fig.add_subplot(211)
	ax1.plot(time, flux, 'ko', markersize=1.5)
	ax1.set_ylabel('SAP Flux [e$^-$ / s]', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font) 
	ax1.tick_params(labelsize=20)
	ax1.set_title('T = {:.2f} ; Rs = {:.3f} ;  Period: {} days \n TOI: {} ;  TIC ID: {}'.format(Tmag, Rs, period, TOI, TIC), **axis_font)

	#Bottom subplot: Folded LC
	ax2 = fig.add_subplot(212)

	ax2.plot(phase, flux, 'bo', markersize='1.5')
	ax2.set_ylabel('SAP Flux [e$^-$ / s]', **axis_font)
	ax2.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax2.tick_params(labelsize=20)
	ax2.set_xlabel('Phase', **axis_font)
	
	plt.tight_layout()
	
	if save:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc.png'.format(TOI, TIC))
		plt.close()
	
	else:
		plt.show()
		

if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-fn', '--filename', type=str, nargs='*')
	parser.add_argument('-s', '--save', action='store_true')
	parser.add_argument('-sec', '--sector', type=int)
	
	args = parser.parse_args()
	
	filenames = args.filename
	save = args.save
	sector = args.sector
	
	axis_font = {'fontname':'DejaVu Sans', 'size':'30'}
	df = pandas.read_csv('/home/astro/phrvdf/tess_data_alerts/TOIs_20181016.csv', index_col='tic_id')		#.csv file containing info on parameters (period, epoch, ID, etc.) of all TOIs
	length = len(df.iloc[0])
	
	period = 2.5
	if sector == 1:
		epoch = 1342.19
	else:
		epoch = 1371.13
	
	for i in range(len(filenames)):
		
		hdr = fits.open(filenames[i])[0].header
		
		if hdr['ORIGIN'] == 'NASA/Ames':
			
			TIC = hdr['TICID']							#Loads TIC ID number from the fits header
			Tmag = hdr['TESSMAG']						#Loads TESS mag from fits header
			Rs = hdr['RADIUS']							#Loads stellar radius [solar radius] from fits header
			
			df2 = df.loc[TIC]
			
			if len(df2) == length:
				
				T_dur = df2.loc['Duration']				#Transit duration [hours]
				TOI = df2.loc['toi_id']            	  	#TIC ID for the object - used for plot title
				comments = df2.loc['Comment']			#Any existing comments on the object
			
				print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
				print "Existing comments on this object are: {}".format(comments)
				
				time, flux, flux_RW, time_RW = tess_LC_dataload_spoc(filenames[i])
				
				phase, phase_days = phase_fold(time, epoch, period)
				
				if save:
					plot_RW_desat_period_spoc(time, flux, phase, TOI, TIC, period, Tmag, Rs, save=True)

				elif save == False:
					plot_RW_desat_period_spoc(time, flux, phase, TOI, TIC, period, Tmag, Rs, save=False)

				print TIC, TOI
	
			else:
				for j in range(len(df2)):
					
					df3 = df2.iloc[j]
					
					T_dur = df3.loc['Duration']				#Transit duration [hours]
					TOI = df3.loc['toi_id']            	  	#TIC ID for the object - used for plot title
					comments = df3.loc['Comment']			#Any existing comments on the object
			
					print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
					print "Existing comments on this object are: {}".format(comments)
						
					time, flux, flux_RW, time_RW = tess_LC_dataload_spoc(filenames[i])
				
					phase, phase_days = phase_fold(time, epoch, period)
							
					if save == True:
						plot_RW_desat_period_spoc(time, flux, phase, TOI, TIC, period, Tmag, Rs, save=True)
				
					elif save == False:
						plot_RW_desat_period_spoc(time, flux, phase, TOI, TIC, period, Tmag, Rs, save=False)

					print TIC, TOI
			
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
