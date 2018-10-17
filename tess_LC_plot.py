'''Module containg various useful functions for Ed's exoplanet/TESS/NGTS PhD'''

from __future__ import division
import batman
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize
import pandas
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
	
def TIC_byID(ID):
	
	from astroquery.mast import Catalogs
	
	catTable = Catalogs.query_criteria(ID=ID, catalog="Tic")
	return catTable
		
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

	FLUX = DATA_WHOLE['PDCSAP_FLUX']  		#PDC corrected flux from target star
	FLUX_ERR = DATA_WHOLE['PDCSAP_FLUX_ERR']#Error in PDC corrected flux
	raw_flux = DATA_WHOLE['SAP_FLUX']		#Simple Aperture Photometry flux from target star
	
	#Load Quality flags in to remove flagged data points
	flags = DATA_WHOLE['QUALITY']
	flag_indices = np.where(flags > 0)

	flux_flagremoved = np.delete(FLUX, flag_indices)
	fluxerr_flagremoved = np.delete(FLUX_ERR, flag_indices)
	time_flagremoved = np.delete(time, flag_indices)
	raw_flux_flagremoved = np.delete(raw_flux, flag_indices)

	#Remove time points during central gap
	null_indices = np.where(np.isnan(time_flagremoved))
	time_nullremoved = np.delete(time_flagremoved, null_indices)
	flux_nullremoved = np.delete(flux_flagremoved, null_indices)
	fluxerr_nullremoved = np.delete(fluxerr_flagremoved, null_indices)
	raw_flux_nullremoved = np.delete(raw_flux_flagremoved, null_indices)
	
	return time_nullremoved, flux_nullremoved, fluxerr_nullremoved, time, raw_flux_nullremoved
	
def tess_LC_dataload_qlp(file_name):
	
	#Load FITS file
	hdul = fits.open(file_name)           		#Read in FITS file
	hdr = hdul[0].header                  		#Primary FITS header
	DATA_WHOLE = hdul[1].data             		#Extracts whole data

	#Extract desired columns
	time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
	FLUX = DATA_WHOLE['SAP_FLUX']
	
	#Load Quality flags in to remove flagged data points
	flags = DATA_WHOLE['QUALITY']
	flag_indices = np.where(flags > 0)

	flux_flagremoved = np.delete(FLUX, flag_indices)
	#fluxerr_flagremoved = np.delete(FLUX_ERR, flag_indices)
	time_flagremoved = np.delete(time, flag_indices)

	#Remove time points during central gap
	null_indices = np.where(np.isnan(time_flagremoved))
	time_nullremoved = np.delete(time_flagremoved, null_indices)
	flux_nullremoved = np.delete(flux_flagremoved, null_indices)
	
	return time_nullremoved, flux_nullremoved

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
	
def normalise_LC(flux, phase, period, Tdur):

	#Now we want to mask out the transit points, to do statistices on the rest of the LC
	transit_indices = np.where(np.abs(phase * period) <= T_dur / (2 * 24))	#Array indices of all phase/flux values during the transit
	FLUX_OOT = np.delete(flux, transit_indices)				#"Out Of Transit" flux values

	median = np.median(FLUX_OOT)										#median of all out-of-transit flux values

	flux_normalised = flux / median
	
	return flux_normalised
	
def plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, save=False):
	
	fig = plt.figure(figsize=(30, 20))
	
	#Top Subplot: Raw SAP Flux
	ax1 = fig.add_subplot(311)
	ax1.plot(time, raw_flux, 'ko', markersize=1.5)
	ax1.set_ylabel('Raw SAP Flux [e$^-$ / s]', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font) 
	ax1.tick_params(labelsize=20)
	ax1.set_title('T = {} ; Rs = {} ;  Period: {} days \n TOI: {} ;  TIC ID: {}'.format(Tmag, Rs, period, TOI, TIC), **axis_font)


	#Middle subplot: Unfolded LC
	ax2 = fig.add_subplot(312)
	ax2.plot(time, flux_normalised, 'ko', markersize='1.5')
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax2.tick_params(labelsize=20)
	#Bottom subplot: Folded LC
	ax3 = fig.add_subplot(313)

	ax3.plot(phase, flux_normalised, 'bo', markersize='1.5')
	#ax3.plot(phase_ordered, flux_model, 'r-')
	ax3.set_ylabel('Relative Flux', **axis_font)
	ax3.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax3.tick_params(labelsize=20)
	ax3.set_xlabel('Phase', **axis_font)
	
	plt.tight_layout()
	
	if save == True:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc.png'.format(TOI, TIC))
		plt.close()
	
	else:
		plt.show()

def plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs, save=False):
	
	fig = plt.figure(figsize=(30, 20))
				
	#Top subplot: Unfolded LC
	ax1 = fig.add_subplot(211)
	ax1.plot(time, flux, 'ko', markersize='3')
	ax1.set_ylabel('Relative Flux', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax1.tick_params(labelsize=20)
	ax1.set_title('T = {} ; Rs = {} ;  Period: {} days \n TOI: {} ;  TIC ID: {}'.format(Tmag, Rs, period, TOI, TIC), **axis_font)

	#Bottom subplot: Folded LC
	ax2 = fig.add_subplot(212)
	ax2.plot(phase, flux, 'bo', markersize='3')
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.tick_params(labelsize=20)
	ax2.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax2.set_xlabel('Phase', **axis_font)
				
	plt.tight_layout()
	
	if save == True:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc.png'.format(TOI, TIC))
	
	else:
		plt.show()

if __name__ == "__main__":
		
	parser = argparse.ArgumentParser()
	parser.add_argument('-fn', '--filename', type=str, nargs='*')
	parser.add_argument('-s', '--save', action='store_true')
	
	args = parser.parse_args()
	
	filenames = args.filename
	save = args.save
	
	axis_font = {'fontname':'DejaVu Sans', 'size':'30'}
	df = pandas.read_csv('/home/astro/phrvdf/tess_data_alerts/TOIs_20181016.csv', index_col='tic_id')		#.csv file containing info on parameters (period, epoch, ID, etc.) of all TOIs
	length = len(df.iloc[0])
		
	for i in range(len(filenames)):
		
		hdul = fits.open(filenames[i])
		hdr = hdul[0].header
		
		if hdr['ORIGIN'] == 'NASA/Ames':
		
			TIC = hdr['TICID']							#Loads TIC ID number from the fits header
			Tmag = hdr['TESSMAG']						#Loads TESS mag from fits header
			Rs = hdr['RADIUS']							#Loads stellar radius [solar radius] from fits header
			
			df2 = df.loc[TIC]
			
			if len(df2) == length:
				
				epoch = df2.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
				period = df2.loc['Period']			  	#Orbital Period [days]
				T_dur = df2.loc['Duration']				#Transit duration [hours]
				TOI = df2.loc['toi_id']            	  	#TIC ID for the object - used for plot title
				comments = df2.loc['Comment']			#Any existing comments on the object
			
				print "Epoch of first transit is {} [BJD - 2457000]".format(epoch)
				print "Orbital period is {} days".format(period)
				print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
				print "Existing comments on this object are: {}".format(comments)
						
				time, flux, fluxerr, time_whole, raw_flux = tess_LC_dataload_spoc(filenames[i])
			
				phase, phase_days = phase_fold(time, epoch, period)
			
				flux_normalised = normalise_LC(flux, phase, period, T_dur)
				
				if save == True:
					plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, save=True)
				
				elif save == False:
					plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, save=False)

				print TIC, TOI
				
			else:
				for j in range(len(df2)):
					
					df3 = df2.iloc[j]
					
					epoch = df3.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
					period = df3.loc['Period']			  	#Orbital Period [days]
					T_dur = df3.loc['Duration']				#Transit duration [hours]
					TOI = df3.loc['toi_id']            	  	#TIC ID for the object - used for plot title
					comments = df3.loc['Comment']			#Any existing comments on the object
			
					print "Epoch of first transit is {} [BJD - 2457000]".format(epoch)
					print "Orbital period is {} days".format(period)
					print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
					print "Existing comments on this object are: {}".format(comments)
						
					time, flux, fluxerr, time_whole, raw_flux = tess_LC_dataload_spoc(filenames[i])
			
					phase, phase_days = phase_fold(time, epoch, period)
			
					flux_normalised = normalise_LC(flux, phase, period, T_dur)
				
					if save == True:
						plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, save=True)
				
					elif save == False:
						plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, save=False)

					print TIC, TOI
			
		elif hdr['ORIGIN'] == 'MIT/QLP':
			
			TIC = np.int(hdr['TIC'])
			catTable = TIC_byID(TIC)
			
			Tmag, Rs = catTable['Tmag'][0], catTable['rad'][0]
			
			df2 = df.loc[TIC]

			
			if len(df2) == length:
				
				epoch = df2.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
				period = df2.loc['Period']			  	#Orbital Period [days]
				T_dur = df2.loc['Duration']				#Transit duration [hours]
				TOI = df2.loc['toi_id']            	  	#TIC ID for the object - used for plot title
				comments = df2.loc['Comment']			#Any existing comments on the object
			
				print "Epoch of first transit is {} [BJD - 2457000]".format(epoch)
				print "Orbital period is {} days".format(period)
				print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
				print "Existing comments on this object are: {}".format(comments)
			

				time, flux = tess_LC_dataload_qlp(filenames[i])
			
				phase, phase_days = phase_fold(time, epoch, period)
				
				flux_normalised = normalise_LC(flux, phase, period, T_dur)
				
				if save == True:
					plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs, save=True)
					
				elif save == False:
					plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs, save=False)
				
				print TIC, TOI
		
			else:
				for j in range(len(df2)):
					
					df3 = df2.iloc[j]
					
					epoch = df3.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
					period = df3.loc['Period']			  	#Orbital Period [days]
					T_dur = df3.loc['Duration']				#Transit duration [hours]
					pipeline = df3.loc['src']				#Pipeline used to reduce data - so that can call correct columns
					TOI = df3.loc['toi_id']            	  	#TIC ID for the object - used for plot title
					comments = df3.loc['Comment']			#Any existing comments on the object
			
					print "Epoch of first transit is {} [BJD - 2457000]".format(epoch)
					print "Orbital period is {} days".format(period)
					print "Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.)
					print "Pipeline used to process data is {}".format(pipeline)
					print "Existing comments on this object are: {}".format(comments)
						
					time, flux, time_whole, raw_flux = tess_LC_dataload_qlp(filenames[i])
			
					phase, phase_days = phase_fold(time, epoch, period)
			
					flux_normalised = normalise_LC(flux, phase, period, T_dur)
				
					if save == True:
						plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs, save=True)
					
					elif save == False:
						plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs, save=False)
				
					print TIC, TOI
					
				

	
	
			
			

		
		
	
	
	
	
	
	
		




		
	
	
	
	
