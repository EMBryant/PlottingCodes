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

def light_curve_model(t, rp, a, t0 = 0., per = 1., inc = 90., ecc = 0., w = 90., limb_dark = "uniform", u = []):
	'''Python Function to use the batman package to produce a model lightcurve using the following input parameters:
			t: 						numpy array containing the time values at which to calculate the model LC
			rp: 					radius of the planet, in units of stellar radii
			a						semi-major axis of the planet's orbit, in units of stellar radii
			t0 = 0.					time of the transit centre
			per = 1.				orbital period
			inc = 90.				orbital inclination [degrees]
			ecc = 0.				orbital eccentricity
			w = 90. 				longitude of periastron [degrees]
			limb_dark = "uniform"	limb darkening mechanism 
			u = []					limb darkening coefficients
			
		NOTE: Units of t, t0, and per are not set to be a specific unit but must all be consistent
		
		OUTPUT:
			flux_model :  a numpy array of the same dimensions as t, containing the relative for the model LC at each 
							time point in t'''

	params = batman.TransitParams()					#a batman instance containing all the input parameters
	params.t0 = t0									#time of the transit centre
	params.per = per								#orbital period
	params.rp = rp									#radius of the planet, in units of stellar radii
	params.a = a									#semi-major axis of the planet's orbit, in units of stellar radii
	params.inc = inc								#orbital inclination [degrees]
	params.ecc = ecc								#orbital eccentricity
	params.w = w									#longitude of periastron [degrees]
	params.u = u									#limb darkening coefficients
	params.limb_dark = limb_dark					#limb darkening mechanism 
	
	m = batman.TransitModel(params, t)				#Initialises the transit model
	flux_model = m.light_curve(params)				#Calculates the light curve flux
	
	return flux_model
		
def light_curve_chi_squared_calc(params_initial, t, FLUX, FLUX_ERR):
	'''Python Function to calculate the chi_squared value for a given model light curve found using the batman package given the initial parameter inputs:
	
		INPUTS:
			params_initial:			numpy array containing initial parameter guesses
			t: 						numpy array containing the time values at which to calculate the model LC and test the parameters
			FLUX:					numpy array containing the flux data value at each point in 't'
			FLUX_ERR:				numpy array containing the errors in each value in 'FLUX'
						
		NOTE: Units of t, t0, and per are not set to be a specific unit but must all be consistent
		
		OUTPUT:
			chi_2 :		a chi squared values for the parameters'''
	
	params = batman.TransitParams()					#a batman instance containing all the input parameters
	params.t0 = 0.									#time of the transit centre
	params.per = 1.									#orbital period
	params.rp = params_initial[0]					#radius of the planet, in units of stellar radii
	params.a = params_initial[1]					#semi-major axis of the planet's orbit, in units of stellar radii
	params.inc = 90.					#orbital inclination [degrees]
	params.ecc = 0.					#orbital eccentricity
	params.w = 90.					#longitude of periastron [degrees]
	params.u = []									#limb darkening coefficients
	params.limb_dark = "uniform"					#limb darkening mechanism 
			
	m = batman.TransitModel(params, t)				#Initialises the transit model
	flux_model_initial = m.light_curve(params)		#Calculates the light curve flux
	
	chis_2 = np.zeros_like(t)
	for i in range(len(t)):
		chis_2[i] = (FLUX[i] - flux_model_initial[i])**2 / FLUX_ERR[i]**2
	
	chi_2 = np.sum(chis_2) / (len(chis_2) - 1)
	
	return chi_2
	
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

