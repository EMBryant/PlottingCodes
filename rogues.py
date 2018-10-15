'''Module containg various useful functions for Ed's exoplanet/TESS/NGTS PhD'''

from __future__ import division
import batman
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize

def phase_fold(time, epoch, period, max_phase):
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

def light_curve_parameters_solve(params_initial):
	'''Function to find the best fit parameters to describe a light curve, given an initial set of parameters:
	
			INPUTS:
				params_initial:			a batman.TransitParams instance containing initial guess parameter set
			
			OUTPUTS:
				params_best:			a batman.TransitParams instance containing the best fit parameter set
				chi_squared_best:		a chi squared value for the best fit parameters '''
				
	
def TIC_byID(ID):
	
	from astroquery.mast import Catalogs
	
	catTable = Catalogs.query_criteria(ID = ID, catalog="Tic")
	return catTable
		
	
	
	
	
	
