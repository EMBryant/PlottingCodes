'''Module containg various useful functions for Ed's exoplanet/TESS/NGTS PhD'''

from __future__ import division
import batman
import numpy as np
from matplotlib import pyplot as plt

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

	params = batman.TransitParams()
	params.t0 = t0
	params.per = per
	params.rp = rp
	params.a = a
	params.inc = inc
	params.ecc = ecc
	params.w = w
	params.u = u
	params.limb_dark = limb_dark
	
	m = batman.TransitModel(params, t)
	flux_model = m.light_curve(params)
	
	return flux_model
