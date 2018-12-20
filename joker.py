#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  joker.py
#  
#
#  29th November 2018
#  Edward Bryant <phrvdf@monju.astro.warwick.ac.uk>
#  
#  This script is to take BRIGHTEN working group NGTS observations of TOIs and plot them  
#  overlaid with a batman model LC which has been obtained from a least squares fit to the TESS LC

import argparse
import numpy as np
import batman as bm
from scipy.stats import chisquare, sem
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit
from astropy.io import fits as pyfits
from lightkurve import TessLightCurve
import pandas

def ParseArgs():
    '''
    Function to parse command line arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--fn_tess', type=str, help="Name of TESS LC data file")
    parser.add_argument('--fn_ngts', type=str, nargs='*', help="Name of NGTS LC data file")
    parser.add_argument('--ngts_id', type=str, nargs='*', help="Action ID of the NGTS obs")
    parser.add_argument('--pipe', type=str, default='spoc', help="Pipeline used to reduce data (spoc or qlp)")
    parser.add_argument('-lk', '--lightkurve', action='store_true', help="Add this to use lightkurve to flatten the flux time series")
    parser.add_argument('--wl', type=int)
    parser.add_argument('--prim', action='store_true', help="Add this to indicate primary planet in multi-planet system")
    parser.add_argument('--sec', action='store_true', help="Add this to indicate secondary planet in multi-planet system")
    parser.add_argument('--tert', action='store_true', help="Add this to indicate tertiary planet in multi-planet system")
    return parser.parse_args()
    
def lc_min(params, phase, flux, err):
	'''
	Function which calculates Chi2 value for a given set of input parameters.
	Function to be minimized to find best fit parameters
	'''
	
	#Define the system parameters for the batman LC model
	pm = bm.TransitParams()
	
	pm.t0 = params['t0'].value   #Time of transit centre
#	pm.t0 = 0.   #Time of transit centre
	pm.per = 1.                  #Orbital period = 1 as phase folded
	pm.rp = params['rp'].value   #Ratio of planet to stellar radius
	pm.a = params['a'].value             #Semi-major axis (units of stellar radius)
	pm.inc = params['inc'] .value          #Orbital Inclination [deg]
	pm.ecc = 0.                  #Orbital eccentricity (fix circular orbits)
	pm.w = 90.                   #Longitude of periastron [deg] (unimportant as circular orbits)
	pm.u = [0.3, 0.1]            #Stellar LD coefficients
	pm.limb_dark="quadratic"     #LD model
	
	#Initialize the batman LC model and compute model LC
	m = bm.TransitModel(pm, phase)
	f_model = m.light_curve(pm)
	residuals = (flux - f_model)**2/err**2
	return residuals
	
def lc_t0_ngts_min(params, time, flux, err, pm_ngts):
    '''
    Function that will be mnimized to fit the transit centre time for the 
    NGTS obs, given the TESS Model LC
    '''
    
    pm_ngts.t0 = params['t0'].value   #Time pf transit centre
    m = bm.TransitModel(pm_ngts, time)
    f_model = m.light_curve(pm_ngts)
    residuals = (flux - f_model)**2/err**2
    return residuals
	
def lc_bin(time, flux, err, bin_width):
	'''
	Function to bin the data into bins of a given width. time and bin_width 
	must have the same units
	'''
	
	edges = np.arange(np.min(time), np.max(time), bin_width)
	dig = np.digitize(time, edges)
	time_binned = (edges[1:] + edges[:-1]) / 2
	flux_binned = np.array([np.nan if len(flux[dig == i]) == 0 else flux[dig == i].mean() for i in range(1, len(edges))])
#	err_binned = np.array([np.nan if len(flux[dig == i]) == 0 else sem(flux[dig == i]) for i in range(1, len(edges))])
	err_binned = np.array([np.nan if len(flux[dig == i]) == 0 else np.sqrt(np.sum(err[dig == i]**2))/len(err[dig == i]) for i in range(1, len(edges))])
	time_bin = time_binned[~np.isnan(err_binned)]
	err_bin = err_binned[~np.isnan(err_binned)]
	flux_bin = flux_binned[~np.isnan(err_binned)]	
	
	return time_bin, flux_bin, err_bin	

def tess_LC_dataload_spoc(file_name):
	'''Loads TESS LC data for a given object and uses quality flags etc. to remove bad points 
		'''
	#Load FITS file
	hdul = pyfits.open(file_name)               #Read in FITS file
	hdr = hdul[0].header                  		#Primary FITS header
	
	TIC = hdr['TICID']							#Loads TIC ID number from the fits header
#	Tmag = hdr['TESSMAG']						#Loads TESS mag from fits header
#	Rs = hdr['RADIUS']							#Loads stellar radius [solar radius] from fits header
#	sector = hdr['SECTOR']						#Loads observing sector data was obtained in
			
	
	DATA_WHOLE = hdul[1].data             		#Extracts whole data

	#Extract desired columns
	time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
	#time_corr = DATA_WHOLE['TIMECORR']    		#Time correction: time - time_corr gives light arrival time at spacecraft

	FLUX = DATA_WHOLE['PDCSAP_FLUX']  		#PDC corrected flux from target star
	FLUX_ERR = DATA_WHOLE['PDCSAP_FLUX_ERR']#Error in PDC corrected flux
	
	#Load Quality flags in to remove flagged data points
	flags = DATA_WHOLE['QUALITY']
	flag_indices = np.where(flags > 0)

	flux_flagremoved = np.delete(FLUX, flag_indices)
	fluxerr_flagremoved = np.delete(FLUX_ERR, flag_indices)
	time_flagremoved = np.delete(time, flag_indices)

	#Remove time points during central gap
	null_indices = np.where(np.isnan(time_flagremoved))
	time_nullremoved = np.delete(time_flagremoved, null_indices)
	flux_nullremoved = np.delete(flux_flagremoved, null_indices)
	fluxerr_nullremoved = np.delete(fluxerr_flagremoved, null_indices)
	
	return time_nullremoved, flux_nullremoved, fluxerr_nullremoved, TIC#, Tmag, Rs
	
def tess_LC_dataload_qlp(file_name):
	
	#Load FITS file
	hdul = pyfits.open(file_name)           		#Read in FITS file
	hdr = hdul[0].header                  		#Primary FITS header
	DATA_WHOLE = hdul[1].data             		#Extracts whole data
	
	TIC = np.int(hdr['TIC'])						#Loads TIC ID number from the fits header
	
	#Extract desired columns
	time = DATA_WHOLE['TIME']			  		#Time [BJD - 2457000]
	FLUX = DATA_WHOLE['SAP_FLUX']
	
	#Load Quality flags in to remove flagged data points
	flags = DATA_WHOLE['QUALITY']
	flag_indices = np.where(flags > 0)

	flux_flagremoved = np.delete(FLUX, flag_indices)
	time_flagremoved = np.delete(time, flag_indices)

	#Remove time points during central gap
	null_indices = np.where(np.isnan(time_flagremoved))
	time_nullremoved = np.delete(time_flagremoved, null_indices)
	flux_nullremoved = np.delete(flux_flagremoved, null_indices)
	
	return time_nullremoved, flux_nullremoved, TIC
	
def lc_detrend(time, flux, err, per, epoc, tdur):
    '''
    Function to flatten LCs. Removes transits and then fits a moving 
    average trend using a Savitzky-Golay filter - implemented using lightkurve
    '''
    
    flux_trend = np.zeros_like(flux) + flux
    phase_trend = (time - epoc) / per
    
    n_trans = np.int((time[-1] - epoc) / per + 1)
    for i in range(n_trans):
        trans = np.array((time[np.abs(phase_trend - i) < tdur/24/per], flux[np.abs(phase_trend - i) < tdur/24/per]))
        length = np.int(len(trans[0])/4)
		
        m = (np.mean(trans[1, -1*length:]) - np.mean(trans[1, :length])) / (np.mean(trans[0, -1*length:]) - np.mean(trans[0, :length]))
        c = np.mean(trans[1, :length]) - m*np.mean(trans[0, :length])
		
        flux_trend[np.where(np.abs(phase_trend - i) < tdur/24/per)] = m*trans[0] + c
	
    flat_lc, trend_lc = TessLightCurve(time, flux_trend).flatten(window_length=101, return_trend=True)
    
    flux_flat = flux/trend_lc.flux
    err_flat = err/trend_lc.flux
    return flux_flat, err_flat, trend_lc.flux, flux_trend

def ngts_LC_dataload(file_name):
	
    DATA = np.loadtxt(file_name)
    time, flux, err = DATA[:, 2], DATA[:, 3], DATA[:, 4]
    return time, flux, err
    

def lc_fit(time, flux, err, period, epoc, tdur, objid):
    '''
    Function to fit a batman model to an input lc datafile to find the best 
    fit system parameters
    '''
        
    phase = ((time - epoc)/period)%1  #Convert time values in to phases
    phase = np.array([p-1 if p > 0.8 else p for p in phase], dtype=float)
    
    p_fit = phase[phase < 0.2]  #Crop phase and flux arrays to only contain values
    f_fit = flux[phase < 0.2]   #in range (-0.2 ,  0.2)
    e_fit = err[phase < 0.2]
    
    transit_indices = np.where(np.abs(p_fit * period) <= tdur / (2 * 24))	#Array indices of all phase/flux values during the transit
    FLUX_OOT = np.delete(f_fit, transit_indices)				#"Out Of Transit" flux values

    median = np.median(FLUX_OOT)										#median of all out-of-transit flux values

    f_fit /= median
    e_fit /= median
    
    params=Parameters()         #Parameter instance to hold fit parameters
    params.add('rp', value=0.05, min=0., max=1.)    #Planet:Star radius ratio
    params.add('a', value=10., min=0., max=100.)    #Semi-major axis
    params.add('inc', value=89., min=60., max=90.)  #Orbital inclination
    params.add('t0', value=0.0, min=-0.005, max=0.005) #Transit centre time
    
    res = minimize(lc_min, params, args=(p_fit, f_fit, e_fit), method='leastsq') #perform minimization
    chi2 = np.sum(res.residual) / res.nfree
    t0_best, rp_best, a_best, inc_best = res.params['t0'].value, res.params['rp'].value, res.params['a'].value, res.params['inc'].value
#    rp_best, a_best, inc_best = res.params['rp'].value, res.params['a'].value, res.params['inc'].value
#    t0_best = 0.
    print('Best fit parameters: rp={:.8f}; a={:.8f}; inc={:.8f}; t0={:.8f}'.format(rp_best, a_best, inc_best, t0_best))
   
    print('Minimization result: {}: {}; chi2={:.4f}'.format(res.success, res.message, chi2))
    
    #Produce a best fit model using the minimization results
    pm_best = bm.TransitParams()
    
    pm_best.t0 = t0_best                #Time of transit centre
    
    pm_best.per = 1.               #Orbital period = 1 as phase folded
    pm_best.rp = rp_best           #Ratio of planet to stellar radius
    pm_best.a = a_best             #Semi-major axis (units of stellar radius)
    pm_best.inc = inc_best         #Orbital Inclination [deg]
    pm_best.ecc = 0.               #Orbital eccentricity (fix circular orbits)
    pm_best.w = 90.                #Longitude of periastron [deg] (unimportant as circular orbits)
    pm_best.u = [0.1, 0.3]         #Stellar LD coefficients
    pm_best.limb_dark="quadratic"  #LD model
    
    p_best = np.linspace(-0.2, 0.2, 10000)     #Produce a model LC using 
    m_best = bm.TransitModel(pm_best, p_best)  #the output best fit parameters
    f_best = m_best.light_curve(pm_best)
    
    if len(np.where(f_best < 1)[0]) > 0:
        p1 = p_best[np.where(f_best < 1)[0][0]]    #Phase of first contact
        p4 = p_best[np.where(f_best < 1)[0][-1]]   #Phase of final contact
        
        t_dur = (p4 - p1) * period *24             #Transit duration [hours]
        t_depth = (1 - f_best.min()) * 100         #Transit depth [percent]
        
    else:
        t_dur = np.nan
        t_depth = np.nan
		
    b_best = a_best * np.cos(inc_best * np.pi / 180)  #Impact Parameter
    
    return p_best, f_best, p_fit, f_fit, e_fit, rp_best, a_best, inc_best, t0_best, t_dur, t_depth, b_best, res, chi2
 
def lc_fit_lk(time, flux, err, period, epoc, tdur, objid):
    '''
    Function to fit a batman model to an input lc datafile to find the best 
    fit system parameters
    '''
    
    phase = ((time - epoc)/period)%1  #Convert time values in to phases
    phase = np.array([p-1 if p > 0.8 else p for p in phase], dtype=float)
    
    p_fit = phase[phase < 0.2]  #Crop phase and flux arrays to only contain values
    f_fit = flux[phase < 0.2]   #in range (-0.2 ,  0.2)
    e_fit = err[phase < 0.2]
    
    params=Parameters()         #Parameter instance to hold fit parameters
    params.add('rp', value=0.05, min=0., max=1.)    #Planet:Star radius ratio
    params.add('a', value=10., min=0., max=100.)    #Semi-major axis
    params.add('inc', value=89., min=60., max=90.)  #Orbital inclination
    params.add('t0', value=0.0, min=-0.005, max=0.005) #Transit centre time
    
    res = minimize(lc_min, params, args=(p_fit, f_fit, e_fit), method='leastsq') #perform minimization
    chi2 = np.sum(res.residual) / res.nfree
    t0_best, rp_best, a_best, inc_best = res.params['t0'].value, res.params['rp'].value, res.params['a'].value, res.params['inc'].value

    print('Best fit parameters: rp={:.8f}; a={:.8f}; inc={:.8f}; t0={:.8f}'.format(rp_best, a_best, inc_best, t0_best))
   
    print('Minimization result: {}: {}; chi2={:.4f}'.format(res.success, res.message, chi2))
    
    #Produce a best fit model using the minimization results
    pm_best = bm.TransitParams()
    
    pm_best.t0 = t0_best                #Time of transit centre
    
    pm_best.per = 1.               #Orbital period = 1 as phase folded
    pm_best.rp = rp_best           #Ratio of planet to stellar radius
    pm_best.a = a_best             #Semi-major axis (units of stellar radius)
    pm_best.inc = inc_best         #Orbital Inclination [deg]
    pm_best.ecc = 0.               #Orbital eccentricity (fix circular orbits)
    pm_best.w = 90.                #Longitude of periastron [deg] (unimportant as circular orbits)
    pm_best.u = [0.1, 0.3]         #Stellar LD coefficients
    pm_best.limb_dark="quadratic"  #LD model
    
    p_best = np.linspace(-0.2, 0.2, 10000)     #Produce a model LC using 
    m_best = bm.TransitModel(pm_best, p_best)  #the output best fit parameters
    f_best = m_best.light_curve(pm_best)
    
    if len(np.where(f_best < 1)[0]) > 0:
        p1 = p_best[np.where(f_best < 1)[0][0]]    #Phase of first contact
        p4 = p_best[np.where(f_best < 1)[0][-1]]   #Phase of final contact
        
        t_dur = (p4 - p1) * period *24             #Transit duration [hours]
        t_depth = (1 - f_best.min()) * 100         #Transit depth [percent]
        
    else:
        t_dur = np.nan
        t_depth = np.nan
		
    b_best = a_best * np.cos(inc_best * np.pi / 180)
    
    return p_best, f_best, p_fit, f_fit, e_fit, rp_best, a_best, inc_best, t0_best, t_dur, t_depth, b_best, res, chi2

def t0_fit_ngts(time, flux, err, pm_ngts):
    '''
    Function to fit the transit centre of the NGTS obs, given the TESS
    model LC
    '''
    
    params=Parameters()
    params.add('t0', value=pm_ngts.t0, min=pm_ngts.t0 - 0.005*pm_ngts.per, max=pm_ngts.t0 + 0.005*pm_ngts.per)
    res = minimize(lc_t0_ngts_min, params, args=(time, flux, err, pm_ngts), method='leastsq')
    chi2 = np.sum(res.residual) / res.nfree
    t0_best = res.params['t0'].value

    print('Best fit NGTS t0 variation = {:.8f}'.format(t0_best - pm_ngts.t0))
    print('Minimization result: {}: {}; chi2={:.4f}'.format(res.success, res.message, chi2))
    
    return t0_best
    
def rms(time, flux, model):
    
    err = np.ones(len(time))
    t30, f30, e30 = lc_bin(time, flux, err, 30/1440)
    t30_model, f30_model, e30_model = lc_bin(time, model, err, 30/1440)    
    rms30 = np.std(f30/f30_model) * 10**6     #rms of NGTS obs, averaged over 30 mins, in ppm
    
    return rms30
  
if __name__ == '__main__':
    
    args = ParseArgs()
    
    ###################################################################
    ############### TESS DATA #########################################
    ###################################################################
    
    df = pandas.read_csv('/home/astro/phrvdf/tess_data_alerts/TOIs_20181016.csv', index_col='tic_id')		#.csv file containing info on parameters (period, epoch, ID, etc.) of all TOIs
    length=len(df.iloc[0])
    
    if args.pipe == 'spoc':
		
        time_tess, flux_tess, err_tess, tic = tess_LC_dataload_spoc(args.fn_tess)
        
        df2 = df.loc[tic]
        
        if args.prim:
            df2 = df2.iloc[0]
		
        elif args.sec:
            df2 = df2.iloc[1]
			
        elif args.tert:
            df2 = df2.iloc[2]
        
        epoc = df2.loc['Epoc']
        per = df2.loc['Period']
        tdur = df2.loc['Duration']
        toi = df2.loc['toi_id']
        
        print('Running fit for TIC {} (TOI-00{})'.format(tic, toi))
        
        print('Epoc: {}; Per: {}; tdur: {}'.format(epoc, per, tdur))
        
        
        if args.lightkurve == True:
            
            flux_flat, err_flat, trend_flux, flux_trend = lc_detrend(time_tess, flux_tess, err_tess, per, epoc, tdur)
            
            p_best, f_best, p_fit, f_fit, e_fit, rp_best, a_best, inc_best, t0_best, t_dur, t_depth, b_best, res, chi2 = lc_fit_lk(time_tess, flux_flat, err_flat, per, epoc, tdur, tic)
    
        else: p_best, f_best, p_fit, f_fit, e_fit, rp_best, a_best, inc_best, t0_best, t_dur, t_depth, b_best, res, chi2 = lc_fit(time_tess, flux_tess, err_tess, per, epoc, tdur, tic)
    
        #Produce binned data set for plotting
        bw = 5 / (1440*per)                    #Bin width - 10 mins in units of phase
        p_bin_tess, f_bin_tess, e_bin_tess = lc_bin(p_fit, f_fit, e_fit, bw)
    
    elif args.pipe == 'qlp':
		
        time_tess, flux_tess, tic = tess_LC_dataload_qlp(args.fn_tess)
		
        epoc = df.loc[tic, 'Epoc']
        per = df.loc[tic, 'Period']
        tdur = df.loc[tic, 'Duration']
        toi = df.loc[tic, 'toi_id']
        
        print('Running fit for TIC {} (TOI-00{})'.format(tic, toi))
        
        print('Epoc: {}; Per: {}; tdur: {}'.format(epoc, per, tdur))
        
        p_best, f_best, p_fit, f_fit, e_fit, rp_best, a_best, inc_best, t0_best, t_dur, t_depth, b_best, res, chi2 = lc_fit(time_tess, flux_tess, np.zeros(len(flux_tess)) + 0.001, per, epoc, tdur, tic)
    
    ###################################################################
    ############### NGTS DATA #########################################
    ###################################################################
    
    #Produce a batman model for the NGTS data
    
    pm_ngts = bm.TransitParams()
    
    pm_ngts.t0 = epoc + t0_best*per
    pm_ngts.per = per
    pm_ngts.rp = rp_best
    pm_ngts.a = a_best
    pm_ngts.inc = inc_best
    pm_ngts.ecc = 0.
    pm_ngts.w = 90.
    pm_ngts.u = [0.1, 0.3]
    pm_ngts.limb_dark = "quadratic"
    
    
    time_ngts1, flux_ngts1, err_ngts1 = ngts_LC_dataload(args.fn_ngts[0])
    tbin_ngts1, fbin_ngts1, ebin_ngts1 = lc_bin(time_ngts1, flux_ngts1, err_ngts1, 5/1440)
    
    t_offset1 = np.int(time_ngts1[0]) 
    
    time_model1 = time_ngts1 - 2457000
    
    m_ngts1 = bm.TransitModel(pm_ngts, time_model1)
    flux_model_ngts1 = m_ngts1.light_curve(pm_ngts)
    
    rms30_ngts1 = rms(time_ngts1, flux_ngts1, flux_model_ngts1)
    
    if len(args.fn_ngts) == 2:
        time_ngts2, flux_ngts2, err_ngts2 = ngts_LC_dataload(args.fn_ngts[1])
        tbin_ngts2, fbin_ngts2, ebin_ngts2 = lc_bin(time_ngts2, flux_ngts2, err_ngts2, 5/1440)
        
        time_ngts = np.append(time_ngts1, time_ngts2)
        flux_ngts = np.append(flux_ngts1, flux_ngts2)
        err_ngts = np.append(err_ngts1, err_ngts2)
        
        t_offset2 = np.int(time_ngts2[0])
        time_model2 = time_ngts2 - 2457000
        
        m_ngts2 = bm.TransitModel(pm_ngts, time_model2)
        flux_model_ngts2 = m_ngts2.light_curve(pm_ngts)
        
        rms30_ngts2 = rms(time_ngts2, flux_ngts2, flux_model_ngts2)
        
        phase_ngts = ((time_ngts - 2457000 - (epoc+t0_best*per))/per)%1
        for j in range(len(phase_ngts)):
            if phase_ngts[j] > 0.5: phase_ngts[j] -= 1.
        pbin_ngts, fbin_ngts, ebin_ngts = lc_bin(phase_ngts, flux_ngts, err_ngts, 5/(1440*per))
		
    elif len(args.fn_ngts) == 3:
        time_ngts2, flux_ngts2, err_ngts2 = ngts_LC_dataload(args.fn_ngts[1])
        tbin_ngts2, fbin_ngts2, ebin_ngts2 = lc_bin(time_ngts2, flux_ngts2, err_ngts2, 5/1440)
        
        t_offset2 = np.int(time_ngts2[0])
        time_model2 = time_ngts2 - 2457000
        
        m_ngts2 = bm.TransitModel(pm_ngts, time_model2)
        flux_model_ngts2 = m_ngts2.light_curve(pm_ngts)
        
        rms30_ngts2 = rms(time_ngts2, flux_ngts2, flux_model_ngts2)
        
        time_ngts3, flux_ngts3, err_ngts3 = ngts_LC_dataload(args.fn_ngts[2])
        tbin_ngts3, fbin_ngts3, ebin_ngts3 = lc_bin(time_ngts3, flux_ngts3, err_ngts3, 5/1440)
		
        t_offset3 = np.int(time_ngts3[0])
        time_model3 = time_ngts3 - 2457000
        
        m_ngts3 = bm.TransitModel(pm_ngts, time_model3)
        flux_model_ngts3 = m_ngts3.light_curve(pm_ngts)
        
        rms30_ngts3 = rms(time_ngts3, flux_ngts3, flux_model_ngts3)
        
        time_ngts = np.append(time_ngts1, np.append(time_ngts2, time_ngts3))
        flux_ngts = np.append(flux_ngts1, np.append(flux_ngts2, flux_ngts3))
        err_ngts = np.append(err_ngts1, np.append(err_ngts2, err_ngts3))      
        
        phase_ngts = ((time_ngts - 2457000 - (epoc+t0_best*per))/per)%1
        for j in range(len(phase_ngts)):
            if phase_ngts[j] > 0.5: phase_ngts[j] -= 1.
        pbin_ngts, fbin_ngts, ebin_ngts = lc_bin(phase_ngts, flux_ngts, err_ngts, 5/(1440*per))
		
#    time_model = time_ngts - 2457000
    
    
#    m_ngts = bm.TransitModel(pm_ngts, time_model)
#    flux_model_ngts = m_ngts.light_curve(pm_ngts)
    
#    pm_ngts.t0 = t0_fit_ngts(time_ngts, flux_ngts, err_ngts, pm_ngts)
#    flux_model_ngts2 = m_ngts.light_curve(pm_ngts)
    
    
    #################Calculate RMS averaged over 30 mins################
    
#    t30_ngts, f30_ngts, e30_ngts = lc_bin(time_ngts, flux_ngts, err_ngts, 30/1440)
#    t30_model, f30_model, e30_model = lc_bin(time_model1, flux_model_ngts1, err_ngts, 30/1440)
    
#    rms30_ngts = np.std(f30_ngts/f30_model) * 10**6     #rms of NGTS obs, averaged over 30 mins, in ppm
    
    if args.pipe == 'spoc':
        rms30_tess = pyfits.open(args.fn_tess)[1].header['CDPP0_5']
    
    elif args.pipe == 'qlp':
        pm_ngts.t0 = epoc
        m_tess = bm.TransitModel(pm_ngts, np.array(time_tess, dtype=float))
        f_model_tess = m_tess.light_curve(pm_ngts)
        
        rms30_tess = np.std(flux_tess/f_model_tess) * 10**6
    

    ####################################################################
    ######################### PLOTTING #################################
    ####################################################################

    axis_font = {'fontname':'DejaVu Sans', 'size':'15'}
    
    #Produce plot of data and best fit model LC
    fig = plt.figure(figsize=(9, 7.5))
    
    #Plot the TESS phase folded LC
    if len(args.fn_ngts) == 1:
        ax1 = fig.add_subplot(211)
        ax1.set_xlabel('Phase', **axis_font)

    else:
        ax1 = fig.add_subplot(311)
    if args.pipe == 'spoc': 
        ax1.plot(p_fit, f_fit, marker='o', color='gray', linestyle='none', markersize=2)
        ax1.errorbar(p_bin_tess, f_bin_tess, yerr=e_bin_tess, marker='o', linestyle='none', color='red', markersize=5)
    elif args.pipe == 'qlp':
        ax1.plot(p_fit, f_fit, 'ko', markersize=5)
    ax1.plot(p_best, f_best, 'g--', linewidth=2)
    
    ax1.set_ylabel('Flux - TESS LC', **axis_font)
    ax1.set_title('TOI-00{};   TIC  {};  RMS 30 min TESS: {:.1f} ppm \n Depth: {:.5f}%;  Duration: {:5f} hours; Epoc; {:.6f} + {:.8f} \n (Rp/Rs): {:.6f};  b: {:.4f}; chi2: {:.8f}; Period: {:.6f} days'.format(toi, tic, rms30_tess, t_depth, t_dur, epoc, t0_best, rp_best, b_best, chi2, per), **axis_font)


    #Plot the NGTS DATA
    
    #Only one NGTS night
    if len(args.fn_ngts) == 1:
		    
        ax2 = fig.add_subplot(212)
    
        ax2.plot(time_ngts1 - t_offset1, flux_ngts1, marker='o', color='gray', linestyle='none', markersize=2)
        ax2.errorbar(tbin_ngts1 - t_offset1, fbin_ngts1, yerr=ebin_ngts1, marker='o', linestyle='none', color='red', markersize=5)
        ax2.plot(time_ngts1 - t_offset1, flux_model_ngts1, 'g--', linewidth=2)
    
        ax2.set_xlabel('BJD_TDB - {} +'.format(t_offset1), **axis_font)
        ax2.set_ylabel('Flux - NGTS LC', **axis_font)
        ax2.set_title('NGTS Action ID: {};  RMS 30 min NGTS: {:.1f} ppm'.format(args.ngts_id, rms30_ngts1))
    
    #Two NGTS nights
    elif len(args.fn_ngts) == 2:
		
        ax2 = fig.add_subplot(312)
		
        ax2.plot(phase_ngts, flux_ngts, marker='o', color='gray', linestyle='none', markersize=2)
        ax2.errorbar(pbin_ngts, fbin_ngts, yerr=ebin_ngts, marker='o', linestyle='none', color='red', markersize=5)
        ax2.plot(p_best, f_best, 'g--', linewidth=2)
        ax2.set_xlabel('Phase', **axis_font)
        ax2.set_ylabel('FLUX - NGTS LC', **axis_font)

		
        ax3 = fig.add_subplot(325)
    
        ax3.plot(time_ngts1 - t_offset1, flux_ngts1, marker='o', color='gray', linestyle='none', markersize=2)
        ax3.errorbar(tbin_ngts1 - t_offset1, fbin_ngts1, yerr=ebin_ngts1, marker='o', linestyle='none', color='red', markersize=5)
        ax3.plot(time_ngts1 - t_offset1, flux_model_ngts1, 'g--', linewidth=2)
    
        ax3.set_xlabel('BJD_TDB - {} +'.format(t_offset1), **axis_font)
        ax3.set_ylabel('Flux - NGTS LC', **axis_font)
        ax3.set_title('NGTS Action ID: {};  RMS 30 min NGTS: {:.1f} ppm'.format(args.ngts_id[0], rms30_ngts1))
    
		
        ax4 = fig.add_subplot(326)
    
        ax4.plot(time_ngts2 - t_offset2, flux_ngts2, marker='o', color='gray', linestyle='none', markersize=2)
        ax4.errorbar(tbin_ngts2 - t_offset2, fbin_ngts2, yerr=ebin_ngts2, marker='o', linestyle='none', color='red', markersize=5)
        ax4.plot(time_ngts2 - t_offset2, flux_model_ngts2, 'g--', linewidth=2)
    
        ax4.set_xlabel('BJD_TDB - {} +'.format(t_offset2), **axis_font)
        ax4.set_title('NGTS Action ID: {};  RMS 30 min NGTS: {:.1f} ppm'.format(args.ngts_id[1], rms30_ngts2))
        
    #Three NGTS nights
    elif len(args.fn_ngts) == 3:
		
        ax2 = fig.add_subplot(312)
		
        ax2.plot(phase_ngts, flux_ngts, marker='o', color='gray', linestyle='none', markersize=2)
        ax2.errorbar(pbin_ngts, fbin_ngts, yerr=ebin_ngts, marker='o', linestyle='none', color='red', markersize=5)
        ax2.plot(p_best, f_best, 'g--', linewidth=2)
        ax2.set_xlabel('Phase', **axis_font)
        ax2.set_ylabel('FLUX - NGTS LC', **axis_font)

		
        ax3 = fig.add_subplot(337)
    
        ax3.plot(time_ngts1 - t_offset1, flux_ngts1, marker='o', color='gray', linestyle='none', markersize=2)
        ax3.errorbar(tbin_ngts1 - t_offset1, fbin_ngts1, yerr=ebin_ngts1, marker='o', linestyle='none', color='red', markersize=5)
        ax3.plot(time_ngts1 - t_offset1, flux_model_ngts1, 'g--', linewidth=2)
    
        ax3.set_xlabel('BJD_TDB - {} +'.format(t_offset1), **axis_font)
        ax3.set_ylabel('Flux - NGTS LC', **axis_font)
        ax3.set_title('NGTS Action ID: {};  RMS 30 min NGTS: {:.1f} ppm'.format(args.ngts_id[0], rms30_ngts1))
    
		
        ax4 = fig.add_subplot(338)
    
        ax4.plot(time_ngts2 - t_offset2, flux_ngts2, marker='o', color='gray', linestyle='none', markersize=2)
        ax4.errorbar(tbin_ngts2 - t_offset2, fbin_ngts2, yerr=ebin_ngts2, marker='o', linestyle='none', color='red', markersize=5)
        ax4.plot(time_ngts2 - t_offset2, flux_model_ngts2, 'g--', linewidth=2)
    
        ax4.set_xlabel('BJD_TDB - {} +'.format(t_offset2), **axis_font)
        ax4.set_title('NGTS Action ID: {};  RMS 30 min NGTS: {:.1f} ppm'.format(args.ngts_id[1], rms30_ngts2))
        
        ax5 = fig.add_subplot(339)
    
        ax5.plot(time_ngts3 - t_offset3, flux_ngts3, marker='o', color='gray', linestyle='none', markersize=2)
        ax5.errorbar(tbin_ngts3 - t_offset3, fbin_ngts3, yerr=ebin_ngts3, marker='o', linestyle='none', color='red', markersize=5)
        ax5.plot(time_ngts3 - t_offset3, flux_model_ngts3, 'g--', linewidth=2)
    
        ax5.set_xlabel('BJD_TDB - {} +'.format(t_offset3), **axis_font)
        ax5.set_title('NGTS Action ID: {};  RMS 30 min NGTS: {:.1f} ppm'.format(args.ngts_id[2], rms30_ngts3))
      
    plt.show()
