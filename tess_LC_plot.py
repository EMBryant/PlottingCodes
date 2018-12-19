'''Module containg various useful functions for Ed's exoplanet/TESS/NGTS PhD'''
import batman as bm
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize as mini
from scipy.stats import sem
import pandas
import argparse
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
from lightkurve import TessLightCurve
	
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
	
	#Remove remaining nans
	fluxerr_nullremoved = fluxerr_nullremoved[~np.isnan(flux_nullremoved)]
	time_nullremoved = time_nullremoved[~np.isnan(flux_nullremoved)]
	raw_flux_nullremoved = raw_flux_nullremoved[~np.isnan(flux_nullremoved)]
	flux_nullremoved = flux_nullremoved[~np.isnan(flux_nullremoved)]
	
	
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
	
def normalise_LC(flux, err, phase, period, Tdur):

	#Now we want to mask out the transit points, to do statistices on the rest of the LC
	transit_indices = np.where(np.abs(phase * period) <= T_dur / (2 * 24))	#Array indices of all phase/flux values during the transit
	FLUX_OOT = np.delete(flux, transit_indices)				#"Out Of Transit" flux values

	median = np.median(FLUX_OOT)										#median of all out-of-transit flux values

	flux_normalised = flux / median
	err_normalised = err / median
	
	return flux_normalised, err_normalised
	
def lc_model(t, rp, a, t0 = 0., per = 1., inc = 90., ecc = 0., w = 90., limb_dark = "quadratic", u = [0.1, 0.3]):
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

	params = bm.TransitParams()						#a batman instance containing all the input parameters
	params.t0 = t0									#time of the transit centre
	params.per = per								#orbital period
	params.rp = rp									#radius of the planet, in units of stellar radii
	params.a = a									#semi-major axis of the planet's orbit, in units of stellar radii
	params.inc = inc								#orbital inclination [degrees]
	params.ecc = ecc								#orbital eccentricity
	params.w = w									#longitude of periastron [degrees]
	params.u = u									#limb darkening coefficients
	params.limb_dark = limb_dark					#limb darkening mechanism 
	
	m = bm.TransitModel(params, t)					#Initialises the transit model
	flux_model = m.light_curve(params)				#Calculates the light curve flux
	
	return flux_model
	
def fit_phasefolded(X0, phase, flux, err):
	
#	t0 = X0[0]
	rp = X0[0]
	a = X0[1]
	inc = X0[2]
	ecc = X0[3]
	w = X0[4]
	u1 = X0[5]
	u2 = X0[6]
	
	flux_model = lc_model(phase, rp, a, t0 = 0., per = 1., inc = inc, ecc = ecc, w = w, limb_dark = "quadratic", u = [u1, u2])

	chi_vals = np.sqrt((flux - flux_model)**2 / err**2)
	
	fit_val = np.sum(chi_vals) / (len(chi_vals) - 1)
	
	return fit_val

def best_fit_LC_solve(phase, flux, per, rp, a, t0=0., inc=89., ecc=0., w=90., ld="quadratic", u=[0.1, 0.3]):
	
	bins = [p * per * 24 * 60 // 10 for p in phase]
	length = np.int(max(bins) - min(bins))
	print(length)
	p_bin = np.zeros(length)
	f_bin = np.zeros(length)
	e_bin = np.zeros(length)
	sigma = np.std(flux)
	for i in range(length):
#		if len(np.where(bins == i + min(bins))[0]) == 1:
#			p_bin[i] = phase[i]
#			f_bin[i] = flux[i]
#			e_bin[i] = sigma
		
		if len(np.where(bins == i + min(bins))[0]) > 0:
			p_bin[i] = np.mean(phase[np.where(bins == i + min(bins))[0]])
			f_bin[i] = np.mean(flux[np.where(bins == i + min(bins))[0]])
			e_bin[i] = sem(flux[np.where(bins == i + min(bins))[0]])
		else:
			p_bin[i] = -1000
			f_bin[i] = -1000
			e_bin[i] = -1000
	
	p_bin = np.delete(p_bin, np.where(p_bin < -900))
	f_bin = np.delete(f_bin, np.where(f_bin < -900))
	e_bin = np.delete(e_bin, np.where(e_bin < -900))
	
	print(np.where(np.isnan(e_bin)))
	
	res = mini(fit_phasefolded, [rp, a, inc, ecc, w, u[0], u[1]], args=(p_bin, f_bin, e_bin))
	for i in range(20):
		if res.success == False:
			res = mini(fit_phasefolded, res.x, args=(p_bin, f_bin, e_bin))
		else: res = res
		
#		print(res.success)
#	t0_best = res.x[0]
	rp_best = res.x[0]
	a_best = res.x[1]
	inc_best = res.x[2]
	ecc_best = res.x[3]
	w_best = res.x[4]
	u1_best, u2_best = res.x[5], res.x[6]
	
	flux_best = lc_model(p_bin, rp_best, a_best, t0 = 0., per = 1., inc = inc_best, ecc = ecc_best, w = w_best, limb_dark = "quadratic", u = [u1_best, u2_best])
	
	return flux_best, p_bin, f_bin, e_bin, res.x, res.fun


def best_fit_LC_solve_qlp(phase, flux, rp, a, t0=0., inc=90., ecc=0., w=90., ld="quadratic", u=[0.1, 0.3]):
	
	phase = np.array(phase, dtype='float64')
	
	res = mini(fit_phasefolded, [rp, a, inc, ecc, w, u[0], u[1]], args=(phase, flux, np.zeros(len(flux))+np.std(flux)))
	
	for i in range(20):
		if res.success == False:
			res = mini(fit_phasefolded, res.x, args=(phase, flux, np.zeros(len(flux))+np.std(flux)))
		else: res = res
		
		#print(res.success)
#	t0_best = res.x[0]
	rp_best = res.x[0]
	a_best = res.x[1]
	inc_best = res.x[2]
	ecc_best = res.x[3]
	w_best = res.x[4]
	u1_best, u2_best = res.x[5], res.x[6]
	
	flux_best = lc_model(phase, rp_best, a_best, t0 = 0., per = 1., inc = inc_best, ecc = ecc_best, w = w_best, limb_dark = "quadratic", u = [u1_best, u2_best])
	
	return flux_best, res.x, res.fun


	
def plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, sector, epoc):
	
	n_trans = np.int((time[-1] - epoc) / period) + 1
	trans_markers = np.zeros(n_trans)
	for i in range(len(trans_markers)):
		trans_markers[i] = epoc + i*period
	
	fig = plt.figure(figsize=(30, 20))
	
	#Top Subplot: Raw SAP Flux
	ax1 = fig.add_subplot(311)
	ax1.plot(time, raw_flux, 'ko', markersize=1.5)
	ax1.plot(trans_markers, np.ones(n_trans)*np.min(raw_flux)-(np.median(raw_flux)-np.min(raw_flux)), 'r^', markersize=25)
	ax1.set_ylabel('Raw SAP Flux [e$^-$ / s]', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font) 
	ax1.tick_params(labelsize=20)
	ax1.set_title('T = {:.3f} ; Rs = {} \n Period: {:.4f} days; Epoc: {:.4f} \n TOI: {} ;  TIC ID: {} ;  Sector: {}'.format(Tmag, Rs, period, epoc, TOI, TIC, sector), **axis_font)


	#Middle subplot: Unfolded LC
	ax2 = fig.add_subplot(312)
	ax2.plot(time, flux_normalised, 'ko', markersize='1.5')
	ax2.plot(trans_markers, np.ones(n_trans)*np.min(flux_normalised)-(1 - np.min(flux_normalised)), 'r^', markersize=25)
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax2.tick_params(labelsize=20)
	
	#Bottom subplot: Folded LC
	ax3 = fig.add_subplot(313)

	ax3.plot(phase, flux_normalised, 'bo', markersize='1.5')
	ax3.set_ylabel('Relative Flux', **axis_font)
	ax3.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax3.tick_params(labelsize=20)
	ax3.set_xlabel('Phase', **axis_font)
	
	plt.tight_layout()
	
	if save:
		if lightkurve == True:
			plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}_lk.png'.format(TOI, TIC, sector))
		else:	
			plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}.png'.format(TOI, TIC, sector))
		plt.close()
	
	else:
		plt.show()
	
def plot_LC_model(flux_normalised, flux_model, phase, phase_model, f_bin, e_bin, TOI, TIC, period, Tmag, Rs, sector, rp, a, inc, fit_val):
	
	fig = plt.figure(figsize=(30, 20))
	
	#Bottom subplot: Folded LC
	ax1 = fig.add_subplot(111)

	ax1.plot(phase, flux_normalised, 'ko', markersize='1.5')
	ax1.errorbar(phase_model, f_bin, yerr=e_bin, marker='o', color='red', linestyle='none', markersize=3)
	ax1.plot(phase_model, flux_model, 'b-')
	ax1.set_ylabel('Relative Flux', **axis_font)
	ax1.set_xticks([-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
	ax1.tick_params(labelsize=20)
	ax1.set_xlabel('Phase', **axis_font)
	ax1.set_title('T = {:.2f}; Rs = {:.2f}; P = {:.2f} days; TIC: {}; Sector: {} \n Rp = {:.4f} Rs; a = {:.3f} Rs;  i = {:.3f} deg;  $\chi^2$ = {:.3f}'.format(Tmag, Rs, period, TIC, sector, rp, a, inc, fit_val), **axis_font)
	
#	plt.tight_layout()
	
	if save:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}_withmodel.png'.format(TOI, TIC, sector))
		plt.close()
	
	else:
		plt.show()

def plot_LC_model_qlp(flux_normalised, flux_model, phase, phase_model, TOI, TIC, period, Tmag, Rs, sector, rp, a, inc, fit_val):
	
	fig = plt.figure(figsize=(30, 20))
	
	#Bottom subplot: Folded LC
	ax1 = fig.add_subplot(111)

	ax1.plot(phase, flux_normalised, 'ko', markersize='3')
	ax1.plot(phase_model, flux_model, 'bo')
	ax1.set_ylabel('Relative Flux', **axis_font)
	ax1.set_xticks([-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
	ax1.tick_params(labelsize=20)
	ax1.set_xlabel('Phase', **axis_font)
	ax1.set_title('T = {:.2f}; Rs = {:.2f}; P = {:.2f} days; TIC: {}; Sector: {} \n Rp = {:.4f} Rs; a = {:.3f} Rs;  i = {:.3f} deg;  $\chi^2$ = {:.3f}'.format(Tmag, Rs, period, TIC, sector, rp, a, inc, fit_val), **axis_font)
	
#	plt.tight_layout()
	
	if save:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}_withmodel.png'.format(TOI, TIC, sector))
		plt.close()
	
	else:
		plt.show()


def plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs, sector):
	
	fig = plt.figure(figsize=(30, 20))
				
	#Top subplot: Unfolded LC
	ax1 = fig.add_subplot(211)
	ax1.plot(time, flux, 'ko', markersize='3')
	ax1.set_ylabel('Relative Flux', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax1.tick_params(labelsize=20)
	ax1.set_title('T = {} ; Rs = {} ;  Period: {} days \n TOI: {} ;  TIC ID: {} ; Sector: {}'.format(Tmag, Rs, period, TOI, TIC, sector), **axis_font)

	#Bottom subplot: Folded LC
	ax2 = fig.add_subplot(212)
	ax2.plot(phase, flux, 'bo', markersize='3')
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.tick_params(labelsize=20)
	ax2.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax2.set_xlabel('Phase', **axis_font)
				
	plt.tight_layout()
	
	if save:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}.png'.format(TOI, TIC, sector))
	
	else:
		plt.show()


def plot_LC_qlp_model(time, phase, flux, flux_model, phase_model, TOI, TIC, period, Tmag, Rs, sector):
	
	fig = plt.figure(figsize=(30, 20))
				
	#Top subplot: Unfolded LC
	ax1 = fig.add_subplot(211)
	ax1.plot(time, flux, 'ko', markersize='3')
	ax1.set_ylabel('Relative Flux', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax1.tick_params(labelsize=20)
	ax1.set_title('T = {} ; Rs = {} ;  Period: {} days \n TOI: {} ;  TIC ID: {} ; Sector: {}'.format(Tmag, Rs, period, TOI, TIC, sector), **axis_font)

	#Bottom subplot: Folded LC
	ax2 = fig.add_subplot(212)
	ax2.plot(phase, flux, 'bo', markersize='3')
	ax2.plot(phase_model, flux_model, 'r-')
	ax2.set_ylabel('Relative Flux', **axis_font)
	ax2.tick_params(labelsize=20)
	ax2.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
	ax2.set_xlabel('Phase', **axis_font)
				
	plt.tight_layout()
	
	if save:
		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}_withmodel.png'.format(TOI, TIC, sector))
	
	else:
		plt.show()

if __name__ == "__main__":
		
	parser = argparse.ArgumentParser()
	parser.add_argument('-fn', '--filename', type=str, nargs='*')
	parser.add_argument('-s', '--save', action='store_true')
	parser.add_argument('-lk', '--lightkurve', action='store_true')
	parser.add_argument('-wl', '--window_length', type=int)
	parser.add_argument('-m', '--model', action='store_true')
	parser.add_argument('-r', '--radius', type=float)
	parser.add_argument('-a', '--smaxis', type=float)
	
	args = parser.parse_args()
	
	filenames = args.filename
	save = args.save
	lightkurve = args.lightkurve
	windowlength = args.window_length
	mod = args.model
	rp, a = args.radius, args.smaxis
	
	axis_font = {'fontname':'DejaVu Sans', 'size':'20'}
	df = pandas.read_csv('/home/astro/phrvdf/tess_data_alerts/toi_ephems_20181212.csv', index_col='tic_id')		#.csv file containing info on parameters (period, epoch, ID, etc.) of all TOIs
	length = len(df.iloc[0])
		
	for i in range(len(filenames)):
		
		hdul = fits.open(filenames[i])
		hdr = hdul[0].header
		
		if hdr['ORIGIN'] == 'NASA/Ames':
		
			TIC = hdr['TICID']							#Loads TIC ID number from the fits header
			Tmag = hdr['TESSMAG']						#Loads TESS mag from fits header
			Rs = hdr['RADIUS']							#Loads stellar radius [solar radius] from fits header
			sector = hdr['SECTOR']						#Loads observing sector data was obtained in
			
			df2 = df.loc[TIC]
			
			if len(df2) == length:
				
				epoch = df2.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
				period = df2.loc['Period']			  	#Orbital Period [days]
				T_dur = df2.loc['Duration']				#Transit duration [hours]
				TOI = df2.loc['toi_id']            	  	#TIC ID for the object - used for plot title
#				tdepth = df2.loc['Transit Depth']
#				comments = df2.loc['Comment']			#Any existing comments on the object
			
				print("Epoch of first transit is {} [BJD - 2457000]".format(epoch))
				print("Orbital period is {} days".format(period))
				print("Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.))
#				print("Existing comments on this object are: {}".format(comments))
						
				time, flux, fluxerr, time_whole, raw_flux = tess_LC_dataload_spoc(filenames[i])
				
				if lightkurve == True:
					lc = TessLightCurve(time, flux)
					flat_lc = lc.flatten(window_length = windowlength)
					
					time = flat_lc.time
					flux = flat_lc.flux
					
				phase, phase_days = phase_fold(time, epoch, period)
			
				flux_normalised, err_normalised = normalise_LC(flux, fluxerr, phase, period, T_dur)
				
				if mod:
					flux_best, p_bin, f_bin, e_bin, best_fit_params, fit_val = best_fit_LC_solve(phase, flux_normalised, period, rp, a, t0=0., inc=89., ecc=0., w=90., ld="quadratic", u=[0.1, 0.3])
					
					print(best_fit_params, fit_val)
					
					plot_LC_model(flux_normalised, flux_best, phase, p_bin, f_bin, e_bin, TOI, TIC, period, Tmag, Rs, sector, best_fit_params[0], best_fit_params[1], best_fit_params[2], fit_val)
				else:
					plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, sector, epoch)

				print(TIC, TOI)
				
			else:
				for j in range(len(df2)):
					
					df3 = df2.iloc[j]
					
					epoch = df3.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
					period = df3.loc['Period']			  	#Orbital Period [days]
					T_dur = df3.loc['Duration']				#Transit duration [hours]
					TOI = df3.loc['toi_id']            	  	#TIC ID for the object - used for plot title
#					comments = df3.loc['Comment']			#Any existing comments on the object
			
					print("Epoch of first transit is {} [BJD - 2457000]".format(epoch))
					print("Orbital period is {} days".format(period))
					print("Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.))
#					print("Existing comments on this object are: {}".format(comments))
					
					time, flux, fluxerr, time_whole, raw_flux = tess_LC_dataload_spoc(filenames[i])
			
					phase, phase_days = phase_fold(time, epoch, period)
			
					flux_normalised, err_normalised = normalise_LC(flux, fluxerr, phase, period, T_dur)
					
					if mod:
						phase_sort = np.sort(phase)
						phase_model = np.array(phase_sort, dtype='float64')
						flux_model = lc_model()
						plot_LC_spoc_model(time, raw_flux, flux_normalised, flux_model, phase, phase_model, TOI, TIC, period, Tmag, Rs)
					else:
						plot_LC_spoc(time, raw_flux, flux_normalised, phase, TOI, TIC, period, Tmag, Rs, sector, epoch)

					print(TIC, TOI)
			
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
			
				print("Epoch of first transit is {} [BJD - 2457000]".format(epoch))
				print("Orbital period is {} days".format(period))
				print("Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.))
				print("Existing comments on this object are: {}".format(comments))
				

				time, flux = tess_LC_dataload_qlp(filenames[i])
			
				phase, phase_days = phase_fold(time, epoch, period)
								
				if mod:
					flux_best, best_fit_params, fit_val = best_fit_LC_solve_qlp(phase, flux, rp, a, t0=0., inc=89., ecc=0., w=90., ld="quadratic", u=[0.1, 0.3])
					
					s = phase.argsort()
					flux_model = flux_best[s]
					phase_model = phase[s]
					
					print(best_fit_params, fit_val)
					
					plot_LC_model_qlp(flux, flux_model, phase, phase_model, TOI, TIC, period, Tmag, Rs, 1, best_fit_params[0], best_fit_params[1], best_fit_params[2], fit_val)
				else:
					plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs)
				
				print(TIC, TOI)
		
			else:
				for j in range(len(df2)):
					
					df3 = df2.iloc[j]
					
					epoch = df3.loc['Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
					period = df3.loc['Period']			  	#Orbital Period [days]
					T_dur = df3.loc['Duration']				#Transit duration [hours]
					TOI = df3.loc['toi_id']            	  	#TIC ID for the object - used for plot title
					comments = df3.loc['Comment']			#Any existing comments on the object
			
					print("Epoch of first transit is {} [BJD - 2457000]".format(epoch))
					print("Orbital period is {} days".format(period))
					print("Transit duration is {} hours ({} days)".format(T_dur, T_dur/24.))
					print("Existing comments on this object are: {}".format(comments))
					
					time, flux, time_whole, raw_flux = tess_LC_dataload_qlp(filenames[i])
			
					phase, phase_days = phase_fold(time, epoch, period)
			
					if mod:
						phase_sort = np.sort(phase)
						phase_model = np.array(phase_sort, dtype='float64')
						flux_model = lc_model()
						plot_LC_qlp_model(time, phase, flux, flux_model, phase_model, TOI, TIC, period, Tmag, Rs)
					else:
						plot_LC_qlp(time, phase, flux, TOI, TIC, period, Tmag, Rs)
				
					print(TIC, TOI)
					
				

	
	
			
			

		
		
	
	
	
	
	
	
		




		
	
	
	
	
