'''Module containg various useful functions for Ed's exoplanet/TESS/NGTS PhD'''
import numpy as np
from matplotlib import pyplot as plt
import argparse
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

def plot_LC_spoc(time, raw_flux, flux, TIC, Tmag, Rs):
	
	fig = plt.figure(figsize=(15, 10))
	
	#Top Subplot: Raw SAP Flux
	ax1 = fig.add_subplot(211)
	ax1.plot(time, raw_flux, 'ko', markersize=1.5)
	ax1.set_ylabel('Raw SAP Flux [e$^-$ / s]', **axis_font)
	ax1.set_xlabel('Time [BJD - 2457000]', **axis_font) 
	ax1.tick_params(labelsize=10)
	ax1.set_title('Tmag = {:.3f} ; Rs = {};  TIC ID: {}'.format(Tmag, Rs, TIC), **axis_font)


	#Middle subplot: Unfolded LC
	ax2 = fig.add_subplot(212)
	ax2.plot(time, flux, 'ko', markersize='1.5')
	ax2.set_ylabel('PDC Flux [e$^-$ / s]', **axis_font)
	ax2.set_xlabel('Time [BJD - 2457000]', **axis_font)
	ax2.tick_params(labelsize=10)
	
#	if save:
#		plt.savefig('/home/astro/phrvdf/tess_data_alerts/tess_LC_plots/tess_{}_{}_lc_sector{}.png'.format(TOI, TIC, sector))
#		plt.close()
#	else:
	plt.show()
	

if __name__ == "__main__":
		
    parser = argparse.ArgumentParser()
    parser.add_argument('-fn', '--filename', type=str, nargs='*')
    parser.add_argument('-s', '--save', action='store_true')
	
    args = parser.parse_args()
	
    filenames = args.filename
    save = args.save
	
    axis_font = {'fontname':'DejaVu Sans', 'size':'15'}
		
    for i in range(len(filenames)):
		
        hdul = fits.open(filenames[i])
        hdr = hdul[0].header
		
        TIC = hdr['TICID']							#Loads TIC ID number from the fits header
        Tmag = hdr['TESSMAG']						#Loads TESS mag from fits header
        Rs = hdr['RADIUS']							#Loads stellar radius [solar radius] from fits header
        sector = hdr['SECTOR']						#Loads observing sector data was obtained in
			
        time, flux, fluxerr, time_whole, raw_flux = tess_LC_dataload_spoc(filenames[i])
        flux_normalised = flux/np.median(flux)
		
        plot_LC_spoc(time, raw_flux, flux_normalised, TIC, Tmag, Rs)

        print(TIC)
					
				

	
	
			
			

		
		
	
	
	
	
	
	
		




		
	
	
	
	
