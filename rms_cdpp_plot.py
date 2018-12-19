'''Script to plot RMS CDPP on 0.5 hr time scale against TESS magnitude 
	for all spoc lightcurves'''
	
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import argparse

def data_load(filename):
	hdul = fits.open(filename)
	
	rms_05 = hdul[1].header['CDPP0_5']
	rms_1 = hdul[1].header['CDPP1_0']
	rms_2 = hdul[1].header['CDPP2_0']
	
	Tmag = hdul[0].header['TESSMAG']
	
	return rms_05, rms_1, rms_2, Tmag
	
if __name__ == "__main__":
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-fn', '--filename', type=str, nargs='*')

	args = parser.parse_args()
	
	fn = args.filename

	axis_font = {'fontname':'DejaVu Sans', 'size':'20'}
	
	rms5 = np.zeros(len(fn))
	rms1 = np.zeros(len(fn))
	rms2 = np.zeros(len(fn))
	Tmags = np.zeros(len(fn))
	
	for i in range(len(fn)):
		
		rms5[i], rms1[i], rms2[i], Tmags[i] = data_load(fn[i])
	
	plt.figure()
	
	plt.plot(Tmags, rms5, 'ko', markersize=5, label='TESS - 0.5 Hour')
	plt.plot(14.04, 2198.769, 'ro', markersize=5, label='NGTS - 0.5 Hour')
#	plt.plot(Tmags, rms1, 'bo', markersize=5, label='TESS - 1 Hour')
#	plt.plot(Tmags, rms2, 'go', markersize=5, label='TESS - 2 Hour')
	plt.xlabel('TESS Magnitude', **axis_font)
	plt.ylabel('RMS [ppm]', **axis_font)
	plt.title('RMS CDPP - multiple timescales', **axis_font)
	
	plt.legend(fontsize=20)
	
	plt.show()
	
