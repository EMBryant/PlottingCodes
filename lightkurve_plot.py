from __future__ import division
from lightkurve import TessLightCurveFile
import matplotlib.pyplot as plt
import sys
import pandas
import rogues

file_name = sys.argv[1]
TOI = float(sys.argv[2])
df = pandas.read_csv('../TESS/TOIs_Sec1_20180905.csv', index_col='toi_id')		#.csv file containing info on parameters (period, epoch, ID, etc.) of all TOIs


epoch = df.loc[TOI, 'Epoc'] 		      	#Time of first transit centre [BJD - 2457000]
period = df.loc[TOI, 'Period']			  	#Orbital Period [days]
T_dur = df.loc[TOI, 'Duration']				#Transit duration [hours]
pipeline = df.loc[TOI, 'src']				#Pipeline used to reduce data - so that can call correct columns
TIC = df.loc[TOI, 'tic_id']              	#TIC ID for the object - used for plot title
comments = df.loc[TOI, 'Comment']			#Any existing comments on the object

lcf = TessLightCurveFile(file_name)
sap_flux = lcf.SAP_FLUX
lc = lcf.PDCSAP_FLUX

flat_lc = lc.flatten(window_length=1001)

phase, phase_days = rogues.phase_fold(flat_lc.time, epoch, period, 0.75)

axis_font = {'fontname':'DejaVu Sans', 'size':'20'}


fig = plt.figure()

ax1 = fig.add_subplot(311)

ax1.plot(sap_flux.time, sap_flux.flux, 'ko', markersize=1.5)
ax1.set_ylabel('Raw SAP Flux [e$^-$ / s]', **axis_font)
ax1.set_xlabel('Time [BJD - 2457000]', **axis_font) 
ax1.set_title('TOI: {} ;  TIC ID: {} ;  Period: {} days'.format(TOI, TIC, period), **axis_font)

ax2 = fig.add_subplot(312)

ax2.plot(flat_lc.time, flat_lc.flux, 'ko', markersize=1.5)
ax2.set_ylabel('Relative Flux', **axis_font)
ax2.set_xlabel('Time [BJD - 2457000]', **axis_font)

ax3 = fig.add_subplot(313)

ax3.plot(phase, flat_lc.flux, 'ko', markersize=1.5)
ax3.set_ylabel('Relative Flux', **axis_font)
ax3.set_xticks([-0.25, 0.0, 0.25, 0.5, 0.75])
ax3.set_xlabel('Phase', **axis_font)

plt.show()
