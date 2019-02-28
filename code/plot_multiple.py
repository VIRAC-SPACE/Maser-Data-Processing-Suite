#!/usr/bin/env python

import numpy
import scipy
import sys
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
from numpy import ma
from numpy import genfromtxt
import time
import math
import pickle
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
from matplotlib.ticker import MaxNLocator
from matplotlib import rcParams
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

#matplotlib.rcParams["text.latex.preamble"].append(r'\usepackage[dvips]{graphicx}\usepackage{xfrac}')







fin = open(sys.argv[1],"r")

max_flux_limit = 9999

JD_shift = 0 # SHIFT INPUT DATE TO THIS VALUE
jd_min = float(sys.argv[4])
jd_max = float(sys.argv[5])
vmin = float(sys.argv[6]) # km/s
vmax = float(sys.argv[7]) # km/s
jd_first = sys.argv[8]
source_name = sys.argv[9]

vrange = vmax-vmin

sources_vrange = genfromtxt('DB_vrange.csv', names=True, delimiter=',', dtype=None, encoding=None)

i = 0
found_vrange = False
found_ind = -1
for source in sources_vrange:
	if(str(source[0]) == source_name):
		print('We found vrange for source '+source[0]+': from '+str(source[1])+' to '+str(source[2]))
		found_vrange = True
		found_ind = i
		vmin = source[1]
		vmax = source[2]
	i = i + 1



velocity = []
observed_flux = []
observed_time = []

print('Reading data...\n')

while True:
	line = fin.readline()
	if(len(line)) == 0:
		break
	if ( line.find("#")==-1 ):
		line = line.split()
		velocity.append( float(line[0]) )
		if (float(line[1]) > max_flux_limit ):
			observed_flux.append( max_flux_limit )
		else:
			if float(line[1]) < 0:
				observed_flux.append( 0.001 )
			else:
				observed_flux.append( float(line[1])) 
    
		observed_time.append( float(line[2])+JD_shift )
#		if (float(line[2]) == 316) and ( float(line[0])== -9.41):
#			print float(line[1])	
			
fin.close()
velocity = numpy.array(velocity)
observed_flux = numpy.array(observed_flux)
observed_time = numpy.array(observed_time)


print('Initialize arrays...\n')

x = numpy.arange( vmin, vmax, 0.01 )
y = numpy.arange( 0, jd_max, 0.5 )


X, Y = numpy.meshgrid(x,y)


print('Regridding data...\n')


#Z = matplotlib.mlab.griddata(velocity,observed_time,observed_flux, X, Y, interp='linear')

Z = scipy.interpolate.griddata((velocity, observed_time), observed_flux, (X[:], Y[:]), method='linear')

plt.figure(figsize=(4.2, 8), dpi=100)

plt.xlabel('Velocity (km/s)')
plt.ylabel('JD (days) - '+jd_first) # Zero = 2455479.5

axes = plt.gca()
axes.set_xlim([vmin,vmax])
axes.set_ylim([jd_min,jd_max])

#axes.set_aspect(0.05)

print('Plot data...\n')

#Z = numpy.where(Z > 1e-5, Z, 1e-30)


# Draw contours
#CS = plt.contourf(X, Y, Z, 500, cmap='jet')

CS = plt.pcolor(X, Y, 10*scipy.log10(Z), cmap='spectral', vmin=-0.5, vmax=2)

cbar = plt.colorbar(CS)
cbar.set_clim(vmin=0) # ,vmax=max_flux_limit
cbar.ax.set_ylabel(r'$Flux~(\mathrm{Jy})$')
cbar.locator = MaxNLocator( nbins = 50 )

text_file = open(sys.argv[2], "r")
days = text_file.readline().rstrip().split(',')
days = [float(i) for i in days]
dates = text_file.readline().rstrip().split(',')
dates_days = text_file.readline().rstrip().split(',')
dates_days = [float(i) for i in dates_days]

#print(days)
#print(dates)
#print(dates_days)

text_file.close()

#days = [0.0000,3.9326,7.9618,10.9563,14.9840,24.9847,29.9583,33.9125,37.8570,41.8479,45.9153,53.8861,54.8799,58.8104,62.9507,68.8340,73.8083,79.8076,83.8729,104.8327,117.8257,126.7785,130.7771,138.5660,146.6479,152.5063,159.5368,163.4070,166.6410,167.4979,171.5035,174.6201,179.4854,183.4924,185.6410,186.4875,188.5500]
#
#dates = ["2018-07-31 21:11","2018-08-27 17:55","2018-09-17 20:23","2018-10-29 17:33","2018-12-10 13:07","2018-12-31 09:31","2019-01-18 12:57"]
#
#dates_days = [14.9840,41.8479,62.9507,104.8327,146.6479,167.4979,185.6410]
#
#
for iday in days:
#	#CS = plt.plot([vmin, -4.6], [iday+0.0, iday+0.0], color='w', linewidth=0.1, alpha=0.7)
	#print(iday)
	CS = plt.plot([vmax-math.fabs(vmax-vmin)/8.0, vmax], [iday+0.0, iday+0.0], color='w', linewidth=0.5, alpha=0.9)
#

for i in range(0,len(dates_days)-1):
	iday = dates_days[i]
	axes.text(vmin+0.2, iday, dates[i], fontsize=5, color='white', horizontalalignment='left', alpha= 1.0)

axes.text(0.01,1.01, dates[len(dates_days)-1], fontsize=5, color='black', horizontalalignment='left', alpha= 1.0, transform=axes.transAxes)
#
#
#skip_days = [6.0986,26.1021,127.9042,139.4798,146.4507,153.0514]
#skip_days_len = len(skip_days)
#for x in range(0, skip_days_len, 2):
#	print(skip_days[x],skip_days[x+1])
#	plt.Rectangle((vmin, skip_days[x]),math.fabs(vmax-vmin),math.fabs(skip_days[x+1]-skip_days[x]),edgecolor='None', facecolor='white', alpha=0.9)
#
#

print('Saving figure...\n')

#saving to png file
plt.savefig("maps/" + sys.argv[3], dpi=700, format='png')
#plt.show()




