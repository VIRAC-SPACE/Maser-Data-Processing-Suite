import numpy as np
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
import datetime

rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True

from help import *

files = ["cepa_00_01_03_15_Nov_2019_IRBENE_418.dat",
]

hardCodedDates = ["00_01_03_15_Nov_2019",
]
format = '%H_%M_%S_%d_%b_%Y'

dates = [convertDatetimeObjectToMJD(datetime.datetime.strptime(d, format)) for d in hardCodedDates]

dir = '/home/janis/Documents/maser/output/cepa/6668/'
print("creating data...")
velocity = list()
observed_flux = list()
observed_time = list()
for file in files:
    data = np.fromfile(dir + file, dtype="float128", count=-1, sep=" ").reshape((file_len(dir + file),4))
    velocity.append(correctNumpyReadData(data[:, [0]]))
    observed_flux.append(correctNumpyReadData(data[:, [3]]))
    observed_time.append([dates[files.index(file)]] * len(correctNumpyReadData(data[:, [0]])))
    
velocity = np.array(velocity).reshape(1, 2047)[0]
observed_flux = np.array(observed_flux).reshape(1, 2047)[0]
observed_time = np.array(observed_time).reshape(1, 2047)[0]
print("Data created")

print("creating mesh grid ...")
x = np.arange(np.min(velocity), np.max(velocity), 0.01)
y = np.arange(0, np.max(observed_time), 0.5)
X, Y = np.meshgrid(x,y)
print("mesh grid created")

print("interpolation ...")
Z = scipy.interpolate.griddata((velocity, observed_time), observed_flux, (X[:], Y[:]), method='linear')
print("interpolation done")

plt.figure(figsize=(4.2, 8), dpi=100)

plt.xlabel('Velocity (km/s)')
plt.ylabel('JD (days) - ')

axes = plt.gca()
axes.set_xlim([np.min(velocity),np.max(velocity)])
axes.set_ylim([np.min(observed_time), np.max(observed_time)])

CS = plt.contourf(X, Y, Z, 500, cmap='jet')

#CS = plt.pcolor(X, Y, scipy.log10(Z), cmap='spectral', vmin=-0.5, vmax=2)

cbar = plt.colorbar(CS)
cbar.set_clim(vmin=0) # ,vmax=max_flux_limit
cbar.ax.set_ylabel(r'$Flux~(\mathrm{Jy})$')
cbar.locator = MaxNLocator( nbins = 50 )

plt.show()


