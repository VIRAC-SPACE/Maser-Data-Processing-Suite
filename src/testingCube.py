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

files = ["cepa_18_14_07_22_Dec_2018_IRBENE16_290.dat", 
"cepa_19_08_58_22_Dec_2018_IRBENE16_291.dat",
"cepa_20_03_49_22_Dec_2018_IRBENE16_292.dat",
"cepa_20_58_40_22_Dec_2018_IRBENE16_293.dat",
"cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat",
"cepa_00_48_39_23_Dec_2018_IRBENE16_295.dat",
"cepa_01_43_29_23_Dec_2018_IRBENE16_296.dat",
"cepa_02_38_20_23_Dec_2018_IRBENE16_297.dat",
"cepa_03_33_10_23_Dec_2018_IRBENE16_298.dat",
"cepa_04_28_01_23_Dec_2018_IRBENE16_299.dat",
"cepa_05_22_51_23_Dec_2018_IRBENE16_300.dat",
"cepa_06_17_42_23_Dec_2018_IRBENE16_301.dat",
"cepa_07_12_32_23_Dec_2018_IRBENE16_302.dat",
"cepa_08_07_23_23_Dec_2018_IRBENE16_303.dat",
"cepa_16_57_44_23_Dec_2018_IRBENE16_304.dat",
"cepa_17_52_37_23_Dec_2018_IRBENE16_305.dat",
"cepa_18_55_58_23_Dec_2018_IRBENE16_306.dat",
"cepa_19_50_50_23_Dec_2018_IRBENE16_307.dat",
"cepa_20_45_43_23_Dec_2018_IRBENE16_308.dat",
"cepa_21_40_34_23_Dec_2018_IRBENE16_309.dat"]

hardCodedDates = ["18_14_07_22_Dec_2018",
"19_08_58_22_Dec_2018",
"20_03_49_22_Dec_2018",
"20_58_40_22_Dec_2018",
"23_53_48_22_Dec_2018",
"00_48_39_23_Dec_2018",
"01_43_29_23_Dec_2018",
"02_38_20_23_Dec_2018",
"03_33_10_23_Dec_2018",
"04_28_01_23_Dec_2018",
"05_22_51_23_Dec_2018",
"06_17_42_23_Dec_2018",                  
"07_12_32_23_Dec_2018",                  
"08_07_23_23_Dec_2018",                  
"16_57_44_23_Dec_2018",
"17_52_37_23_Dec_2018",
"18_55_58_23_Dec_2018",
"19_50_50_23_Dec_2018",
"20_45_43_23_Dec_2018",                   
"21_40_34_23_Dec_2018"]
format = '%H_%M_%S_%d_%b_%Y'

dates = [convertDatetimeObjectToMJD(datetime.datetime.strptime(d, format)) for d in hardCodedDates]

dir = '/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/output/'
print("creating data...")
velocity = list()
observed_flux = list()
observed_time = list()
for file in files:
    data = np.fromfile(dir + file, dtype="float128", count=-1, sep=" ").reshape((file_len(dir + file),4))
    velocity.append(correctNumpyReadData(data[:, [0]]))
    observed_flux.append(correctNumpyReadData(data[:, [3]]))
    observed_time.append([dates[files.index(file)]] * len(correctNumpyReadData(data[:, [0]])))
    
velocity = np.array(velocity).reshape(1, 40940)[0]
observed_flux = np.array(observed_flux).reshape(1, 40940)[0]
observed_time = np.array(observed_time).reshape(1, 40940)[0]
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


