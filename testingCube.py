import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.visualization import (MinMaxInterval, SqrtStretch, ImageNormalize)
from astropy.time import Time
import os
import datetime
from datacube import Datacube
import xarray

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

source = 'cepa'
dir = '/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/output/'

velocitieList = list()
amplitudeList = list()
dateList = list()
for outputFile in os.listdir(dir):
    if outputFile.startswith(source):
        file = dir + outputFile
        data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),4))
        velocitieList.append(data[:, [0]][0])
        amplitudeList.append(data[:, [3]][0])
        time = " ".join(outputFile.split(".")[0].split("_")[1:7])
        time=datetime.datetime.strptime(time, "%H %M %S %d %b %Y")
        dateList.append(time)

def convertDatetimeObjectToMJD(time):
            time=time.isoformat()
            t=Time(time, format='isot')
            return t.mjd

#velocitieList = [v[0] for v in velocitieList]      
dateList = [convertDatetimeObjectToMJD(date) for date in dateList]



'''
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#psm = axs.pcolormesh(velocitieList, dateList, amplitudeList)
X, Y = np.meshgrid(velocitieList, dateList)
ax.plot_surface(X, Y, np.array(amplitudeList), rstride=10, cstride=10)
#fig.colorbar(psm)
plt.show()
'''


'''
image = [dateList, velocitieList]

# Create an ImageNormalize object
norm = ImageNormalize(image)


# or equivalently using positional arguments
# norm = ImageNormalize(image, MinMaxInterval(), SqrtStretch())

# Display the image
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
im = ax.imshow(image, origin='lower', norm=norm)
fig.colorbar(im)
plt.show()
'''

