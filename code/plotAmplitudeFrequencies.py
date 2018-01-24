#! /usr/bin/python
from __future__ import division

import os
import sys
import numpy as np
from matplotlib import pyplot as plt, axis
from matplotlib.widgets import *
from time import strptime
import scipy.constants
from astropy.modeling import models, fitting
from astropy.modeling.polynomial import Chebyshev1D
from astropy.convolution import Gaussian1DKernel, convolve
import peakutils
import pylab

from experimentsLogReader import ExperimentLogReader

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def usage():
    print ('Usage: ' + sys.argv[0] + ' log file' + 'data file')

def calibration(station, Tsys):
    scale = 1
    if station == "IRBENE":
        scale = 12 
        
    elif station == "IRBENE16":
        scale = 26
    return scale*Tsys
    
def dopler(ObservedFrequency, velocityReceiver):
    c = scipy.constants.speed_of_light
    f0 = 6668519200 # Hz 
    velocitySoure = (-((ObservedFrequency/f0)-1)*c + (velocityReceiver * 1000))/1000
    return velocitySoure

def is_outlier(points, thresh=4.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score < thresh
    
if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    
    experimentName = sys.argv[2].split("_")[0].split("/")[1]
    
    data = np.fromfile(sys.argv[2], dtype="float64", count=-1, sep=" ") .reshape((file_len(sys.argv[2]),5))
    data = np.delete(data, (0), axis=0) #izdzes masiva primo elementu
    
    #slikto punktu izdzesana
    outliersMask = is_outlier(data[:, [0]])
    data = data[outliersMask]
    
    dataPoints = data.shape[0]
    
    logs  = ExperimentLogReader(sys.argv[1]).getLgs()
    scanNumber = sys.argv[2].split(".")[0].split("_")[1][1:len(sys.argv[2])]
    scan = logs[scanNumber]
    Systemtemperature1u = float(scan["Systemtemperature"][0])
    Systemtemperature9u = float(scan["Systemtemperature"][1])
    
    location = logs["location"]
    
    xdata = data[:, [0]]
    y1data = data[:, [1]] * calibration(location, Systemtemperature1u)
    y2data = data[:, [2]] * calibration(location, Systemtemperature9u)
    
    y1array = np.zeros(dataPoints)
    y2array = np.zeros(dataPoints)
    
    for j in range(0,dataPoints):
        y1array[j] = y1data[j]
    
    for k in range(0,dataPoints):
        y2array[k] = y2data[k]
    
    g1 = Gaussian1DKernel(stddev=np.std(y1array), x_size=19, mode='center', factor=100)
    g2 = Gaussian1DKernel(stddev=np.std(y2array), x_size=19, mode='center', factor=100)
    
    z1 = convolve(y1array, g1, boundary='extend')
    z2 = convolve(y2array, g2, boundary='extend')
    
    #avota izgriesana
    middle = int(dataPoints/2) #vidusunks
    a = int(middle*0.88)
    b = int(middle*1.075)
    
    #galu nogriesana
    m = 0
    n = dataPoints
    
    #1u
    fig = pylab.gcf()
    fig.canvas.set_window_title("Data filtering for experiment " +  experimentName)
    fig.set_size_inches(10.5, 10.5)
    plt.suptitle("source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber), fontsize=16)
    plt.subplot(121)
    plt.subplots_adjust(bottom=0.3, wspace = 0.35)
    plt.plot(xdata, z1, 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('Frequency Mhz')
    plt.ylabel ('Flux density (Jy)')
    plt.legend(loc=2)
    
    #pievieno papildus asi data punktiem
    plt.twiny()
    plt.xlabel("Data points")
    plt.tick_params(axis="x")
    plt.xticks(range(0, dataPoints + 512, 512))
    
    plt.title("1u Polarization",  y=1.08) 
    
    #9u
    plt.subplot(122)
    plt.plot(xdata, z2, 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('Frequency Mhz')
    plt.ylabel ('Flux density (Jy)')
    plt.legend(loc=2)
    
    #pievieno papildus asi data punktiem
    plt.twiny()
    plt.xlabel("Data points")
    plt.tick_params(axis="x")
    plt.xticks(range(0, dataPoints + 512, 512))
    plt.title("9u Polarization",  y=1.08) 
    
    #sliders
    mAxes = plt.axes([0.10, 0.15, 0.65, 0.03])
    nAxes  = plt.axes([0.10, 0.10, 0.65, 0.03])

    mSlider=Slider(mAxes, 'M', m, a-1, valinit=m)
    nSlider=Slider(nAxes , 'N', b-1, n, valinit=b-1)

    def update(val):
        global m,n
        m=int(mSlider.val)
        n=int(nSlider.val)
    
    mSlider.on_changed(update)
    nSlider.on_changed(update)
    
    plt.show()
    
    print m,n
    
    timeStr = scan['startTime'].replace(":", " ")
    dateStrList = scan['dates'].split()
    dateStrList[1] = strptime(dateStrList[1],'%b').tm_mon
    dateStr = str(dateStrList[2]) + " " + str(dateStrList[1]) + " " + str(dateStrList[0])
    RaStr = " ".join(scan["Ra"])
    DecStr = " ".join(scan["Dec"])
    dopsetPar= dateStr + " " + timeStr + " " + RaStr + " " + DecStr
    
    os.system("code/dopsetpy_v1.5 "+dopsetPar)
    
    # dopsetpy parametru nolasisana
    with open('lsrShift.dat') as openfileobject:
        for line in openfileobject:
            Header = line.split(';')
            vards = Header[0]
            if vards == "Date":
                dateStr = Header[1]
            elif vards == "Time":
                laiks = Header[1]
            elif vards == "RA":
                RaStr = Header[1]
            elif vards == "DEC":
                DecStr = Header[1]
            elif vards == "Source":
                Source = Header[1]
            elif vards == "LSRshift":
                lsrShift = Header[1]
            elif vards == "MJD":
                mjd = Header[1]
                print "MJD: \t",mjd
            elif vards == "Vobs":
                Vobs = Header[1]
                print "Vobs: \t",Vobs
            elif vards == "AtFreq":
                AtFreq = Header[1]
                print "At Freq: \t",AtFreq
            elif vards == "FreqShift":
                FreqShift = Header[1]
                print "FreqShift: \t",FreqShift
            elif vards == "VelTotal":
                VelTotal = float(Header[1])
                print "VelTotal: \t",VelTotal
            #Header +=1

    Vobs = float(Vobs)
    lsrCorr = float(lsrShift)*1.e6 # for MHz
      
    FreqStart = scan["FreqStart"] 
    
    #Parveido frekvenci par atrumu
    x = dopler((xdata + FreqStart) * (10 ** 6), VelTotal)
    
    xarray = np.zeros(dataPoints)
    
    for i in range(0,dataPoints):
        xarray[i] = x[i]
    
    #polinomu apreikinasana
    
    # Fit the data using a Chebyshev astro py
    ceb = Chebyshev1D(9, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
    fit_ceb = fitting.LevMarLSQFitter()
    
    print m,a, 
    print b,n
    
    #m = 900
    #n = 3200
    
    ### u1
    ceb_1 = fit_ceb(ceb, np.append(xarray[m:a], xarray[b:n]),  np.append(z1[m:a],z1[b:n]))
   
    ### u9
    ceb_2 = fit_ceb(ceb, np.append(xarray[m:a], xarray[b:n]),  np.append(z2[m:a],z2[b:n]))
     
  
    #Polinom ploting 
      
    #1u
    fig = pylab.gcf()
    fig.canvas.set_window_title("Polynomial")
    plt.plot(xarray[m:n], ceb_1(xarray[m:n]), 'r', label='Chebyshev polynomial')
    plt.plot(np.append(xarray[m:a], xarray[b:n]), np.append(z1[m:a], z1[b:n]), 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("1u Polarization for source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber))  
    plt.show()
    
    #9u
    fig = pylab.gcf()
    fig.canvas.set_window_title("Polynomial")
    plt.plot(xarray[m:n], ceb_2(xarray[m:n]), 'r', label='Chebyshev polynomial')
    plt.plot(np.append(xarray[m:a], xarray[b:n]), np.append(z2[m:a], z2[b:n]), 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("9u Polarization for source " + scan["sourceName"].split(",")[0] + " scan  " + str(scanNumber))
    plt.show()
     
    #Local maximum ploting
    
    thres=0.1
    
    y1values = z1[m:n] - ceb_1(xarray[m:n])
    y2values = z2[m:n] - ceb_2(xarray[m:n])
    
    #indexsu apreikinasana
    indexes_for_ceb = peakutils.indexes(y1values, thres=thres, min_dist=10)
    indexes_for_ceb2 = peakutils.indexes(y2values, thres=thres, min_dist=10)
    
    #1u
    fig = pylab.gcf()
    fig.canvas.set_window_title("Local Maximums")

    plt.plot(xarray[m:n], y1values, 'b', label='Signal - polynomial')  
    plt.plot(xarray[m:n][indexes_for_ceb], y1values[indexes_for_ceb], 'dr', label="Local Maximums for signal")
    
    ax = fig.add_subplot(111)
    for xy in zip(xarray[m:n][indexes_for_ceb], y1values[indexes_for_ceb]):                        
        ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
    
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("1u Polarization for source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber))  
    plt.show()
    
    #9u
    fig = pylab.gcf()
    fig.canvas.set_window_title("Local Maximums")
    
    plt.plot(xarray[m:n], y2values, 'b', label='Signal - polynomial')
    plt.plot(xarray[m:n][indexes_for_ceb2], y2values[indexes_for_ceb2], 'dr', label="Local Maximums for signal")
    
    xa = fig.add_subplot(111)
    for yx in zip(xarray[m:n][indexes_for_ceb2], y2values[indexes_for_ceb2]):                        
        xa.annotate('(%.2f, %.1f)' % yx, xy=yx, textcoords='data')
    
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("9u Polarization for source " + scan["sourceName"].split(",")[0] + " scan  " + str(scanNumber))
    plt.show()
    
    sys.exit(0)