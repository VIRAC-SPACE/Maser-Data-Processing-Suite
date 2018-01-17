#! /usr/bin/python
from __future__ import division

import os
import sys
import numpy as np
from numpy.polynomial.chebyshev import chebfit
from matplotlib import pyplot as plt
from time import strptime
import scipy.constants
from astropy.modeling import models, fitting
from astropy.modeling.polynomial import Chebyshev1D
import peakutils
from peakutils.plot import plot as pplot

from experimentsLogReader import ExperimentLogReader

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def usage():
    print ('Usage: '+sys.argv[0]+' log file' + 'data file' +  " point counts")

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
    #print("velocitySoure ", velocitySoure, " ObservedFrequency ", ObservedFrequency  )
    return velocitySoure

def is_outlier(points, thresh=3.5):
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
      
    data = np.fromfile(sys.argv[2], dtype="float64", count=-1, sep=" ") .reshape((file_len(sys.argv[2]),5))
    data = np.delete(data, (0), axis=0)
    
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
    
    plt.plot(data[:, [0]], data[:, [1]] * calibration(location, Systemtemperature1u), 'ro')
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("1u")
    plt.show()
    
    plt.plot(data[:, [0]], data[:, [2]] * calibration(location, Systemtemperature9u), 'ro')
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("9u")
    plt.show()
    
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
    
    m = 0
    n = dataPoints
     
    FreqStart = scan["FreqStart"] 
  
    x = dopler((data[:, [0]] + FreqStart) * (10 ** 6), VelTotal)
    xlist = list()
    
    y1 = data[:, [1]] * calibration(location, Systemtemperature1u)
    y1list = list()
    
    y2 = data[:, [2]] * calibration(location, Systemtemperature1u)
    y2list = list()
    
    xarray = np.ones(dataPoints)
    y1array = np.ones(dataPoints)
    y2array = np.zeros(dataPoints)
    
    for i in range(0,dataPoints):
        xarray[i] = x[i]
    
    y1array = np.zeros(dataPoints)
    for j in range(0,dataPoints):
        y1array[j] = y1[j]
    
    for k in range(0,dataPoints):
        y2array[k] = y2[k]
    
    middle = int(dataPoints/2)
    
    a = int(middle*0.88)
    b = int(middle*1.075)
    
    ### u1
    
    # polyfit
    z = np.polyfit(np.append(xarray[0:a],xarray[b:dataPoints]), np.append(y1array[0:a],y1array[b:dataPoints]), 9)
    p = np.poly1d(z)
    
    # Fit the data using a box model
    t_init = models.Trapezoid1D(amplitude=1., x_0=0., width=1., slope=0.5)
    fit_t = fitting.LevMarLSQFitter()
    t1 = fit_t(t_init, np.append(xarray[0:a],xarray[b:dataPoints]),  np.append(y1[0:a],y1[b:dataPoints]) )
    
    # Fit the data using a Gaussian
    g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g1 = fit_g(g_init, np.append(xarray[0:a],xarray[b:dataPoints]),  np.append(y1[0:a],y1[b:dataPoints]))
    
    # Fit the data using a Chebyshev astro py
    ceb = Chebyshev1D(9, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
    fit_ceb = fitting.LevMarLSQFitter()
    ceb_1 = fit_ceb(ceb, np.append(xarray[0:a],xarray[b:dataPoints]),  np.append(y1[0:a],y1[b:dataPoints]))
    #cb1 = y2array -  ceb_1(y2array)
    
    # Fit the data using a Chebyshev numpy
    #ceb1 = chebfit(np.append(xarray[0:a],xarray[b:dataPoints]), np.append(y1[0:a],y1[b:dataPoints]), 9, rcond=None, full=False, w=None)
    #c1 = np.poly1d(ceb1)
   
    ### u9
    
    # polyfit
    z2 = np.polyfit(np.append(xarray[0:a],xarray[b:dataPoints]), np.append(y2array[0:a],y2array[b:dataPoints]), 9)
    p2 = np.poly1d(z2)
    
    # Fit the data using a box model
    t_init = models.Trapezoid1D(amplitude=1., x_0=0., width=1., slope=0.5)
    fit_t = fitting.LevMarLSQFitter()
    t2 = fit_t(t_init, np.append(xarray[0:a],xarray[b:dataPoints]),  np.append(y2[0:a],y2[b:dataPoints]))
    
    # Fit the data using a Gaussian
    g_init = models.Gaussian1D(amplitude=1., mean=0, stddev=1.)
    fit_g = fitting.LevMarLSQFitter()
    g2 = fit_g(g_init, np.append(xarray[0:a],xarray[b:dataPoints]), np.append(y2[0:a],y2[b:dataPoints]))
    
    # Fit the data using a Chebyshev astro py
    ceb = Chebyshev1D(9, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
    fit_ceb = fitting.LevMarLSQFitter()
    ceb_2 = fit_ceb(ceb, np.append(xarray[0:a],xarray[b:dataPoints]),  np.append(y2[0:a],y2[b:dataPoints]))
     
    # Fit the data using a Chebyshev numpy
    #ceb2 = chebfit(np.append(xarray[0:a],xarray[b:dataPoints]), np.append(y2[0:a],y2[b:dataPoints]), 9, rcond=None, full=False, w=None)
    #c2 = np.poly1d(ceb2)
    
    #Polinom ploting 
      
    #1u
    plt.figure()
    #plt.plot(xarray, p(xarray), 'r', label='Poly Fit')
    #plt.plot(xarray, t1(xarray), 'b', label='Trapezoid')
    #plt.plot(xarray, g1(xarray), 'g', label='Gaussian')
    plt.plot(xarray, ceb_1(xarray), 'y', label='Chebyshev_1')
    #plt.plot(np.append(xarray[0:a],xarray[b:dataPoints]), c1(np.append(xarray[0:a],xarray[b:dataPoints])), 'y', label='Chebyshev_2')
    #plt.plot(x, y1, 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend(loc=2)
    plt.title("1u Polarization for source " + scan["sourceName"] + " scan " + str(scanNumber))   
    plt.show()
    
    #9u
    plt.figure()
    #plt.plot(xarray, p2(xarray), 'r', label='Poly Fit')
    #plt.plot(xarray, t2(xarray), 'b', label='Trapezoid')
    #plt.plot(xarray, g2(xarray), 'g', label='Gaussian')
    plt.plot(xarray, ceb_2(xarray), 'y', label='Chebyshev_1')
    #plt.plot(np.append(xarray[0:a],xarray[b:dataPoints]), c2(np.append(xarray[0:a],xarray[b:dataPoints])), 'y', label='Chebyshev_2')
    #plt.plot(x, y2, 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend(loc=2)
    plt.title("9u Polarization for source " + scan["sourceName"] + " scan  " + str(scanNumber))
    plt.show()
    
    #Local maximum ploting
    
    #indexes_for_ploy_fitt = peakutils.indexes(y1array - p(xarray), thres=0.02/max(y1array - p(xarray)), min_dist=100)
    #indexes_for_trap = peakutils.indexes(y1array - t1(xarray), thres=0.02/max(y1array - t1(xarray)), min_dist=3)
    #indexes_for_gauss = peakutils.indexes(y1array - g1(xarray), thres=0.02/max(y1array - g1(xarray)), min_dist=3)
    indexes_for_ceb = peakutils.indexes(y1array - ceb_1(xarray), thres=0.02/max(y1array - ceb_1(xarray)), min_dist=50)
    
    #peaks_x = peakutils.interpolate(xarray, y1array, ind=indexes_for_ploy_fitt)
    
    #1u
    fig = plt.figure()
    #plt.plot(xarray, y1array - p(xarray), 'r', label='Poly Fit')
    #plt.plot(xarray, y1array - t1(xarray), 'b', label='Trapezoid')
    #plt.plot(xarray, y1array - g1(xarray), 'g', label='Gaussian')
    plt.plot(xarray, y1array - ceb_1(xarray), 'y', label='Chebyshev_1')
    
    #plt.plot(xarray[indexes_for_ploy_fitt], y1array[indexes_for_ploy_fitt], 'ro', label="Local Maximums for poly fit")
    
    '''
    ax = fig.add_subplot(111)
    for xy in zip(xarray[indexes_for_ploy_fitt], y1array[indexes_for_ploy_fitt]):                        # <--
        ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
    '''
       
    plt.plot(xarray[indexes_for_ceb], y1array[indexes_for_ceb], 'yo', label="Local Maximums for Chebyshev")
    
    ax = fig.add_subplot(111)
    for xy in zip(xarray[indexes_for_ceb], y1array[indexes_for_ceb]):                        
        ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
    
    #pplot(xarray, y1array, indexes_for_ploy_fitt)
    #pplot(xarray, y1array, indexes_for_trap)
    #pplot(xarray, y1array, indexes_for_gauss)
    #pplot(xarray, y1array, indexes_for_ceb)
    
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend(loc=2)
    plt.title("1u Polarization for source " + scan["sourceName"] + " scan " + str(scanNumber))   
    plt.show()
    
    #indexes_for_ploy_fitt2 = peakutils.indexes(y1array - p2(xarray), thres=0.02/max(y1array - p2(xarray)), min_dist=3)
    #indexes_for_trap2 = peakutils.indexes(y1array - t2(xarray), thres=0.02/max(y1array - t2(xarray)), min_dist=3)
    #indexes_for_gauss2 = peakutils.indexes(y1array - g2(xarray), thres=0.02/max(y1array - g2(xarray)), min_dist=3)
    indexes_for_ceb2 = peakutils.indexes(y1array - ceb_2(xarray), thres=0.02/max(y1array - ceb_2(xarray)), min_dist=50)
    
    #9u
    plt.figure()
    #plt.plot(xarray, y2array - p2(xarray), 'r', label='Poly Fit')
    #plt.plot(xarray, y2array - t2(xarray), 'b', label='Trapezoid')
    #plt.plot(xarray, y2array - g2(xarray), 'g', label='Gaussian')
    plt.plot(xarray, y2array - ceb_2(xarray), 'y', label='Chebyshev_1')
    
    #pplot(xarray, y1array, indexes_for_ploy_fitt2)
    #pplot(xarray, y1array, indexes_for_trap2)
    #pplot(xarray, y1array, indexes_for_gauss2)
    pplot(xarray, y1array, indexes_for_ceb2)
    
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend(loc=2)
    plt.title("9u Polarization for source " + scan["sourceName"] + " scan  " + str(scanNumber))
    plt.show()
    
    
    sys.exit(0)