#! /usr/bin/python
from __future__ import division

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
#from PyAstronomy import pyasl
from time import strptime
import scipy.constants
#import jdcal
#from jdcal import gcal2jd
import datetime

from math import pi
from numpy import cos, sin

from experimentsLogReader import ExperimentLogReader

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def usage():
    print ('Usage: '+sys.argv[0]+' log file' + 'data file')
    

    
def dopler(ObservedFrequency, velocityReceiver):
    c = scipy.constants.speed_of_light
    f0 = 6668519200 # Hz 
    lo = 6100000000 # Hz
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
    
    logs  = ExperimentLogReader(sys.argv[1]).getLgs()
    scanNumber = sys.argv[2].split(".")[0].split("_")[1][1:len(sys.argv[2])]
    scan = logs[scanNumber]
    Systemtemperature1u = float(scan["Systemtemperature"][0])
    Systemtemperature9u = float(scan["Systemtemperature"][1])
    
    
    timeStr = scan['startTime'].replace(":", " ")
    dateStrList = scan['dates'].split()
    dateStrList[1] = strptime(dateStrList[1],'%b').tm_mon
    dateStr = str(dateStrList[2]) + " " + str(dateStrList[1]) + " " + str(dateStrList[0])
    RaStr = " ".join(scan["Ra"])
    DecStr = " ".join(scan["Dec"])
    dopsetPar= dateStr + " " + timeStr + " " + RaStr + " " + DecStr
    print scan
    print dopsetPar
    
    outliersMask = is_outlier(data[:, [0]])
    data = data[outliersMask]
      
    scale1U = 1
    scale9U = 1
    
    location = logs["location"]
    
    if location == "IRBENE":
        scale1U = 12
        scale9U = 12
        
    elif location == "IRBENE16":
        scale1U = 26
        scale9U = 26
    ''' 
    plt.scatter(data[:, [0]], data[:, [1]] * Systemtemperature1u * scale1U)
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("1u")
    plt.show()
    
    plt.scatter(data[:, [0]], data[:, [2]] * Systemtemperature9u * scale9U)
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("9u")
    plt.show()
     
    longitude = 21.8605
    latitude = 57.5546
    altitude = 87.30
    
    ra = float(scan["Ra"][0]) + float(scan["Ra"][1])/60.0 + float(scan["Ra"][2])/3600.0
    dec = float(scan["Dec"][0]) + float(scan["Dec"][1])/60.0 + float(scan["Dec"][2])/3600.0
    
    date = scan["dates"]
    print("data", date)
    year = int(date.split(" ")[2])
    day = int(date.split(" ")[0])
    month = strptime(date.split(" ")[1],'%b').tm_mon
    startTime = scan["startTime"]
    print("startTime ", startTime)
    hours = int(startTime.split(":")[0])
    minutes = int(startTime.split(":")[1])
    seconds = int(startTime.split(":")[2])
    
    #JDN = (1461 * (year + 4800 + (month - 14)/12))/4 +(367 * (month - 2 - 12 * ((month - 14)/12)))/12 - (3 * ((year + 4900 + (month - 14)/12)/100))/4 + day - 32075
    JDN = sum(gcal2jd(year, month, day))
    jd = JDN + hmsm_to_days(hours, minutes, minutes ) # ((hours-12.0)/24.0) + (minutes/1440.0) + (seconds/86400.0)
    
    print("ra ", ra)
    print("dec ", dec)
    corr, hjd = pyasl.helcorr(longitude, latitude, altitude, \
            ra, dec, jd, debug=True)
   
    print("ra ", scan["Ra"])
    print("dec ", scan["Dec"])
    print("Barycentric correction [km/s]: ", corr)
    print("Heliocentric Julian day: ", hjd)
    print("hmsm_to_days ", hmsm_to_days(hours, minutes, minutes ) )
    
    FreqStart = float(scan["FreqStart"])
    #FreqStart = 6668.83
    print("FreqStart", FreqStart)
    print("from logs ", scan)
    
    
    heli, bary = pyasl.baryvel(jd, deq=2000.0)

    print("Earth's velocity at JD: ", jd)
    print("Heliocentric velocity [km/s]: ", heli)
    print("Barycentric velocity [km/s] : ", bary)
    
    vh, vb = pyasl.baryCorr(jd, ra, dec, deq=2000.0)
    print("Barycentric velocity of Earth toward Sirius: ", vb)
    print("super ", x_keckhelio(ra, dec, epoch=2000.0, jd=jd, tai=None, longitude=longitude, latitude=latitude, altitude=altitude, obs='irbene'))
    
    dvelh, dvelb = pyasl.baryvel(jd, 2000.0)
    print ("dvelh ", dvelh)
    print ("dvelb ", dvelb)
    
    plt.scatter(dopler((data[:, [0]] + FreqStart) * (10 ** 6), corr),      data[:, [1]] * Systemtemperature1u * scale1U)
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend("1u")
    plt.show()
    
    plt.scatter(dopler((data[:, [0]] + FreqStart) * (10 ** 6), x_keckhelio(ra, dec, epoch=2000.0, jd=jd, tai=None, longitude=longitude, latitude=latitude, altitude=altitude, obs='irbene')),     data[:, [2]] * Systemtemperature9u * scale9U)
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend("9u")
    plt.show()
    '''
    sys.exit(0)