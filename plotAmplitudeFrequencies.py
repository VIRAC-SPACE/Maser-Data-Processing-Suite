#! /usr/bin/python
from __future__ import division

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from time import strptime
import scipy.constants
from scipy import interpolate
from scipy.interpolate import interp1d

from experimentsLogReader import ExperimentLogReader

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def usage():
    print ('Usage: '+sys.argv[0]+' log file' + 'data file')

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
    
    #slito punktu izdzesana
    outliersMask = is_outlier(data[:, [0]])
    data = data[outliersMask]
    
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
    
    os.system("./dopsetpy_v1.5 "+dopsetPar)
    
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
    x = dopler((data[:, [0]] + FreqStart) * (10 ** 6), VelTotal)
   
    y = data[:, [1]] * calibration(location, Systemtemperature1u)
    x = np.arange(0,y.size )
    print x.size, y.size
    f = interp1d([1,2,3],[1,4,6])
    
    plt.plot(dopler((data[:, [0]] + FreqStart) * (10 ** 6), VelTotal), data[:, [1]] * calibration(location, Systemtemperature1u), 'ro')
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend("1u")
    plt.show()
    
    plt.plot(dopler((data[:, [0]] + FreqStart) * (10 ** 6), VelTotal), data[:, [2]] * calibration(location, Systemtemperature1u), 'ro')
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend("9u")
    plt.show()
    
    sys.exit(0)