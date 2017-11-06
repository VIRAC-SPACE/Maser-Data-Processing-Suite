#! /usr/bin/python

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from PyAstronomy import pyasl
from time import strptime
import scipy.constants
import jdcal
from jdcal import gcal2jd
import datetime

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
    
if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    
    logs  = ExperimentLogReader(sys.argv[1]).getLgs()
    
    data = np.fromfile(sys.argv[2], dtype="float64", count=-1, sep=" ") .reshape((file_len(sys.argv[2]),5))
    data = np.delete(data, (0), axis=0)
    
    scanNumber = sys.argv[2].split(".")[0].split("_")[1][1:len(sys.argv[2])]
    scan = logs[scanNumber]
    Systemtemperature1u = float(scan["Systemtemperature"][0])
    Systemtemperature9u = float(scan["Systemtemperature"][1])
    
    scale1U = 1
    scale9U = 1
    
    location = logs["location"]
    
    if location == "IRBENE":
        scale1U = 12
        scale9U = 12
        
    elif location == "IRBENE16":
        scale1U = 26
        scale9U = 26
    
    plt.plot(data[:, [0]], data[:, [1]] * Systemtemperature1u * scale1U)
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("1u")
    plt.show()
    
    plt.plot(data[:, [0]], data[:, [2]] * Systemtemperature9u * scale9U)
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("9u")
    plt.show()
     
    longitude = 21.8605
    latitude = 57.5546
    altitude = 87.30
    
    ra = float(scan["Ra"][0]) + float(scan["Ra"][1])/60 + float(scan["Ra"][2])/3600
    dec = float(scan["Dec"][0]) + float(scan["Dec"][1])/60 + float(scan["Dec"][2])/3600
    
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
    jd = JDN + ((hours-12.0)/24.0) + (minutes/1440.0) + (seconds/86400.0)
    
    corr, hjd = pyasl.helcorr(longitude, latitude, altitude, \
            ra, dec, jd, debug=True)
    
    print("ra ", scan["Ra"])
    print("dec ", scan["Dec"])
    print("Barycentric correction [km/s]: ", corr)
    print("Heliocentric Julian day: ", hjd)
    
    #FreqStart = float(scan["FreqStart"])
    FreqStart = 6668.83
    print("FreqStart", FreqStart)
    plt.plot(dopler((data[:, [0]] + FreqStart) * (10 ** 6), corr),      data[:, [1]] * Systemtemperature1u * scale1U)
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend("1u")
    plt.show()
    
    plt.plot(dopler((data[:, [0]] + FreqStart) * (10 ** 6), corr),     data[:, [2]] * Systemtemperature9u * scale9U)
    plt.grid(True)
    plt.xlabel('velocity')
    plt.legend("9u")
    plt.show()
    
    sys.exit(0)