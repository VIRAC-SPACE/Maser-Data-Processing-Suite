#! /usr/bin/python
from __future__ import division

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib.widgets import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from Tkinter import *
import Tkinter as tk
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

    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score < thresh

    
def frame(parent, sides,**options):
    Width=sides[0]
    Height=sides[1]
    f=Frame(width=Width,height=Height,**options)
    f.grid(row=0, column=0)
    return (f)
    
class MaserPlot(Frame):
    def __init__(self,  window, xdata, ydataU1, ydataU9, dataPoints):
        Frame.__init__(self)
        self.window = window
        self.Frame = frame(window,[1000,1000,1000,1000],background="gray")
        self.startDataPlotButton = Button (self.Frame, text="Plot data points", command=self.plotDataPoints)
        self.startDataPlotButton.grid(row=0, column=0)
        
        self.xdata = xdata
        self.ydataU1= ydataU1
        self.ydataU9 = ydataU9
        
        self.dataPoints = dataPoints
        middle = int(self.dataPoints/2) #vidusunks
        self.a = int(middle*0.88)
        self.b = int(middle*1.075)
        self.m = 0
        self.n = self.dataPoints
        
        y1array = np.zeros(dataPoints)
        y2array = np.zeros(dataPoints)
        
        self.xarray = np.zeros(dataPoints)
    
        for i in range(0,dataPoints):
            self.xarray[i] = self.xdata[i]
        
        for j in range(0,dataPoints):
            y1array[j] = self.ydataU1[j]
        
        for k in range(0,dataPoints):
            y2array[k] = self.ydataU9[k]
            
        g1 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
        g2 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
    
        self.z1 = convolve(y1array, g1, boundary='extend')
        self.z2 = convolve(y2array, g2, boundary='extend')
        
    def plotDataPoints (self):
        self.startDataPlotButton.destroy()
        self.createPolynomialButton = Button (self.Frame, text="Create Polynomial", command=self.plotPolynomial)
        self.createPolynomialButton.grid(row=0, column=2)
        
        self.points_1u = list()
        self.points_9u = list()
        
        plt.suptitle("source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber), fontsize=16)
        
        #u1
        self.fig1 = Figure(figsize=(5,5))
        self.graph1 = self.fig1.add_subplot(111)
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.Frame)
        self.canvas1.show()
        self.fig1.set_canvas(self.canvas1)
        self.canvas1.get_tk_widget().grid(row=1, column=0)
        self.canvas1.mpl_connect('pick_event', self.onpickU1)
        
        self.graph1.plot(self.xarray, self.z1, 'ko', label='Data Points', markersize=1,  picker=5)  
        self.graph1.grid(True)
        self.graph1.set_xlabel('Frequency Mhz')
        self.graph1.set_ylabel ('Flux density (Jy)')
        self.graph1.legend(loc=2)
        
        #pievieno papildus asi data punktiem
        self.second_x_ass = self.graph1.twiny()
        self.second_x_ass.set_xlabel("Data points")
        self.graph1.tick_params(axis="x")
        self.second_x_ass.set_xticks(range(0, self.dataPoints + 512, 1024))
    
        self.graph1.set_title("1u Polarization",  y=1.08) 
        
        #u9
        self.fig2 = Figure(figsize=(5,5))
        self.graph2 = self.fig2.add_subplot(111)
        self.canvas2 = FigureCanvasTkAgg(self.fig2, master=self.Frame)
        self.canvas2.show()
        self.fig2.set_canvas(self.canvas2)
        self.canvas2.get_tk_widget().grid(row=1, column=1)
        self.canvas2.mpl_connect('pick_event', self.onpickU9)
        
        self.graph2.plot(self.xarray, self.z2, 'ko', label='Data Points', markersize=1, picker=5)  
        self.graph2.grid(True)
        self.graph2.set_xlabel('Frequency Mhz')
        self.graph2.set_ylabel ('Flux density (Jy)')
        self.graph2.legend(loc=2)
        
        #pievieno papildus asi data punktiem
        self.second_x_ass_2 = self.graph2.twiny()
        self.second_x_ass_2.set_xlabel("Data points")
        self.graph2.tick_params(axis="x")
        self.second_x_ass_2.set_xticks(range(0, self.dataPoints + 512, 1024))
    
        self.graph2.set_title("9u Polarization",  y=1.08) 
        
        #sliders
        mSlider = Scale(self.Frame, from_= self.m, to = self.a-1, orient=HORIZONTAL, label="M", length=500, variable=self.m)
        mSlider.grid(row=2, column=0)
        self.m = mSlider.get() 
        
        nSlider = Scale(self.Frame, from_ = self.b-1 , to = self.n, orient=HORIZONTAL, label="N", length=500, variable=self.n)
        nSlider.grid(row=2, column=1)
        self.n = nSlider.get() 
    
    def onpickU1(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.graph1.plot(p[0][0], p[0][1], 'ro', markersize=1, picker=5)
        self.points_1u.append(p[0])
        self.canvas1.draw()
        
    def onpickU9(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.graph2.plot(p[0][0], p[0][1], 'ro', markersize=1, picker=5)
        self.points_9u.append(p[0])
        self.canvas2.draw()
    
    def plotPolynomial(self):
        self.canvas1.destroy()
        self.canvas2.destroy()
        
        self.xarray_u1 = self.xarray
        self.xarray_u9 = self.xarray
        
        for p_u1 in self.points_1u:
            self.xarray_u1.delete(p_u1[0])
            self.z1.delete(p_u1[1])
            
        for p_u9 in self.points_9u:
            self.xarray_u9.delete(p_u9[0])
            self.z2.delete(p_u9[1])
        
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
        
    def _quit(self):
        self.Frame.destroy()
        self.window.destroy()
        
def getData(dataFileName, location, Systemtemperature1u, Systemtemperature9u):
    data = np.fromfile(dataFileName, dtype="float64", count=-1, sep=" ") .reshape((file_len(dataFileName),5))
    data = np.delete(data, (0), axis=0) #izdzes masiva primo elementu
    
    dataPoints = data.shape[0]    
    
    xdata = data[:, [0]]
    y1data = data[:, [1]] * calibration(location, Systemtemperature1u)
    y2data = data[:, [2]] * calibration(location, Systemtemperature9u)  
    
    return (xdata, y1data, y2data, dataPoints)

def getLogs(logfileName, dataFileName): 
    logs  = ExperimentLogReader("logs/" + logfileName, "prettyLogs/").getLgs()
    scanNumber = dataFileName.split(".")[0].split("_")[1][1:len(dataFileName)]
    scan = logs[scanNumber]
    
    Systemtemperature1u = float(scan["Systemtemperature"][0])
    Systemtemperature9u = float(scan["Systemtemperature"][1])
    
    location = logs["location"]
    
    return (Systemtemperature1u, Systemtemperature9u, location)

def main(logFileName, corData):
    
    #get Data and Logs
    Systemtemperature1u, Systemtemperature9u, location = getLogs(logFileName, corData)
    xdata, y1data, y2data, dataPoints = getData(corData, location, Systemtemperature1u, Systemtemperature9u)
    
    #Create App
    window = tk.Tk() 
    ploting = MaserPlot(window, xdata, y1data, y2data, dataPoints)
    ploting.mainloop()
    sys.exit(0)

if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    
    logFileName = sys.argv[1]
    corData = sys.argv[2]
    
    main(logFileName, corData)
    
    '''
  
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
    
    ### u1
    ceb_1 = fit_ceb(ceb, np.append(xarray[m:a], xarray[b:n]),  np.append(z1[m:a],z1[b:n]))
   
    ### u9
    ceb_2 = fit_ceb(ceb, np.append(xarray[m:a], xarray[b:n]),  np.append(z2[m:a],z2[b:n]))
  
    #Polinom ploting 
      
    #1u
    fig = pylab.gcf()
    fig.canvas.set_window_title("Polynomial for experiment " +  experimentName)
    fig.set_size_inches(10.5, 10.5)
    plt.suptitle("source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber), fontsize=16)
    plt.subplot(121)
    plt.subplots_adjust(wspace = 0.35)
    plt.plot(xarray[m:n], ceb_1(xarray[m:n]), 'r', label='Chebyshev polynomial')
    plt.plot(np.append(xarray[m:a], xarray[b:n]), np.append(z1[m:a], z1[b:n]), 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("1u Polarization")  
    
    #9u
    plt.subplot(122)
    plt.plot(xarray[m:n], ceb_2(xarray[m:n]), 'r', label='Chebyshev polynomial')
    plt.plot(np.append(xarray[m:a], xarray[b:n]), np.append(z2[m:a], z2[b:n]), 'ko', label='Data Points')
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("9u Polarization")
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
    fig.canvas.set_window_title("Local Maximums for experiment " +  experimentName)
    fig.set_size_inches(10.5, 10.5)
    plt.suptitle("source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber), fontsize=16)
    plt.subplot(121)
    plt.subplots_adjust(wspace = 0.35)

    plt.plot(xarray[m:n], y1values, 'b', label='Signal - polynomial')  
    plt.plot(xarray[m:n][indexes_for_ceb], y1values[indexes_for_ceb], 'dr', label="Local Maximums for signal")
    
    ax = fig.add_subplot(121)
    for xy in zip(xarray[m:n][indexes_for_ceb], y1values[indexes_for_ceb]):                        
        ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
    
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("1u Polarization")
    
    #9u
    plt.subplot(122)
    plt.plot(xarray[m:n], y2values, 'b', label='Signal - polynomial')
    plt.plot(xarray[m:n][indexes_for_ceb2], y2values[indexes_for_ceb2], 'dr', label="Local Maximums for signal")
    
    xa = fig.add_subplot(122)
    for yx in zip(xarray[m:n][indexes_for_ceb2], y2values[indexes_for_ceb2]):                        
        xa.annotate('(%.2f, %.1f)' % yx, xy=yx, textcoords='data')
    
    plt.grid(True)
    plt.xlabel('Velocity (km sec$^{-1}$)')
    plt.ylabel('Flux density (Jy)')
    plt.legend(loc=2)
    plt.title("9u Polarization")
    plt.show()
    '''
    