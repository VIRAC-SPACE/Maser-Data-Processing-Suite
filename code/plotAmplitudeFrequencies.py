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
    def __init__(self,  window, xdata, ydataU1, ydataU9, dataPoints, scan, Systemtemperature1u, Systemtemperature9u):
        
        #Window
        Frame.__init__(self)
        self.window = window
        self.Frame = frame(window,(1000,1000),background="gray")
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
        
        self.scan = scan
        self.Systemtemperature1u = Systemtemperature1u
        self.Systemtemperature9u = Systemtemperature9u
        
    def plotDataPoints (self):
        self.startDataPlotButton.destroy()
        self.createPolynomialButton = Button (self.Frame, text="Create Polynomial", command=self.plotPolynomial)
        self.createPolynomialButton.grid(row=0, column=2)
        
        self.points_1u = list()
        self.points_9u = list()
        
        #plt.suptitle("source " + scan["sourceName"].split(",")[0] + " scan " + str(scanNumber), fontsize=16)
        
        #u1
        self.fig1 = Figure(figsize=(6,6))
        self.graph1 = self.fig1.add_subplot(111)
        self.canvas1 = FigureCanvasTkAgg(self.fig1, master=self.Frame)
        self.canvas1.show()#nodzes ieprieksejos grafikus
        #self.fig3.clf()
        #self.fig4.clf()
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
        self.fig2 = Figure(figsize=(6,6))
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
        self.mSlider = Scale(self.Frame, from_= self.m, to = self.a-1, orient=HORIZONTAL, label="M", length=500, variable=self.m)
        self.mSlider.grid(row=2, column=0)
        self.m = self.mSlider.get()
        
        self.nSlider = Scale(self.Frame, from_ = self.b-1 , to = self.n, orient=HORIZONTAL, label="N", length=500, variable=self.n)
        self.nSlider.grid(row=2, column=1)
        self.nSlider.set(self.n)
        self.n = self.nSlider.get() 
    
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
        
        #nodzes ieprieksejos grafikus
        self.fig1.clf()
        self.fig2.clf()
        self.canvas1.get_tk_widget().delete("all")
        self.canvas2.get_tk_widget().delete("all")
        
        self.createPolynomialButton.destroy()
        self.plotLocalMaximumButton = Button (self.Frame, text="Create local maximum", command=self.plotLocalMaximum)
        self.plotLocalMaximumButton.grid(row=0, column=2)
        
        self.m = self.mSlider.get()
        self.n = self.nSlider.get()
        
        self.mSlider.destroy()
        self.nSlider.destroy()
        
        self.xarray_u1 = self.xarray
        self.xarray_u9 = self.xarray
        
        print "pirms dzesanas ", self.xarray_u9.shape[0]
        
        for p_u1 in self.points_1u:
            self.xarray_u1 = np.delete(self.xarray_u1, self.xarray_u1[self.xarray_u1 == p_u1[0]])
            self.z1 = np.delete(self.z1, self.z1[self.z1 == p_u1[1]])
            
        for p_u9 in self.points_9u:
            self.xarray_u9 = np.delete(self.xarray_u9, self.xarray_u9[self.xarray_u9 == p_u9[0]])
            self.z2 = np.delete(self.z2, self.z2[self.z2 == p_u9[1]])
        
            
        self.dataPoints_u1 = self.xarray_u1.shape[0]
        self.dataPoints_u9 = self.xarray_u9.shape[0]
        
        middle_u1 = int(self.dataPoints_u1 / 2) #vidusunks u1
        middle_u9 = int(self.dataPoints_u9 / 2) #vidusunks u9
        
        self.a_u1 = int(middle_u1 * 0.88)
        self.a_u9 = int(middle_u9 * 0.88)
        
        self.b_u1 = int(middle_u1 * 1.075)
        self.b_u9 = int(middle_u9 * 1.075)
        
        print "pec dzesanas ", self.xarray_u9.shape[0]
        
        timeStr = self.scan['startTime'].replace(":", " ")
        dateStrList = self.scan['dates'].split()
        dateStrList[1] = strptime(dateStrList[1],'%b').tm_mon
        dateStr = str(dateStrList[2]) + " " + str(dateStrList[1]) + " " + str(dateStrList[0])
        RaStr = " ".join(self.scan["Ra"])
        DecStr = " ".join(self.scan["Dec"])
        dopsetPar= dateStr + " " + timeStr + " " + RaStr + " " + DecStr
        
        os.system("code/dopsetpy_v1.5 " + dopsetPar)
    
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
          
        FreqStart = self.scan["FreqStart"] 
        
        #Parveido frekvenci par atrumu
        self.x_u1 = dopler((self.xarray_u1 + FreqStart) * (10 ** 6), VelTotal)
        self.x_u9 = dopler((self.xarray_u9 + FreqStart) * (10 ** 6), VelTotal)
        
        # Fit the data using a Chebyshev astro py
        ceb = Chebyshev1D(9, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
        fit_ceb = fitting.LevMarLSQFitter()
        
        ### u1
        self.ceb_1 = fit_ceb(ceb, np.append(self.x_u1[self.m:self.a_u1], self.x_u1[self.b_u1:self.n]),  np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]))
       
        ### u9
        self.ceb_2 = fit_ceb(ceb, np.append(self.x_u9[self.m:self.a_u9], self.x_u9[self.b_u9:self.n]),  np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]))
        
        #u1 plot
        self.fig3 = Figure(figsize=(6,6))
        self.graph3 = self.fig3.add_subplot(111)
        self.canvas1 = FigureCanvasTkAgg(self.fig3, master=self.Frame)
        self.canvas1.show()
        self.fig3.set_canvas(self.canvas1)
        self.canvas1.get_tk_widget().grid(row=1, column=0)
        self.graph3.plot(np.append(self.x_u1[self.m:self.a_u1], self.x_u1[self.b_u1:self.n]), np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]), 'ko', label='Data Points', markersize=1)
        self.graph3.plot(self.x_u1[self.m:self.n], self.ceb_1(self.x_u1[self.m:self.n]), 'r', label='Chebyshev polynomial')  
        self.graph3.grid(True)
        self.graph3.set_xlabel('Velocity (km sec$^{-1}$)')
        self.graph3.set_ylabel ('Flux density (Jy)')
        self.graph3.legend(loc=2)
        
        #u9 plot
        self.fig4 = Figure(figsize=(6,6))
        self.graph4 = self.fig4.add_subplot(111)
        self.canvas2 = FigureCanvasTkAgg(self.fig4, master=self.Frame)
        self.canvas2.show()
        self.fig4.set_canvas(self.canvas2)
        self.canvas2.get_tk_widget().grid(row=1, column=1)
        self.graph4.plot(np.append(self.x_u9[self.m:self.a_u9], self.x_u9[self.b_u9:self.n]), np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]), 'ko', label='Data Points', markersize=1)
        self.graph4.plot(self.x_u9[self.m:self.n], self.ceb_2(self.x_u9[self.m:self.n]), 'r', label='Chebyshev polynomial')
        self.graph4.grid(True)
        self.graph4.set_xlabel('Velocity (km sec$^{-1}$)')
        self.graph4.set_ylabel ('Flux density (Jy)')
        self.graph4.legend(loc=2)
        
        
    def plotLocalMaximum(self):
        #nodzes ieprieksejos grafikus
        #self.fig3.clf()
        #self.fig4.clf()
        self.canvas1.get_tk_widget().delete("all")
        self.canvas2.get_tk_widget().delete("all")
        
        self.plotLocalMaximumButton.destroy()
        self.monitoringButton = Button (self.Frame, text="Add point to monitoring", command=self.createResult)
        self.monitoringButton.grid(row=0, column=2)
        
        thres=0.1
    
        y1values = self.z1[self.m:self.n] - self.ceb_1(self.x_u1[self.m:self.n])
        y2values = self.z2[self.m:self.n] - self.ceb_2(self.x_u9[self.m:self.n])
        
        #indexsu apreikinasana
        indexes_for_ceb = peakutils.indexes(y1values, thres=thres, min_dist=10)
        indexes_for_ceb2 = peakutils.indexes(y2values, thres=thres, min_dist=10)
        
        #u1 plot
        self.fig5 = Figure(figsize=(6,6))
        self.graph5 = self.fig5.add_subplot(111)
        self.canvas3 = FigureCanvasTkAgg(self.fig5, master=self.Frame)
        self.canvas3.show()
        self.fig5.set_canvas(self.canvas3)
        self.canvas3.get_tk_widget().grid(row=1, column=0)
        self.graph5.plot(self.x_u1[self.m:self.n], y1values, 'b', label='Signal - polynomial', markersize=1)  
        self.graph5.plot(self.x_u1[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb], 'dr', label="Local Maximums for signal", markersize=2)
        
        ax = self.fig5.add_subplot(111)
        for xy in zip(self.x_u1[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb]):                        
            ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
        
        self.graph5.grid(True)
        self.graph5.set_xlabel('Velocity (km sec$^{-1}$)')
        self.graph5.set_ylabel ('Flux density (Jy)')
        self.graph5.legend(loc=2)
        
        #u9 plot
        self.fig6 = Figure(figsize=(6,6))
        self.graph6 = self.fig6.add_subplot(111)
        self.canvas4 = FigureCanvasTkAgg(self.fig6, master=self.Frame)
        self.canvas4.show()
        self.fig6.set_canvas(self.canvas4)
        self.canvas4.get_tk_widget().grid(row=1, column=1)
        self.graph6.plot(self.x_u9[self.m:self.n], y2values, 'b', label='Signal - polynomial', markersize=1)  
        self.graph6.plot(self.x_u9[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2], 'dr', label="Local Maximums for signal", markersize=2)
        
        ax = self.fig6.add_subplot(111)
        for xy in zip(self.x_u9[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2]):                        
            ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
        
        self.graph6.grid(True)
        self.graph6.set_xlabel('Velocity (km sec$^{-1}$)')
        self.graph6.set_ylabel ('Flux density (Jy)')
        self.graph6.legend(loc=2)
        
        #mid plot
        avg_x = (self.x_u1[self.m:self.n] + self.x_u9[self.m:self.n]) / 2
        avg_y = (y1values + y2values) / 2
        indexes_for_avg = peakutils.indexes(avg_y, thres=thres, min_dist=10)
        
        self.fig7 = Figure(figsize=(6,6))
        self.graph7 = self.fig7.add_subplot(111)
        self.canvas5 = FigureCanvasTkAgg(self.fig7, master=self.Frame)
        self.canvas5.show()
        self.fig7.set_canvas(self.canvas5)
        self.canvas5.get_tk_widget().grid(row=1, column=2)
        self.graph7.plot(avg_x, avg_y, 'b', label='Signal - polynomial', markersize=1)
        self.graph7.plot(avg_x[indexes_for_avg], avg_y[indexes_for_avg], 'dr', label="Local Maximums for signal", markersize=2)
        
        ax = self.fig7.add_subplot(111)
        for xy in zip(avg_x[indexes_for_avg],  avg_y[indexes_for_avg]):                        
            ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
        
        self.graph7.grid(True)
        self.graph7.set_xlabel('Velocity (km sec$^{-1}$)')
        self.graph7.set_ylabel ('Flux density (Jy)')
        self.graph7.legend(loc=2)
        
    def createResult(self):
        pass
    
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
    
    return (Systemtemperature1u, Systemtemperature9u, location, scan)

def main(logFileName, corData):
    
    #get Data and Logs
    Systemtemperature1u, Systemtemperature9u, location, scan = getLogs(logFileName, corData)
    xdata, y1data, y2data, dataPoints = getData(corData, location, Systemtemperature1u, Systemtemperature9u)
    
    #Create App
    window = tk.Tk() 
    ploting = MaserPlot(window, xdata, y1data, y2data, dataPoints, scan, Systemtemperature1u, Systemtemperature9u)
    ploting.mainloop()
    sys.exit(0)

if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    
    logFileName = sys.argv[1]
    corData = sys.argv[2]
    
    main(logFileName, corData)