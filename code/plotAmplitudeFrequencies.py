#! /usr/bin/python

import os, sys
import numpy as np
from Tkinter import *
import Tkinter as tk
from time import strptime
import scipy.constants
from astropy.modeling import fitting
from astropy.modeling.polynomial import Chebyshev1D
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.interpolate import UnivariateSpline
import peakutils
import json

from ploting import Plot
from experimentsLogReader import ExperimentLogReader

calibrationScales = {"IRBENE":12, "IRBENE16":26}

def usage():
    print ('Usage: ' + sys.argv[0] + ' log file' + 'data file')

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def sortArrayByIndexies(orginalArray, sortIndexies):
    if len(orginalArray) != len(sortIndexies):
        raise Exception("Arrays must be same size")
    else:
        sortedArray = [0]*len(orginalArray)
        for i in range(0,  len(sortIndexies)):
            sortedArray[int(sortIndexies[i])] = orginalArray[i]      
    return sortedArray

def calibration(calibrationScale, Tsys):
    return calibrationScale*Tsys
    
def dopler(ObservedFrequency, velocityReceiver, f0):
    c = scipy.constants.speed_of_light
    #f0 = 6668519200 # Hz 
    velocitySoure = (-((ObservedFrequency/f0)-1)*c + (velocityReceiver * 1000))/1000
    return velocitySoure

def FWHM(x, y, constant):
    spline = UnivariateSpline(x, y-np.max(y)/2, k=3, s=20)
    spline.set_smoothing_factor(0.5)
    root1 = spline.roots()[0] - constant
    root2 = spline.roots()[-1] + constant
    index_1 =  (np.abs(x-root1)).argmin()
    index_2 =  (np.abs(x-root2)).argmin()
    return (index_1, index_2)
    
def is_outlier(points, thresh=4.5):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score < thresh
    
def frame(parent, size, sides, **options):
    Width=size[0]
    Height=size[1]
    f=Frame(parent, width=Width,height=Height,**options)
    f.pack(side = sides)
    return (f)
    
class MaserPlot(Frame):
    def __init__(self,  window, xdata, ydataU1, ydataU9, dataPoints, Systemtemperature1u, Systemtemperature9u, expername, source, location, scan, scanNumber):
        
        #Data init
        Frame.__init__(self)
        self.window = window
        self.xdata = xdata
        self.ydataU1= ydataU1
        self.ydataU9 = ydataU9
        self.location = location
        self.Systemtemperature1u = Systemtemperature1u
        self.Systemtemperature9u = Systemtemperature9u
        self.source = source
        self.expername = expername
        self.scan = scan
        self.scanNumber = scanNumber
        self.dataPoints = dataPoints
        
        #Making sure that data is numpy array
        self.xarray = np.zeros(self.dataPoints)
        self.y1array = np.zeros(self.dataPoints)
        self.y2array = np.zeros(self.dataPoints)
        
        for i in range(0,dataPoints):
            self.xarray[i] = self.xdata[i]
        
        for j in range(0,dataPoints):
            self.y1array[j] = self.ydataU1[j]
        
        for k in range(0,dataPoints):
            self.y2array[k] = self.ydataU9[k]
            
        self.startWindow()
    
    def startWindow(self):
        #default constants
        self.FWHMconstant = 0.3
        self.polynomialOrder = 9
        self.f0 = 6668519200 
        self.calibrationScale = calibrationScales[self.location]
        self.FreqStart = self.scan["FreqStart"]
        
        #start window frame
        self.plotFrame = frame(self.window,(1000,1000), None, background = "gray")
        self.infoFrame = frame(self.window,(1000,1000), None, background = "gray")
        self.startDataPlotButton = Button (self.plotFrame, text="Plot data points", command=self.plotDataPoints)
        self.startChangeData = Button (self.plotFrame, text="Change Data", command=self.changeData)
        self.startDataPlotButton.pack(fill=BOTH)
        self.startChangeData.pack(side=BOTTOM, fill=BOTH)
            
        #infoFrame
        self.window.title("Info")
        infoPanelLabelsText = ["Experiment name: " + self.expername, "Scan number: " + self.scanNumber, "Source: " + self.source, "Station: " + self.location, "Date: " + self.scan["dates"], "Start time: " + self.scan["startTime"], "Stop time: " + self.scan["stopTime"], "System temperature 1u: " + str(self.Systemtemperature1u), "System temperature 9u: " + str(self.Systemtemperature9u), "Frequency Start: " + str(self.scan["FreqStart"]), "f0", "Calibration scale", "FWHM constant", "Polynomial order"]
        infoPanelEntryText = [{"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"defaultValue":self.Systemtemperature1u,"addEntry":True}, {"defaultValue":self.Systemtemperature9u,"addEntry":True}, {"defaultValue":self.scan["FreqStart"],"addEntry":True}, {"defaultValue":self.f0, "addEntry":True}, {"defaultValue":str(self.calibrationScale), "addEntry":True}, {"defaultValue":str(self.FWHMconstant), "addEntry":True}, {"defaultValue":str(self.polynomialOrder), "addEntry":True}]
        
        for i in range(0, len( infoPanelLabelsText)): 
            self.infoLabel = Label(self.infoFrame, text=infoPanelLabelsText[i], anchor=W, justify=LEFT)
            self.infoLabel.pack(side=TOP, fill=BOTH)
            
            if  infoPanelEntryText[i]["addEntry"]:
                self.infoInputField = Entry(self.infoFrame)
                self.infoInputField.insert(0, str(infoPanelEntryText[i]["defaultValue"]))
                self.infoInputField.pack(fill=BOTH)
                
    def changeData(self):
        childs = self.infoFrame.winfo_children()
        newValues = list()
        
        for child in childs:
            if child.winfo_class() == "Entry":
                newValues.append(child.get())
                
        self.Systemtemperature1u = float(newValues[0])
        self.Systemtemperature9u = float(newValues[1])
        self.FreqStart = float(newValues[2])
        self.f0 = float(newValues[3])
        self.calibrationScale = float(newValues[4])
        self.FWHMconstant = float(newValues[5])
        self.polynomialOrder = int(newValues[6])
    
    def calibration(self):
        self.y1array = self.y1array * calibration(self.calibrationScale, self.Systemtemperature1u)
        self.y2array = self.y2array * calibration(self.calibrationScale, self.Systemtemperature9u)
            
        g1 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
        g2 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
    
        self.z1 = convolve(self.y1array, g1, boundary='extend')
        self.z2 = convolve(self.y2array, g2, boundary='extend')
        
        self.a, self.b = FWHM(self.xarray, (self.z1 + self.z2)/2, self.FWHMconstant)
        self.m = 0
        self.n = self.dataPoints
    
    def back(self):
        self.plot_1.removePolt()
        self.plot_2.removePolt()
        self.masterFrame.destroy()
        self.createPolynomialButton.destroy()
        self.backButton.destroy()
        self.mSlider.destroy()
        self.nSlider.destroy()
        self.startWindow()
          
    def plotDataPoints (self):
        self.calibration()
        self.masterFrame = frame(self.window,(1000,1000), LEFT, background = "gray")
        self.window.title("Data points for " + self.expername + " scan " +  self.scanNumber +  " for Source " + self.source)
        self.startChangeData.destroy()
        self.startDataPlotButton.destroy()
        self.infoFrame.destroy()
        self.createPolynomialButton = Button (self.plotFrame, text="Create Polynomial", command=self.plotPolynomial)
        self.createPolynomialButton.pack(side=TOP)
        self.backButton = Button (self.plotFrame, text="back", command=self.back)
        self.backButton.pack(side=TOP)
        
        self.points_1u = list()
        self.points_9u = list()
        
        #u1
        self.plot_1 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_1.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_1.plot(self.xarray, self.z1, 'ko', 'Data Points', 1, 5)
        self.plot_1.addPickEvent(self.onpickU1)
        self.plot_1.addSecondAss("x", "Data points", 0, self.dataPoints + 512, 1024)
        
        #u9
        self.plot_2 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_2.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_2.plot(self.xarray, self.z2, 'ko', 'Data Points', 1, 5)
        self.plot_2.addPickEvent(self.onpickU9)
        self.plot_2.addSecondAss("x", "Data points", 0, self.dataPoints + 512, 1024)
        
        #sliders
        self.mSlider = Scale(self.plotFrame, from_= self.m, to = self.a-1, orient=HORIZONTAL, label="M", length=500, variable=self.m)
        self.mSlider.pack(side=BOTTOM)
        self.m = self.mSlider.get()
        
        self.nSlider = Scale(self.plotFrame, from_ = self.b-1 , to = self.n, orient=HORIZONTAL, label="N", length=500, variable=self.n)
        self.nSlider.pack(side=BOTTOM)
        self.nSlider.set(self.n)
        self.n = self.nSlider.get() 
    
    def onpickU1(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_1.plot(p[0][0], p[0][1], 'ro', None, 1, 5)
        self.points_1u.append(p[0])
        self.plot_1.canvasShow()
        
    def onpickU9(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_2.plot(p[0][0], p[0][1], 'ro', None, 1, 5)
        self.points_9u.append(p[0])
        self.plot_2.canvasShow()
        
    def onpick_maxU1(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.maxu1_index.append(ind[0])
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_5.plot(p[0][0], p[0][1], 'gd', None, 2, 5)
        if  self.maxU1.count(p[0]) == 0:
            self.maxU1.append(p[0])
        self.plot_5.canvasShow()
        
    def onpick_maxU9(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.maxu9_index.append(ind[0])
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_6.plot(p[0][0], p[0][1], 'gd', None, 2, 5)
        if  self.maxU9.count(p[0]) == 0:
            self.maxU9.append(p[0])
        self.plot_6.canvasShow()
        
    def onpick_maxAVG(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.maxavg_index.append(ind[0])
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_7.plot(p[0][0], p[0][1], 'gd', None, 2, 5)
        if  self.avgMax.count(p[0]) == 0:
            self.avgMax.append(p[0])
        self.plot_7.canvasShow()
    
    def plotPolynomial(self):
        self.window.title("Polynomial " + self.expername + " scan " +  self.scanNumber +  " for Source " + self.source)
        #nodzes ieprieksejos grafikus
        self.plot_1.removePolt()
        self.plot_2.removePolt()
        
        self.plot_1.removePickEvent()
        self.plot_2.removePickEvent()
        
        self.createPolynomialButton.destroy()
        self.plotLocalMaximumButton = Button (self.plotFrame, text="Create local maximum", command=self.plotLocalMaximum)
        self.plotLocalMaximumButton.pack()
        
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
        
        self.a_u1, self.b_u1 = FWHM(self.xarray_u1, self.z1, self.FWHMconstant)
        self.a_u9, self.b_u9 = FWHM(self.xarray_u9, self.z2, self.FWHMconstant)
         
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
          
        #Parveido frekvenci par atrumu
        self.x_u1 = dopler((self.xarray_u1 + self.FreqStart) * (10 ** 6), VelTotal, self.f0)
        self.x_u9 = dopler((self.xarray_u9 + self.FreqStart) * (10 ** 6), VelTotal, self.f0)
        
        # Fit the data using a Chebyshev astro py
        ceb = Chebyshev1D(self.polynomialOrder, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
        fit_ceb = fitting.LevMarLSQFitter()
        
        ### u1
        self.ceb_1 = fit_ceb(ceb, np.append(self.x_u1[self.m:self.a_u1], self.x_u1[self.b_u1:self.n]),  np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]))
       
        ### u9
        self.ceb_2 = fit_ceb(ceb, np.append(self.x_u9[self.m:self.a_u9], self.x_u9[self.b_u9:self.n]),  np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]))
        
        #u1 plot
        self.plot_3 = Plot(6,6, self.window, self.plotFrame)
        self.plot_3.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_3.plot(np.append(self.x_u1[self.m:self.a_u1], self.x_u1[self.b_u1:self.n]), np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]), 'ko', 'Data Points', 1, None)
        self.plot_3.plot(self.x_u1[self.m:self.n], self.ceb_1(self.x_u1[self.m:self.n]), 'r', 'Chebyshev polynomial',  1, None)
        
        #u9 plot
        self.plot_4 = Plot(6,6, self.window, self.plotFrame)
        self.plot_4.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_4.plot(np.append(self.x_u9[self.m:self.a_u9], self.x_u9[self.b_u9:self.n]), np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]), 'ko', 'Data Points', 1, None)
        self.plot_4.plot(self.x_u9[self.m:self.n], self.ceb_2(self.x_u9[self.m:self.n]), 'r', 'Chebyshev polynomial', 1, None)
        
    def plotLocalMaximum(self):
        self.window.title("Local maximums " + self.expername + " scan " +  self.scanNumber +  " for Source " + self.source)
        #nodzes ieprieksejos grafikus
        self.plot_3.removePolt()
        self.plot_4.removePolt()
        
        self.plotLocalMaximumButton.destroy()
        self.monitoringButton = Button (self.plotFrame, text="Add points to monitoring", command=self.createResult)
        self.monitoringButton.pack()
        
        thres=0.1
    
        y1values = self.z1[self.m:self.n] - self.ceb_1(self.x_u1[self.m:self.n])
        y2values = self.z2[self.m:self.n] - self.ceb_2(self.x_u9[self.m:self.n])
        
        #indexsu apreikinasana
        indexes_for_ceb = peakutils.indexes(y1values, thres=thres, min_dist=10)
        indexes_for_ceb2 = peakutils.indexes(y2values, thres=thres, min_dist=10)
        
        #u1
        self.plot_5 = Plot(6,6, self.window, self.plotFrame)
        self.plot_5.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_5.plot(self.x_u1[self.m:self.n], y1values, 'b', 'Signal - polynomial', 1, None)
        self.plot_5.plot(self.x_u1[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb], 'dr', "Local Maximums for signal", 2, 5)
        self.plot_5.addPickEvent(self.onpick_maxU1)
        self.plot_5.annotation(self.x_u1[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb])
        
        #u9
        self.plot_6 = Plot(6,6, self.window, self.plotFrame)
        self.plot_6.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_6.plot(self.x_u9[self.m:self.n], y1values, 'b', 'Signal - polynomial', 1, None)
        self.plot_6.plot(self.x_u9[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2], 'dr', "Local Maximums for signal", 2, 5)
        self.plot_6.addPickEvent(self.onpick_maxU9)
        self.plot_6.annotation(self.x_u9[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2])
        
        #mid plot
        avg_x = (self.x_u1[self.m:self.n] + self.x_u9[self.m:self.n]) / 2
        avg_y = (y1values + y2values) / 2
        indexes_for_avg = peakutils.indexes(avg_y, thres=thres, min_dist=10)
        
        self.plot_7 = Plot(6,6, self.window, self.plotFrame)
        self.plot_7.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Average Polarization")
        self.plot_7.plot(avg_x, avg_y, 'b', 'Signal - polynomial', 1, None)
        self.plot_7.plot(avg_x[indexes_for_avg], avg_y[indexes_for_avg], 'dr', "Local Maximums for signal", 2, 5)
        self.plot_7.addPickEvent(self.onpick_maxAVG)
        self.plot_7.annotation(avg_x[indexes_for_avg],  avg_y[indexes_for_avg])
        
        self.maxU1 = list()
        self.maxU9 = list()
        self.avgMax = list()
        self.maxu1_index = list()
        self.maxu9_index = list()
        self.maxavg_index = list()
        
    def createResult(self):
        #remove graph
        self.plot_5.removePolt()
        self.plot_6.removePolt()
        self.plot_7.removePolt()
        
        self.monitoringButton.destroy()
    
        endLabel = Label(master=self.plotFrame, text="Result file creating in progress!")
        endLabel.pack()
        
        max_x_U1 = list()
        max_x_U9 = list()
        max_x_avg = list()
        max_y_U1 = list()
        max_y_U9 = list()
        max_y_avg = list()
        
        for i in range(0, len(self.maxU1)):
            max_x_U1.append(self.maxU1[i][0])
            max_y_U1.append(self.maxU1[i][1])
            
        for i in range(0, len(self.maxU9)):
            max_x_U9.append(self.maxU9[i][0])
            max_y_U9.append(self.maxU9[i][1])
        
        for i in range(0, len(self.avgMax)):
            max_x_avg.append(self.avgMax[i][0])
            max_y_avg.append(self.avgMax[i][1])
        
        #sorting array to match index    
        max_x_U1 = sortArrayByIndexies(max_x_U1, self.maxu1_index)
        max_y_U1 = sortArrayByIndexies(max_y_U1, self.maxu1_index)
        max_x_U9 = sortArrayByIndexies(max_x_U9, self.maxu9_index)
        max_y_U9 = sortArrayByIndexies(max_y_U9, self.maxu9_index)
        max_x_avg = sortArrayByIndexies(max_x_avg, self.maxavg_index)
        max_y_avg = sortArrayByIndexies(max_y_avg, self.maxavg_index)
            
        resultDir = "results/"
        resultFileName = self.source + ".json"
    
        if os.path.isfile(resultDir + resultFileName):
            pass
        else:
            os.system("touch " + resultDir +  resultFileName)
            
            resultFile = open (resultDir +  resultFileName, "w")
            resultFile.write("{ \n" + "\n}")
            resultFile.close()
        
        with open(resultDir + resultFileName) as result_data:    
            result = json.load(result_data)
        
        if self.expername not in result:
            result[self.expername] = dict()
            if self.scanNumber not in result[self.expername]:
                result[self.expername][self.scanNumber] = dict()
                
        if self.scanNumber not in result[self.expername]:
                result[self.expername][self.scanNumber] = dict()
                   
        result[self.expername][self.scanNumber]["startTime"] = self.scan["startTime"]
        result[self.expername][self.scanNumber]["stopTime"] = self.scan["stopTime"]
        result[self.expername][self.scanNumber]["location"] = self.location
        result[self.expername][self.scanNumber]["Date"] = self.scan["dates"]
                
        result[self.expername][self.scanNumber]["velocity_for_polarizationU1"] = max_x_U1
        result[self.expername][self.scanNumber]["velocity_for_polarizationU9"] = max_x_U9
        result[self.expername][self.scanNumber]["velocity_for_polarizationAVG"] = max_x_avg 
                
        result[self.expername][self.scanNumber]["amplitude_for_polarizationU1"] = max_y_U1
        result[self.expername][self.scanNumber]["amplitude_for_polarizationU9"] = max_y_U9
        result[self.expername][self.scanNumber]["amplitude_for_polarizationAVG"] =  max_y_avg
        
        #result[self.expername][self.scanNumber]["index_for_polarizationU1"] =  self.maxu1_index
        #result[self.expername][self.scanNumber]["index_for_polarizationU9"] =  self.maxu9_index
        #result[self.expername][self.scanNumber]["index_for_polarizationAVG"] =  self.maxavg_index
        
        resultFile = open (resultDir +  resultFileName, "w")
        resultFile.write(json.dumps(result, indent=2))
        resultFile.close() 
        
        self. _quit()
        
    def _quit(self):
        self.plotFrame.destroy()
        self.window.destroy()
        
def getData(dataFileName):
    data = np.fromfile(dataFileName, dtype="float64", count=-1, sep=" ") .reshape((file_len(dataFileName),5))
    data = np.delete(data, (0), axis=0) #izdzes masiva primo elementu
    
    dataPoints = data.shape[0]    
    
    xdata = data[:, [0]]
    y1data = data[:, [1]] 
    y2data = data[:, [2]]
    
    return (xdata, y1data, y2data, dataPoints)

def getLogs(logfileName, dataFileName): 
    logs  = ExperimentLogReader("logs/" + logfileName, "prettyLogs/").getLgs()
    scanNumber = dataFileName.split(".")[0].split("_")[1][1:len(dataFileName)]
    scan = logs[scanNumber]
    
    Systemtemperature1u = float(scan["Systemtemperature"][0])
    Systemtemperature9u = float(scan["Systemtemperature"][1])
    
    location = logs["location"]
    source = scan["source"]
    
    return (Systemtemperature1u, Systemtemperature9u, location, source, scan, scanNumber)

def main(logFileName, corData):
    
    #get Data and Logs
    Systemtemperature1u, Systemtemperature9u, location, source, scan, scanNumber = getLogs(logFileName, corData)
    xdata, y1data, y2data, dataPoints = getData(corData)
    expername = logFileName.split(".")[0][:-2]
    
    #Create App
    window = tk.Tk() 
    ploting = MaserPlot(window, xdata, y1data, y2data, dataPoints, Systemtemperature1u, Systemtemperature9u, expername, source, location, scan, scanNumber)
    ploting.mainloop()
    
    sys.exit(0)

if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
        
    logFileName = sys.argv[1]
    corData = sys.argv[2]
    
    main(logFileName, corData)