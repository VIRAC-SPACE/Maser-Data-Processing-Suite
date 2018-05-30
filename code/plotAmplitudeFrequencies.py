#! /usr/bin/python
import os
import sys
import numpy as np
from tkinter import *
import tkinter as tk
from tkinter import font
from tkinter import messagebox
from time import strptime 
import scipy.constants
from astropy.modeling import fitting
from astropy.modeling.polynomial import Chebyshev1D
from astropy.convolution import Gaussian1DKernel, convolve
from scipy.interpolate import UnivariateSpline
import peakutils
import json
import argparse
import configparser

from ploting import Plot
from experimentsLogReader import ExperimentLogReader

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Plots maser plot. ''',
    epilog="""Maser Ploter.""")

    # Positional mandatory arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("logFile", help="Experiment log file name", type=str)
    parser.add_argument("datafile", help="Experiment correlation file name", type=str)

    # Optional arguments
    parser.add_argument("-s", "--single", help="Set RA, DEC, Epoch, Source name", nargs="*", type=str, default=[])

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')

    # Parse arguments
    args = parser.parse_args()

    return args
    
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

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
        
def frame(parent, size, sides, **options):
    Width=size[0]
    Height=size[1]
    f=tk.Frame(parent, width=Width, height=Height, background="light goldenrod", **options)
    f.pack(side = sides)
    return (f)

def indexies(array, value):
    indexs = list()
    for i in range(0, len(array)-1):
        if array[i] == value:
            indexs.append(i)
    return indexs
    
class MaserPlot(Frame):
    def __init__(self,  window, xdata, ydataU1, ydataU9, dataPoints, Systemtemperature1u, Systemtemperature9u, expername, source, location, scan, scanNumber, calibrationScales):
        
        #Data init
        Frame.__init__(self)
        self.font = font.Font(family="Times New Roman", size=20, weight=font.BOLD)
        self.font_2 = font.Font(family="Times New Roman", size=10)
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
        self.calibrationScales = calibrationScales
        
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
        self.state = 0
        
        #default constants
        self.FWHMconstant = 0.3
        self.polynomialOrder = 9
        self.f0 = 6668519200 
        self.calibrationScale = self.calibrationScales[self.location]
        self.FreqStart = self.scan["FreqStart"]
        
        #start window frame
        self.plotFrame = frame(self.window,(1000,1000), LEFT)
        self.infoFrame = frame(self.window,(1000,1000), None)
        self.masterFrame_1 = frame(self.window,(1000,1000), BOTTOM)
        self.masterFrame = frame(self.masterFrame_1,(1000,1000), None)
        
        self.startDataPlotButton = Button (self.infoFrame, text="Plot data points",  command=self.plotDataPoints, activebackground="Green", background="Green", font=self.font)
        self.startChangeData = Button (self.infoFrame, text="Change Data", command=self.changeData, activebackground="Blue", background="Blue", font=self.font)
        self.startDataPlotButton.pack(side=BOTTOM, fill=BOTH)
        self.startChangeData.pack(fill=BOTH)  
        
        #infoFrame
        self.window.title("Info")
        infoPanelLabelsText = ["Experiment name: " + self.expername, "Scan number: " + self.scanNumber, "Source: " + self.source, "Station: " + self.location, "Date: " + self.scan["dates"], "Start time: " + self.scan["startTime"], "Stop time: " + self.scan["stopTime"], "System temperature 1u: " + str(self.Systemtemperature1u), "System temperature 9u: " + str(self.Systemtemperature9u), "Frequency Start: " + str(self.scan["FreqStart"]), "f0", "Calibration scale", "FWHM constant", "Polynomial order"]
        infoPanelEntryText = [{"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"addEntry":False}, {"defaultValue":self.Systemtemperature1u,"addEntry":True}, {"defaultValue":self.Systemtemperature9u,"addEntry":True}, {"defaultValue":self.scan["FreqStart"],"addEntry":True}, {"defaultValue":self.f0, "addEntry":True}, {"defaultValue":str(self.calibrationScale), "addEntry":True}, {"defaultValue":str(self.FWHMconstant), "addEntry":True}, {"defaultValue":str(self.polynomialOrder), "addEntry":True}]
        
        for i in range(0, len(infoPanelLabelsText)): 
            self.infoLabel = Label(self.infoFrame, text=infoPanelLabelsText[i], anchor=W, justify=LEFT, font=self.font_2)
            self.infoLabel.pack(fill=BOTH)
            
            if  infoPanelEntryText[i]["addEntry"]:
                self.infoInputField = Entry(self.infoFrame, font=self.font_2)
                self.infoInputField.insert(0, str(infoPanelEntryText[i]["defaultValue"]))
                self.infoInputField.pack(fill=BOTH)
                
        self.plot_1 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_1.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_1.plot(self.xarray, self.y1array, 'ko', label='Data Points', markersize=1, picker=5)
        
        self.plot_2 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_2.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_2.plot(self.xarray, self.y2array, 'ko', label='Data Points', markersize=1, picker=5)
                
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
        
        messagebox.showinfo("", "Data was changed")
    
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
        if self.state == 0:
            self.plot_3.removePolt()
            self.plot_4.removePolt()
            self.plotFrame.destroy()
            self.masterFrame.destroy()
            self.createPolynomialButton.destroy()
            self.backButton.destroy()
            self.mSlider.destroy()
            self.nSlider.destroy()
            
            del self.plot_3
            del self.plot_4
            del self.plotFrame
            del self.masterFrame
            del self.createPolynomialButton
            del self.backButton
            del self.mSlider
            del self.nSlider
            
            self.startWindow()
            
        elif self.state == 1:
            self.plot_5.removePolt()
            self.plot_6.removePolt()
            self.plotLocalMaximumButton.destroy()
            self.backButton_1.destroy()
            
            del self.plot_5
            del self.plot_6
            del self.plotLocalMaximumButton
            del self.backButton_1
            
            self.plotDataPoints()
            
        elif self.state == 2:
            self.plot_7.removePolt()
            self.plot_8.removePolt()
            self.plot_9.removePolt()
            self.backButton_2.destroy()
            self.monitoringButton.destroy()
            
            del self.plot_7
            del self.plot_8
            del self.plot_9
            del self.backButton_2
            del self.monitoringButton
            
            self.plotPolynomial()
                
    def updateEnd(self, event):
        self.plot_3.plot(self.xarray[int(self.previousM)], self.z1[int(self.previousM)], 'ko', markersize=1)
        self.plot_4.plot(self.xarray[int(self.previousM)], self.z2[int(self.previousM)], 'ko', markersize=1)
        
        self.plot_3.plot(self.xarray[int(self.previousN-1)], self.z1[int(self.previousN-1)], 'ko', markersize=1)
        self.plot_4.plot(self.xarray[int(self.previousN-1)], self.z2[int(self.previousN-1)], 'ko', markersize=1)
        
        self.plot_3.annotation(self.xarray[int(self.previousM)], self.z1[int(self.previousM)], " ")
        self.plot_4.annotation(self.xarray[int(self.previousM)], self.z2[int(self.previousM)], " ")
   
        self.plot_3.annotation(self.xarray[int(self.previousN-1)], self.z1[int(self.previousN-1)], " ")
        self.plot_4.annotation(self.xarray[int(self.previousN-1)], self.z2[int(self.previousN-1)], " ")

        self.plot_3.remannotation()
        self.plot_4.remannotation()

        self.plot_3.annotation(self.xarray[int(self.mSlider.get())], self.z1[int(self.mSlider.get())], "M")
        self.plot_4.annotation(self.xarray[int(self.mSlider.get())], self.z2[int(self.mSlider.get())], "M")
        
        self.plot_3.annotation(self.xarray[int(self.nSlider.get()-1)], self.z1[int(self.nSlider.get()-1)], "N")
        self.plot_4.annotation(self.xarray[int(self.nSlider.get()-1)], self.z2[int(self.nSlider.get()-1)], "N")
        
        self.plot_3.plot(self.xarray[int(self.mSlider.get())], self.z1[int(self.mSlider.get())], 'ro', markersize=1)
        self.plot_4.plot(self.xarray[int(self.mSlider.get())], self.z2[int(self.mSlider.get())], 'ro', markersize=1)
        
        self.plot_3.plot(self.xarray[int(self.nSlider.get()-1)], self.z1[int(self.nSlider.get()-1)], 'ro', markersize=1)
        self.plot_4.plot(self.xarray[int(self.nSlider.get()-1)], self.z2[int(self.nSlider.get()-1)], 'ro', markersize=1)
        
        self.plot_3.canvasShow()
        self.plot_4.canvasShow()
 
        self.previousM = self.mSlider.get()
        self.previousN = self.nSlider.get()
           
    def plotDataPoints (self):
        self.state = 0
        
        try:
            self.plot_1.removePolt()
            self.plot_2.removePolt()
            self.startChangeData.destroy()
            self.startDataPlotButton.destroy()
            self.masterFrame.destroy()
            self.infoFrame.destroy()
            
            del self.plot_1
            del self.plot_2
            del self.startChangeData
            del self.startDataPlotButton
            del self.masterFrame
            del self.infoFrame
            
        except:
            pass
               
        self.calibration()
         
        self.masterFrame = frame(self.masterFrame_1,(1000,1000), BOTTOM)
        self.window.title("Data points for " + self.expername + " scan " +  self.scanNumber +  " for Source " + self.source)
    
        self.createPolynomialButton = Button (self.masterFrame, text="Create Polynomial", command=self.plotPolynomial, activebackground="Green", background="Green", font=self.font)
        self.createPolynomialButton.pack(fill=BOTH)
        self.backButton = Button (self.masterFrame, text="Back", command=self.back, activebackground="Blue", background="Blue", font=self.font)
        self.backButton.pack(fill=BOTH)
        
        self.points_1u = list()
        self.points_9u = list()

        #u1
        self.plot_3 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_3.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_3.plot(self.xarray, self.z1, 'ko', label='Data Points', markersize=1, picker=5)
        self.plot_3.addPickEvent(self.onpickU1)
        self.plot_3.addSecondAxis("x", "Data points", 0, self.dataPoints + 512, 1024)
        #self.plot_1.addSlider([0.10, 0.15, 0.65, 0.03], "m", 10, 0, 5, None)

        #u9
        self.plot_4 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_4.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_4.plot(self.xarray, self.z2, 'ko', label='Data Points', markersize=1, picker=5)
        self.plot_4.addPickEvent(self.onpickU9)
        self.plot_4.addSecondAxis("x", "Data points", 0, self.dataPoints + 512, 1024)
        #self.plot_2.addSlider([0.10, 0.25, 0.65, 0.03], "n", 10, 0, 5, None)
        
        #sliders
        self.previousM = self.m
        self.previousN = self.n -1
        
        self.mSlider = tk.Scale(self.plotFrame, from_= self.m, to = self.a-1, orient=HORIZONTAL, label="M", length=400, variable=self.m, command=self.updateEnd, foreground="Red", highlightcolor="Yellow")
        self.mSlider.pack(side=BOTTOM)
        self.m = self.mSlider.get()
        
        self.nSlider = tk.Scale(self.plotFrame, from_ = self.b-1 , to = self.n, orient=HORIZONTAL, label="N", length=400, variable=self.n, command=self.updateEnd, foreground="Red", highlightcolor="Yellow")
        self.nSlider.pack(side=BOTTOM)
        self.nSlider.set(self.n)
        self.n = self.nSlider.get()
        
    def onpickU1(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_3.plot(p[0][0], p[0][1], 'ro',  markersize=1,  picker=5)
        self.points_1u.append(p[0])
        #self.points_9u.append(p[0])
        self.plot_3.canvasShow()
        
    def onpickU9(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_4.plot(p[0][0], p[0][1], 'ro', markersize=1, picker=5)
        self.points_9u.append(p[0])
        #self.points_1u.append(p[0])
        self.plot_4.canvasShow()
        
    def onpick_maxU1(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.maxu1_index.append(ind[0])
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_7.plot(p[0][0], p[0][1], 'gd', markersize=2, picker=5)
        if  self.maxU1.count(p[0]) == 0:
            self.maxU1.append(p[0])
        self.plot_7.canvasShow()
        
    def onpick_maxU9(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.maxu9_index.append(ind[0])
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_8.plot(p[0][0], p[0][1], 'gd', markersize=2, picker=5)
        if  self.maxU9.count(p[0]) == 0:
            self.maxU9.append(p[0])
        self.plot_8.canvasShow()
        
    def onpick_maxAVG(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        self.maxavg_index.append(ind[0])
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_9.plot(p[0][0], p[0][1], 'gd', markersize=2, picker=5)
        if  self.avgMax.count(p[0]) == 0:
            self.avgMax.append(p[0])
        self.plot_9.canvasShow()
        
    def plotPolynomial(self):
        self.state = 1
        
        try:
            #nodzes ieprieksejos grafikus
            self.plot_3.removePolt()
            self.plot_4.removePolt()
            
            self.plot_3.removePickEvent()
            self.plot_4.removePickEvent()
            
            self.createPolynomialButton.destroy()
        
            self.m = self.mSlider.get()
            self.n = self.nSlider.get()
            
            self.mSlider.destroy()
            self.nSlider.destroy()
            
            self.backButton.destroy()
            self.masterFrame.destroy()
            
            del self.plot_3
            del self.plot_4
            del self.mSlider
            del self.nSlider
            del self.backButton
            del self.masterFrame
             
        except:
            pass
        
        self.masterFrame = frame(self.masterFrame_1,(1000,1000), BOTTOM)
        self.window.title("Polynomial " + self.expername + " scan " +  self.scanNumber +  " for Source " + self.source)
        
        self.plotLocalMaximumButton = Button (self.masterFrame, text="Create local maximum", command=self.plotLocalMaximum, activebackground="Green", background="Green", font=self.font)
        self.plotLocalMaximumButton.pack(fill=BOTH)
        
        self.backButton_1 = Button (self.masterFrame, text="Back", command=self.back, activebackground="Blue", background="Blue", font=self.font)
        self.backButton_1.pack(fill=BOTH)
        
        self.xarray_u1 = self.xarray
        self.xarray_u9 = self.xarray
        
        print ("pirms dzesanas ", self.xarray_u9.shape[0])
        
        bad_u1_indexies = list()
        bad_u9_indexies = list()
        
        for bad in range(0, len(self.points_1u)):
            bad_u1_indexies.append(indexies(self.xarray_u1, self.points_1u[bad][0]))
        for bad in range(0, len(self.points_9u)):
            bad_u9_indexies.append(indexies(self.xarray_u9, self.points_9u[bad][0]))
            
        bad_u1_indexies, bad_u9_indexies
        
        for p_u1 in self.points_1u:
            self.xarray_u1 = np.delete(self.xarray_u1, self.xarray_u1[self.xarray_u1 == p_u1[0]])
            self.z1 = np.delete(self.z1, self.z1[self.z1 == p_u1[1]])
            
        for p_u9 in self.points_9u:
            self.xarray_u9 = np.delete(self.xarray_u9, self.xarray_u9[self.xarray_u9 == p_u9[0]])
            self.z2 = np.delete(self.z2, self.z2[self.z2 == p_u9[1]])
            
        for bad in range(0, len(bad_u9_indexies)):
            self.xarray_u1 = np.delete(self.xarray_u1, self.xarray_u1[bad_u9_indexies[bad]])
            self.z1 = np.delete(self.z1, self.xarray_u9[bad_u9_indexies[bad]])
                                       
        for bad in range(0, len(bad_u1_indexies)):
            self.xarray_u9 = np.delete(self.xarray_u9, self.xarray_u9[bad_u1_indexies[bad]])
            self.z2 = np.delete(self.z2, self.xarray_u9[bad_u1_indexies[bad]])
                                       
        self.dataPoints_u1 = self.xarray_u1.shape[0]
        self.dataPoints_u9 = self.xarray_u9.shape[0]
        
        self.a_u1, self.b_u1 = FWHM(self.xarray_u1, self.z1, self.FWHMconstant)
        self.a_u9, self.b_u9 = FWHM(self.xarray_u9, self.z2, self.FWHMconstant)
         
        print ("pec dzesanas ", self.xarray_u9.shape[0])
        
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
                    print ("MJD: \t", mjd)
                elif vards == "Vobs":
                    Vobs = Header[1]
                    print ("Vobs: \t", Vobs)
                elif vards == "AtFreq":
                    AtFreq = Header[1]
                    print ("At Freq: \t", AtFreq)
                elif vards == "FreqShift":
                    FreqShift = Header[1]
                    print ("FreqShift: \t", FreqShift)
                elif vards == "VelTotal":
                    VelTotal = float(Header[1])
                    print ("VelTotal: \t", VelTotal)
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
        
        # Fit the data using poly fit
        x_u1 = np.append(self.x_u1[self.m:self.a_u1], self.x_u1[self.b_u1:self.n])
        x_u9 = np.append(self.x_u9[self.m:self.a_u9], self.x_u9[self.b_u9:self.n])
        y_u1 = np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n])
        y_u9 = np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n])
        z_u1 = np.polyfit(x_u1, y_u1, self.polynomialOrder)
        z_u9 = np.polyfit(x_u9, y_u9, self.polynomialOrder)
       
        p_u1 = np.poly1d(z_u1)
        p_u9 = np.poly1d(z_u9)
        
        #u1 plot
        self.plot_5 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_5.creatPlot(LEFT, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_5.plot(np.append(self.x_u1[self.m:self.a_u1], self.x_u1[self.b_u1:self.n]), np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]), 'ko', label='Data Points',  markersize=1)
        self.plot_5.plot(self.x_u1[self.m:self.n], self.ceb_1(self.x_u1[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        #self.plot_5.plot(self.x_u1[self.m:self.n], p_u1(self.x_u1[self.m:self.n]), 'b', label='Numpy polyfit', markersize=1)
        
        #u9 plot
        self.plot_6 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_6.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_6.plot(np.append(self.x_u9[self.m:self.a_u9], self.x_u9[self.b_u9:self.n]), np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]), 'ko', label='Data Points',  markersize=1)
        self.plot_6.plot(self.x_u9[self.m:self.n], self.ceb_2(self.x_u9[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        #self.plot_6.plot(self.x_u9[self.m:self.n], p_u9(self.x_u9[self.m:self.n]), 'b', label='Numpy polyfit', markersize=1)
        
    def plotLocalMaximum(self):
        self.state = 2
        #nodzes ieprieksejos grafikus
        try:
            self.plot_5.removePolt()
            self.plot_6.removePolt()
            self.backButton_1.destroy()
            self.plotLocalMaximumButton.destroy()
            self.masterFrame.destroy()
            
            del self.plot_5
            del self.plot_6
            del self.backButton_1
            del self.plotLocalMaximumButton
            del self.masterFrame
            
        except:
            pass
        
        self.masterFrame = frame(self.masterFrame_1,(1000,1000), BOTTOM)
        self.window.title("Local maximums " + self.expername + " scan " +  self.scanNumber +  " for Source " + self.source)
        
        self.monitoringButton = Button (self.masterFrame, text="Add points to monitoring", command=self.createResult, activebackground="Green", background="Green", font=self.font)
        self.monitoringButton.pack(fill=BOTH)
        
        self.backButton_2 = Button (self.masterFrame, text="Back", command=self.back, activebackground="Blue", background="Blue", font=self.font)
        self.backButton_2.pack(fill=BOTH)
        
        thres=0.1
    
        y1values = self.z1[self.m:self.n] - self.ceb_1(self.x_u1[self.m:self.n])
        y2values = self.z2[self.m:self.n] - self.ceb_2(self.x_u9[self.m:self.n])
        
        #indexsu apreikinasana
        indexes_for_ceb = peakutils.indexes(y1values, thres=thres, min_dist=10)
        indexes_for_ceb2 = peakutils.indexes(y2values, thres=thres, min_dist=10)
        
        #u1
        self.plot_7 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_7.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_7.plot(self.x_u1[self.m:self.n], y1values, 'b', label='Signal - polynomial', markersize=1)
        self.plot_7.plot(self.x_u1[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb], 'dr', label="Local Maximums for signal", markersize=2, picker=5)
        self.plot_7.addPickEvent(self.onpick_maxU1)
        self.plot_7.annotations(self.x_u1[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb])
        
        #u9
        self.plot_8 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_8.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_8.plot(self.x_u9[self.m:self.n], y1values, 'b', label='Signal - polynomial', markersize=1)
        self.plot_8.plot(self.x_u9[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2], 'dr', label="Local Maximums for signal", markersize=2, picker=5)
        self.plot_8.addPickEvent(self.onpick_maxU9)
        self.plot_8.annotations(self.x_u9[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2])
        
        #mid plot
        avg_x = (self.x_u1[self.m:self.n] + self.x_u9[self.m:self.n]) / 2
        avg_y = (y1values + y2values) / 2
        indexes_for_avg = peakutils.indexes(avg_y, thres=thres, min_dist=10)
        
        self.plot_9 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_9.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Average Polarization")
        self.plot_9.plot(avg_x, avg_y, 'b', label='Signal - polynomial', markersize=1)
        self.plot_9.plot(avg_x[indexes_for_avg], avg_y[indexes_for_avg], 'dr', label="Local Maximums for signal", markersize=2, picker=5)
        self.plot_9.addPickEvent(self.onpick_maxAVG)
        self.plot_9.annotations(avg_x[indexes_for_avg], avg_y[indexes_for_avg])
        
        self.maxU1 = list()
        self.maxU9 = list()
        self.avgMax = list()
        self.maxu1_index = list()
        self.maxu9_index = list()
        self.maxavg_index = list()
        
    def createResult(self):
        
        try:
            #remove graph
            self.plot_7.removePolt()
            self.plot_8.removePolt()
            self.plot_9.removePolt()
            
            self.monitoringButton.destroy()
            self.masterFrame.destroy()
            
            del self.plot_7
            del self.plot_8
            del self.plot_9
            del self.monitoringButton
            del self.masterFrame
            
        except:
            pass
        
        self.masterFrame = frame(self.window,(1000,1000), BOTTOM)
        endLabel = Label(master=self.masterFrame, text="Result file creating in progress!")
        endLabel.pack()
        
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
        
        self.maxU1.sort(key=lambda tup: tup[0], reverse=True)
        self.maxU9.sort(key=lambda tup: tup[0], reverse=True)
        self.avgMax.sort(key=lambda tup: tup[0], reverse=True)
              
        result[self.expername][self.scanNumber]["startTime"] = self.scan["startTime"]
        result[self.expername][self.scanNumber]["stopTime"] = self.scan["stopTime"]
        result[self.expername][self.scanNumber]["location"] = self.location
        result[self.expername][self.scanNumber]["Date"] = self.scan["dates"]
                
        result[self.expername][self.scanNumber]["polarizationU1"] = self.maxU1
        result[self.expername][self.scanNumber]["polarizationU9"] = self.maxU9
        result[self.expername][self.scanNumber]["polarizationAVG"] = self.avgMax
                
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
    try:
        data = np.fromfile(dataFileName, dtype="float64", count=-1, sep=" ") .reshape((file_len(dataFileName),5))
        
    except IOError as e:
        print ("IO Error",  e)
        sys.exit(1)
            
    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(1)
             
    else:
        data = np.delete(data, (0), axis=0) #izdzes masiva primo elementu
        dataPoints = data.shape[0]
        
        xdata = data[:, [0]]
        y1data = data[:, [1]] 
        y2data = data[:, [2]]
    
        return (xdata, y1data, y2data, dataPoints)

def getLogs(logfileName, dataFileName, singleSourceExperiment, prettyLogsPath): 
    logs  = ExperimentLogReader(logfileName, prettyLogsPath, singleSourceExperiment).getLogs()
    scanNumber = dataFileName.split(".")[0].split("_")[-1][1:len(dataFileName)]
    
    try:
        scan = logs[scanNumber]
        
    except KeyError as e:
        print ("KeyError",  e)
        sys.exit(1)
            
    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(1)
        
    else:
    
        Systemtemperature1u = float(scan["Systemtemperature"][0])
        Systemtemperature9u = float(scan["Systemtemperature"][1])
        
        location = logs["location"]
        source = scan["source"]
       
        return (Systemtemperature1u, Systemtemperature9u, location, source, scan, scanNumber)
    
def main():
    # Parse the arguments
    args = parseArguments()
    configFilePath = str(args.__dict__["config"])
    logFileName = str(args.__dict__["logFile"])
    dataFileName = str(args.__dict__["datafile"])
    singleSourceExperiment = list(args.__dict__["single"])
    
    #Creating config parametrs
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    logPath = config.get('paths', "logPath")
    prettyLogsPath =  config.get('paths', "prettyLogsPath")
    dataFilesPath =  config.get('paths', "dataFilePath")
    
    calibrationScales = {"IRBENE":config.getint('parametrs', "irbene"), "IRBENE16":config.getint('parametrs', "irbene16")}

    #get Data and Logs
    try:
        Systemtemperature1u, Systemtemperature9u, location, source, scan, scanNumber = getLogs(logPath + logFileName, dataFilesPath + dataFileName, singleSourceExperiment, prettyLogsPath)
        expername = logFileName.split(".")[0][:-2]
        xdata, y1data, y2data, dataPoints = getData(dataFilesPath + dataFileName)
    
    except TypeError as e:
        print ("TypeError error:", e)
        sys.exit(1)
         
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise
        sys.exit(1)
        
    else:
        #Create App
        window = tk.Tk()
        window.configure(background='light goldenrod')
        ploting = MaserPlot(window, xdata, y1data, y2data, dataPoints, Systemtemperature1u, Systemtemperature9u, expername, source, location, scan, scanNumber, calibrationScales)
        img = tk.Image("photo", file="viraclogo.png")
        window.call('wm','iconphoto', window._w,img)
        
        ploting.mainloop()
        
        sys.exit(0)

if __name__=="__main__":
    main()
