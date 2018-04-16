#! /usr/bin/python
import sys
import os
import argparse
import configparser
from tkinter import *
import tkinter as tk
from tkinter import font
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
import astropy.modeling as model
from astropy.modeling import fitting
from astropy.modeling.polynomial import Chebyshev1D
from scipy.interpolate import UnivariateSpline
import peakutils
import json

from ploting import Plot

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''plotting tool. ''', epilog="""PRE PLOTTER.""")
    
    # Positional mandatory arguments
    #parser.add_argument("source", help="Experiment source", type=str)
    #parser.add_argument("date", help="Experiment date", type=str)
    #parser.add_argument("logFile", help="Experiment log file name", type=str)
    parser.add_argument("datafile", help="Experiment correlation file name", type=str)

    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration Yaml file", type=str, default="config/config.cfg")
    parser.add_argument("-i", "--interval", help="Set interval", type=float, default=0.9)
    parser.add_argument("-t", "--threshold", help="Set threshold for outlier filter", type=float, default=1.0)
    parser.add_argument("-f", "--filter", help="Set filter default is True if filter is False bad data points is no removed", type=str, default="True")
    parser.add_argument("-s", "--single", help="Set RA, DEC, Epoch, Source name", nargs="*", type=str, default=[])
    # option -s example cepa 225617.90 620149.7 2000.0

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()
    return args

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def frame(parent, size, sides, **options):
    Width=size[0]
    Height=size[1]
    f=tk.Frame(parent, width=Width, height=Height, background="light goldenrod", **options)
    f.pack(side = sides)
    return (f)

def FWHM(x, y, constant):
    max = np.max(y)
    std = np.std(y)
    root1 = max + 2*std  #+ constant
    root2 = max - 2*std  #- constant
    
    index_1 =  (np.abs(y-root1)).argmin() #- constant
    index_2 =  (np.abs(y-root2)).argmin() #- constant
    print ("roots ", root1, root2)
    print ("Indexies ", index_1, index_2)
    #index_1 = 900
    #index_2 = 1500
    return (index_1, index_2)
   
class Analyzer(Frame):
    def __init__(self, window, datafile):
        Frame.__init__(self)
        self.font_2 = font.Font(family="Times New Roman", size=10)
        self.window = window
        self.FWHMconstant = 200
        self.polynomialOrder = 3
        
        try:
            data = np.fromfile(datafile, dtype="float64", count=-1, sep=" ") .reshape((file_len(datafile),3))
        
        except IOError as e:
            print ("IO Error",  e)
            sys.exit(1)
                
        except:
            print("Unexpected error:", sys.exc_info()[0])
            sys.exit(1)
                 
        else:
            self.dataPoints = data.shape[0]
            
            self.xdata = data[:, [0]]
            self.y_u1 = data[:, [1]] 
            self.y_u9 = data[:, [2]]
            
            #Making sure that data is numpy array
            self.xarray = np.zeros(self.dataPoints)
            self.y1array = np.zeros(self.dataPoints)
            self.y2array = np.zeros(self.dataPoints)
            
            for i in range (0, self.dataPoints):
                self.xarray[i] = self.xdata[i]
                self.y1array[i] = self.y_u1[i]
                self.y2array[i] = self.y_u9[i]
        
        self.font = font.Font(family="Times New Roman", size=20, weight=font.BOLD)    
        self.masterFrame = frame(self.window,(1000,1000), RIGHT)
        
        self.plotInitData()
            
    def plotInitData(self):
        
        #self.infoFrame = frame(self.window,(1000,1000), LEFT)
        self.plotFrame = frame(self.window,(1000,1000), None)
        self.plot_1 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_1.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_1.plot(self.xarray, self.y1array, 'ko', label='Data Points', markersize=1)
        
        self.plot_2 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_2.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_2.plot(self.xarray, self.y2array, 'ko', label='Data Points', markersize=1)
        
        self.plotSmoothData = Button (self.masterFrame, text="Smooth Data", command=self.plotSmoothData, activebackground="Blue", background="Blue", font=self.font)
        self.plotSmoothData.pack()
        
        infoPanelLabelsText = ["FWHM constant", "Polynomial order"]
        infoPanelEntryText = [ {"defaultValue":str(self.FWHMconstant), "addEntry":True}, {"defaultValue":str(self.polynomialOrder), "addEntry":True}]
        
        for i in range(0, len(infoPanelLabelsText)): 
            self.infoLabel = Label(self.masterFrame, text=infoPanelLabelsText[i], anchor=W, justify=LEFT, font=self.font_2)
            self.infoLabel.pack(fill=BOTH)
            
            if  infoPanelEntryText[i]["addEntry"]:
                self.infoInputField = Entry(self.masterFrame, font=self.font_2)
                self.infoInputField.insert(0, str(infoPanelEntryText[i]["defaultValue"]))
                self.infoInputField.pack(fill=BOTH)
          
    def plotSmoothData(self):
        self.plotSmoothData.destroy()
        del self.plotSmoothData
        
        self.plotPolinomial = Button (self.masterFrame, text="Plot polinomial", command=self.plotPlonomials, activebackground="Blue", background="Blue", font=self.font)
        self.plotPolinomial.pack()
        
        g1 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
        g2 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
    
        self.z1 = convolve(self.y1array, g1, boundary='extend')
        self.z2 = convolve(self.y2array, g2, boundary='extend')
        
        self.plot_1.removePolt()
        self.plot_2.removePolt()
        del self.plot_1
        del self.plot_2
        
        #self.plotFrame = frame(self.window,(1000,1000), TOP)
        self.plot_3 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_3.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_3.plot(self.xdata, self.z1, 'ko', label='Data Points', markersize=1, picker=5)
        
        self.plot_4 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_4.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_4.plot(self.xdata, self.z2, 'ko', label='Data Points', markersize=1, picker=5)
        
    def plotPlonomials(self):
        self.plot_3.removePolt()
        self.plot_4.removePolt()
        del self.plot_3
        del self.plot_4
        
        self.m = 0
        self.n = self.dataPoints
        
        self.a_u1, self.b_u1 = FWHM(self.xdata, self.z1, self.FWHMconstant)
        self.a_u9, self.b_u9 = FWHM(self.xdata, self.z2, self.FWHMconstant)
        
        # Fit the data using a Chebyshev astro py
        ceb = Chebyshev1D(self.polynomialOrder, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
        fit_ceb = fitting.LevMarLSQFitter()
        
        ### u1
        self.ceb_1 = fit_ceb(ceb, np.append(self.xdata[self.m:self.a_u1], self.xdata[self.b_u1:self.n]),  np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]))
       
        ### u9
        self.ceb_2 = fit_ceb(ceb, np.append(self.xdata[self.m:self.a_u9], self.xdata[self.b_u9:self.n]),  np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]))
        
        #u1 plot
        self.plot_5 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_5.creatPlot(LEFT, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_5.plot(np.append(self.xdata[self.m:self.a_u1], self.xdata[self.b_u1:self.n]), np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]), 'ko', label='Data Points',  markersize=1)
        self.plot_5.plot(self.xdata[self.m:self.n], self.ceb_1(self.xdata[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        #self.plot_5.plot(self.x_u1[self.m:self.n], p_u1(self.x_u1[self.m:self.n]), 'b', label='Numpy polyfit', markersize=1)
        
        #u9 plot
        self.plot_6 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_6.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_6.plot(np.append(self.xdata[self.m:self.a_u9], self.xdata[self.b_u9:self.n]), np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]), 'ko', label='Data Points',  markersize=1)
        self.plot_6.plot(self.xdata[self.m:self.n], self.ceb_2(self.xdata[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        #self.plot_6.plot(self.x_u9[self.m:self.n], p_u9(self.x_u9[self.m:self.n]), 'b', label='Numpy polyfit', markersize=1)
        
def main(): 
    
    args = parseArguments()
    '''
    source = str(args.__dict__["source"])
    date = str(args.__dict__["date"])
    logFile = str(args.__dict__["logFile"])
    interval = float(args.__dict__["interval"])
    threshold = float(args.__dict__["threshold"])
    filter = str(args.__dict__["filter"])
    singleSourceExperiment = list(args.__dict__["single"])
    '''
    datafile = str(args.__dict__["datafile"])
    configFilePath = str(args.__dict__["config"])
    
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    dataFilesPath =  config.get('paths', "dataFilePath")
    prettyLogsPath =  config.get('paths', "prettyLogsPath")
    logPath = config.get('paths', "logPath")
    #logs  = ExperimentLogReader(logPath + logFile, prettyLogsPath, singleSourceExperiment).getLogs()
    
    #Create App
    window = tk.Tk()
    window.configure(background='light goldenrod')
    ploting = Analyzer(window, dataFilesPath + datafile)
    img = tk.Image("photo", file="viraclogo.png")
    window.call('wm','iconphoto', window._w,img)
        
    ploting.mainloop()
        
    sys.exit(0)

if __name__=="__main__":
    main()