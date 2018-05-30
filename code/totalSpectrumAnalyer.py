#! /usr/bin/python
import sys
import os
import argparse
import configparser
from tkinter import *
import tkinter as tk
from tkinter import font, messagebox
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.modeling import fitting
from astropy.modeling.polynomial import Chebyshev1D
from scipy.interpolate import UnivariateSpline
import peakutils
import json
import re

from ploting import Plot

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''plotting tool. ''', epilog="""PRE PLOTTER.""")
    
    # Positional mandatory arguments
    parser.add_argument("datafile", help="Experiment correlation file name", type=str)

    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")

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

def indexies(array, value):
    indexs = list()
    for i in range(0, len(array)-1):
        if array[i] == value:
            indexs.append(i)
    return indexs

def frame(parent, size, sides, **options):
    Width=size[0]
    Height=size[1]
    f=tk.Frame(parent, width=Width, height=Height, background="light goldenrod", **options)
    f.pack(side = sides, fill=BOTH)
    return (f)

def FWHM(x, y, constant):
    spline = UnivariateSpline(x, y-np.max(y)/2, k=3, s=20)
    spline.set_smoothing_factor(0.5)
    root1 = spline.roots()[0] - constant
    root2 = spline.roots()[-1] + constant
    index_1 =  (np.abs(x-root1)).argmin()
    index_2 =  (np.abs(x-root2)).argmin()
    return (index_1, index_2)
    
    '''
    max = np.max(y)
    std = np.std(y)
    root1 = max + 2*std  #+ constant
    root2 = max - 2*std  #- constant
    
    index_1 =  (np.abs(y-root1)).argmin() - constant
    index_2 =  (np.abs(y-root2)).argmin() + constant
    print ("roots ", root1, root2)
    print ("Indexies ", index_1, index_2)
    #index_1 = 900
    #index_2 = 1500
    return (index_1, index_2)
    '''
   
class Analyzer(Frame):
    def __init__(self, window, datafile, resultFilePath):
        Frame.__init__(self)
        self.font_2 = font.Font(family="Times New Roman", size=10)
        self.window = window
        self.FWHMconstant = 1
        self.polynomialOrder = 3
        self.source = re.split("([A-Z, a-z]+)", datafile.split("/")[-1].split(".")[0])[1]
        self.expername = datafile.split("/")[-1].split(".")[0]
        self.date = re.split("([A-Z, a-z]+)", datafile.split("/")[-1].split(".")[0])[2][0:-1]
        self.location = datafile.split("/")[-1].split(".")[0].split("_")[-1]
        self.resultFilePath = resultFilePath
        
        self.infoSet = set()
        
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
            self.m = 0
            self.n = self.dataPoints
            
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
                
            self.xarray =  np.flip(self.xarray,0)
            self.y1array =  np.flip(self.y1array,0)
            self.y2array =  np.flip(self.y2array,0)
            
            self.font = font.Font(family="Times New Roman", size=20, weight=font.BOLD)    
            self.masterFrame = frame(self.window,(1000,1000), RIGHT)
            
            self.plotInitData()
            
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
        
    def changeData(self):
        childs = self.masterFrame.winfo_children()
        newValues = list()
        
        for child in childs:
            if child.winfo_class() == "Entry":
                newValues.append(child.get())
                
        self.FWHMconstant = int(newValues[0])
        self.polynomialOrder = int(newValues[1])
        
        messagebox.showinfo("", "Data was changed")
                
    def plotInitData(self):
        self.window.title("Info")
        
        #self.infoFrame = frame(self.window,(1000,1000), LEFT)
        self.plotFrame = frame(self.window,(1000,1000), None)
        self.startChangeData = Button (self.masterFrame, text="Change Data", command=self.changeData, activebackground="Blue", background="Blue", font=self.font)
        self.startChangeData.pack(fill=BOTH)  
        
        self.plot_1 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_1.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_1.plot(self.xarray, self.y1array, 'ko', label='Data Points', markersize=1)
        
        self.plot_2 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_2.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_2.plot(self.xarray, self.y2array, 'ko', label='Data Points', markersize=1)
        
        self.plotSmoothData = Button (self.masterFrame, text="Smooth Data", command=self.plotSmoothData, activebackground="Green", background="Green", font=self.font)
        self.plotSmoothData.pack(fill=BOTH)
        
        infoPanelLabelsText = ["FWHM constant", "Polynomial order"]
        infoPanelEntryText = [ {"defaultValue":str(self.FWHMconstant), "addEntry":True}, {"defaultValue":str(self.polynomialOrder), "addEntry":True}]
        
        for i in range(0, len(infoPanelLabelsText)):
            
            self.infoLabel = Label(self.masterFrame, text=infoPanelLabelsText[i], anchor=W, justify=LEFT, font=self.font_2)
            self.infoLabel.pack(fill=BOTH)
            self.infoSet.add(self.infoLabel)
            
            if  infoPanelEntryText[i]["addEntry"]:
                self.infoInputField = Entry(self.masterFrame, font=self.font_2)
                self.infoInputField.insert(0, str(infoPanelEntryText[i]["defaultValue"]))
                self.infoInputField.pack(fill=BOTH)
                self.infoSet.add(self.infoInputField)
          
    def plotSmoothData(self):
        self.window.title("Smooth Data")
        
        self.infoLabel.destroy()
        self.startChangeData.destroy()
        del self.startChangeData
        self.plotSmoothData.destroy()
        del self.plotSmoothData
        
        while len(self.infoSet) !=0:
            info_item = self.infoSet.pop()
            info_item.destroy()
            del info_item
            
        del self.infoSet
        
        self.plotPolinomial = Button (self.masterFrame, text="Plot polinomial", command=self.plotPlonomials, activebackground="Green", background="Green", font=self.font)
        self.plotPolinomial.pack(fill=BOTH)
        
        g1 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
        g2 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
    
        self.z1 = convolve(self.y1array, g1, boundary='extend')
        self.z2 = convolve(self.y2array, g2, boundary='extend')
        
        self.plot_1.removePolt()
        self.plot_2.removePolt()
        del self.plot_1
        del self.plot_2
        
        self.points_1u = list()
        self.points_9u = list()
        
        #self.plotFrame = frame(self.window,(1000,1000), TOP)
        self.plot_3 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_3.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_3.plot(self.xarray, self.z1, 'ko', label='Data Points', markersize=1, picker=5)
        self.plot_3.addPickEvent(self.onpickU1)
        #self.plot_3.addSecondAxis("x", "Data points", 0, self.dataPoints + 512, 1024)
        
        self.plot_4 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_4.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_4.plot(self.xarray, self.z2, 'ko', label='Data Points', markersize=1, picker=5)
        self.plot_4.addPickEvent(self.onpickU9)
        #self.plot_4.addSecondAxis("x", "Data points", 0, self.dataPoints + 512, 1024)
        
        self.a, self.b = FWHM(self.xarray, (self.z1 + self.z2)/2, self.FWHMconstant)
        
        #sliders
        self.previousM = self.m
        self.previousN = self.n -1
        
        self.mSlider = tk.Scale(self.masterFrame, from_= self.m, to = self.a-1, orient=HORIZONTAL, label="M", length=400, variable=self.m, command=self.updateEnd, foreground="Red", highlightcolor="Yellow")
        self.mSlider.pack(side=BOTTOM, fill=BOTH)
        self.m = self.mSlider.get()
        
        self.nSlider = tk.Scale(self.masterFrame, from_ = self.b-1 , to = self.n, orient=HORIZONTAL, label="N", length=400, variable=self.n, command=self.updateEnd, foreground="Red", highlightcolor="Yellow")
        self.nSlider.pack(side=BOTTOM, fill=BOTH)
        self.nSlider.set(self.n)
        self.n = self.nSlider.get()
        
    def plotPlonomials(self):
        self.window.title("Polynomial and Data points")
        self.m = self.mSlider.get()
        self.n = self.nSlider.get()
        self.mSlider.destroy()
        self.nSlider.destroy()
        self.plot_3.removePickEvent()
        self.plot_4.removePickEvent()
        self.plot_3.removePolt()
        self.plot_4.removePolt()
        del self.plot_3
        del self.plot_4
        
        self.plotPolinomial.destroy()
        del self.plotPolinomial
        
        self.plotLocalMaximum = Button (self.masterFrame, text="Plot local maximums", command=self.plotLocalMaximum, activebackground="Green", background="Green", font=self.font)
        self.plotLocalMaximum.pack(fill=BOTH)
     
        print ("Before deliting  ", self.xarray.shape[0])
        bad_u1_indexies = list()
        bad_u9_indexies = list()
        
        for bad in range(0, len(self.points_1u)):
            bad_u1_indexies.append(indexies(self.xarray, self.points_1u[bad][0]))
        for bad in range(0, len(self.points_9u)):
            bad_u9_indexies.append(indexies(self.xarray, self.points_9u[bad][0]))
            
        bad_u1_indexies.sort()
        bad_u9_indexies.sort()
        
        bad_list_indexies = bad_u1_indexies + bad_u9_indexies
        bad_list_indexies = np.unique(bad_list_indexies, return_index=False, return_inverse=False, return_counts=False,)
        
        for bad in range(0, len(bad_list_indexies)):
            self.xarray = np.delete(self.xarray, self.xarray[bad_list_indexies[bad]])
            self.z1 = np.delete(self.z1, self.z1[bad_list_indexies[bad]])
            self.z2 = np.delete(self.z2, self.z2[bad_list_indexies[bad]])
            
        print ("After deliting  ", self.xarray.shape[0])
                                       
        self.a_u1, self.b_u1 = FWHM(self.xarray, self.z1, self.FWHMconstant)
        self.a_u9, self.b_u9 = FWHM(self.xarray, self.z2, self.FWHMconstant)
         
        # Fit the data using a Chebyshev astro py
        ceb = Chebyshev1D(self.polynomialOrder, domain=None, window=[-1, 1], n_models=None, model_set_axis=None, name=None, meta=None)
        fit_ceb = fitting.LevMarLSQFitter()
        
        ### u1
        self.ceb_1 = fit_ceb(ceb, np.append(self.xarray[self.m:self.a_u1], self.xarray[self.b_u1:self.n]),  np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]))
       
        ### u9
        self.ceb_2 = fit_ceb(ceb, np.append(self.xarray[self.m:self.a_u9], self.xarray[self.b_u9:self.n]),  np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]))
        
        #u1 plot
        self.plot_5 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_5.creatPlot(LEFT, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_5.plot(np.append(self.xarray[self.m:self.a_u1], self.xarray[self.b_u1:self.n]), np.append(self.z1[self.m:self.a_u1], self.z1[self.b_u1:self.n]), 'ko', label='Data Points',  markersize=1)
        self.plot_5.plot(self.xarray[self.m:self.n], self.ceb_1(self.xarray[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        #self.plot_5.plot(self.x_u1[self.m:self.n], p_u1(self.x_u1[self.m:self.n]), 'b', label='Numpy polyfit', markersize=1)
        
        #u9 plot
        self.plot_6 = Plot(6,6, self.masterFrame, self.plotFrame)
        self.plot_6.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_6.plot(np.append(self.xarray[self.m:self.a_u9], self.xarray[self.b_u9:self.n]), np.append(self.z2[self.m:self.a_u9], self.z2[self.b_u9:self.n]), 'ko', label='Data Points',  markersize=1)
        self.plot_6.plot(self.xarray[self.m:self.n], self.ceb_2(self.xarray[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        #self.plot_6.plot(self.x_u9[self.m:self.n], p_u9(self.x_u9[self.m:self.n]), 'b', label='Numpy polyfit', markersize=1)
        
    def plotLocalMaximum(self):
        self.plotLocalMaximum.destroy()
        del self.plotLocalMaximum
        self.plot_5.removePolt()
        self.plot_6.removePolt()
        del self.plot_5
        del self.plot_6
        self.window.title("Local maximums")
        self.monitoringButton = Button (self.masterFrame, text="Add points to monitoring", command=self.createResult, activebackground="Green", background="Green", font=self.font)
        self.monitoringButton.pack(fill=BOTH)
        
        thres=0.1
        
        #self.z1 = self.z1.reshape(len(self.z1), 1)
        #self.z2 = self.z2.reshape(len(self.z2), 1)
        
        y1values = self.z1[self.m:self.n] - self.ceb_1(self.xarray[self.m:self.n])
        y2values = self.z2[self.m:self.n] - self.ceb_2(self.xarray[self.m:self.n])
          
        #indexsu apreikinasana
        indexes_for_ceb = peakutils.indexes(y1values, thres=thres, min_dist=10)
        indexes_for_ceb2 = peakutils.indexes(y2values, thres=thres, min_dist=10)
        
        self.plot_7 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_7.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization")
        self.plot_7.plot(self.xarray[self.m:self.n], y1values, 'b', label='Signal - polynomial', markersize=1)
        self.plot_7.plot(self.xarray[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb], 'dr', label="Local Maximums for signal", markersize=2, picker=5)
        self.plot_7.addPickEvent(self.onpick_maxU1)
        self.plot_7.annotations(self.xarray[self.m:self.n][indexes_for_ceb], y1values[indexes_for_ceb])
        
        #u9
        self.plot_8 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_8.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization")
        self.plot_8.plot(self.xarray[self.m:self.n], y1values, 'b', label='Signal - polynomial', markersize=1)
        self.plot_8.plot(self.xarray[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2], 'dr', label="Local Maximums for signal", markersize=2, picker=5)
        self.plot_8.addPickEvent(self.onpick_maxU9)
        self.plot_8.annotations(self.xarray[self.m:self.n][indexes_for_ceb2], y2values[indexes_for_ceb2])
        
        #mid plot
        avg_y = (y1values + y2values) / 2
        
        indexes_for_avg = peakutils.indexes(avg_y, thres=thres, min_dist=10)
        
        self.plot_9 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_9.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Average Polarization")
        self.plot_9.plot(self.xarray[self.m:self.n], avg_y, 'b', label='Signal - polynomial', markersize=1)
        self.plot_9.plot(self.xarray[self.m:self.n][indexes_for_avg], avg_y[indexes_for_avg], 'dr', label="Local Maximums for signal", markersize=2, picker=5)
        self.plot_9.addPickEvent(self.onpick_maxAVG)
        self.plot_9.annotations(self.xarray[self.m:self.n][indexes_for_avg], avg_y[indexes_for_avg])
        
        self.maxU1 = list()
        self.maxU9 = list()
        self.avgMax = list()
        self.maxu1_index = list()
        self.maxu9_index = list()
        self.maxavg_index = list()
        
    def createResult(self):
        
        endLabel = Label(master=self.masterFrame, text="Result file creating in progress!")
        endLabel.pack()
        
        resultFileName = self.source + ".json"
    
        if os.path.isfile(self.resultFilePath + resultFileName):
            pass
        else:
            os.system("touch " + self.resultFilePath +  resultFileName)
            
            resultFile = open (self.resultFilePath +  resultFileName, "w")
            resultFile.write("{ \n" + "\n}")
            resultFile.close()
        
        with open(self.resultFilePath + resultFileName) as result_data:    
            result = json.load(result_data)
        
        if self.expername not in result:
            result[self.expername] = dict()
                
        self.maxU1.sort(key=lambda tup: tup[0], reverse=True)
        self.maxU9.sort(key=lambda tup: tup[0], reverse=True)
        self.avgMax.sort(key=lambda tup: tup[0], reverse=True)
              
        result[self.expername]["location"] = self.location
        result[self.expername]["Date"] = self.date
                    
        result[self.expername]["polarizationU1"] = self.maxU1
        result[self.expername]["polarizationU9"] = self.maxU9
        result[self.expername]["polarizationAVG"] = self.avgMax
                
        #result[self.expername][self.scanNumber]["index_for_polarizationU1"] =  self.maxu1_index
        #result[self.expername][self.scanNumber]["index_for_polarizationU9"] =  self.maxu9_index
        #result[self.expername][self.scanNumber]["index_for_polarizationAVG"] =  self.maxavg_index
        
        resultFile = open (self.resultFilePath +  resultFileName, "w")
        resultFile.write(json.dumps(result, indent=2))
        resultFile.close() 
        
        self. _quit()
        
    def _quit(self):
        self.masterFrame.destroy()
        self.plotFrame.destroy()
        self.window.destroy()
        
def main(): 
    
    args = parseArguments()
   
    datafile = str(args.__dict__["datafile"])
    configFilePath = str(args.__dict__["config"])
    
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    dataFilesPath =  config.get('paths', "dataFilePath")
    resultFilePath =  config.get('paths', "resultFilePath")
    
    #Create App
    window = tk.Tk()
    window.configure(background='light goldenrod')
    ploting = Analyzer(window, dataFilesPath + datafile, resultFilePath)
    img = tk.Image("photo", file="viraclogo.png")
    window.call('wm','iconphoto', window._w,img)
        
    ploting.mainloop()
        
    sys.exit(0)

if __name__=="__main__":
    main()
    