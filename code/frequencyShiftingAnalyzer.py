#! /usr/bin/python
import sys
import os
import argparse
import configparser
from tkinter import *
import tkinter as tk
from tkinter import font
import numpy as np
import pandas as pd
from pandas.stats.moments import rolling_mean

from ploting import Plot

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Creates input file for plotting tool. ''', epilog="""PRE PLOTTER.""")
    
    # Positional mandatory arguments
    parser.add_argument("source", help="Experiment source", type=str)
    parser.add_argument("date", help="Experiment date", type=str)

    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration Yaml file", type=str, default="config/config.cfg")
    parser.add_argument("-i", "--interval", help="Set interval", type=float, default=0.9)
    parser.add_argument("-t", "--threshold", help="Set threshold for outlier filter", type=float, default=1.0)
    parser.add_argument("-f", "--filter", help="Set filter default is True if filter is False bad data points is no removed", type=str, default="True")
    parser.add_argument("-p", "--paircount", help="Set pair count", type=int, default=2)

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

def indexies(array, value):
    indexs = list()
    for i in range(0, len(array)-1):
        if array[i] == value:
            indexs.append(i)
    return indexs

def is_outlier(points, threshold):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score < threshold

class Analyzer(Frame):
    def __init__(self, window, source, date, filter, threshold, badPointRange, dataPath):
        Frame.__init__(self)
        self.window = window
        self.source = source
        self.threshold = threshold
        self.filter = filter
        self.badPointRange = badPointRange
        self.date = date
        self.dataFileDir = dataPath + self.source + "/" + self.date
        self.scanPairs = self.createScanPairs()
        self.datPairsCount = len(self.scanPairs)
        self.font = font.Font(family="Times New Roman", size=20, weight=font.BOLD)
        self.font_2 = font.Font(family="Times New Roman", size=10)
        self.masterFrame = frame(self.window,(1000,1000), None)
        self.index = 0
        
        self.__UI__()
        
    def createScanPairs(self):
        print (self.dataFileDir)
        
        dataFiles = list()
        for dataFile in os.listdir(self.dataFileDir):
            dataFiles.append(dataFile)
        
        dataFiles.sort()
    
        scanPairs = list()
        i = 0
        j = 1
        
        for k in range(0, int(len(dataFiles) - len(dataFiles) /2)):
            scanPairs.append((dataFiles[i], dataFiles[j])) 
            i = i + 2
            j = j + 2
        
        return scanPairs
    
    def __getDataForPolarization__(self, data1, data2, filter):
        if filter == True:
            outliersMask_1 = is_outlier(data1, self.threshold)
            outliersMask_2 = is_outlier(data2, self.threshold)
            
            bad_point_index_1 = indexies(outliersMask_1, False)
            bad_point_index_2 = indexies(outliersMask_2, False)
             
            xdata_1_f = data1[:, [0]].tolist()
            xdata_2_f = data2[:, [0]].tolist()
            ydata_1_u1 = data1[:, [1]].tolist()
            ydata_2_u1 = data2[:, [1]].tolist()
            ydata_1_u9 = data1[:, [2]].tolist()
            ydata_2_u9 = data2[:, [2]].tolist()
            
            df_y1_u1 = pd.DataFrame(data=ydata_1_u1)
            df_y1_u9 = pd.DataFrame(data=ydata_1_u9)
            df_y2_u1 = pd.DataFrame(data=ydata_2_u1)
            df_y2_u9 = pd.DataFrame(data=ydata_2_u9)
            
            mean_y1_u1 = np.nan_to_num(df_y1_u1.rolling(window=self.badPointRange, center=True).mean())
            mean_y1_u9 = np.nan_to_num(df_y1_u9.rolling(window=self.badPointRange, center=True).mean())
            mean_y2_u1 = np.nan_to_num(df_y2_u1.rolling(window=self.badPointRange, center=True).mean())
            mean_y2_u9 = np.nan_to_num(df_y2_u9.rolling(window=self.badPointRange, center=True).mean())
            
            mean_y1_u1_2 = np.mean(ydata_1_u1)
            mean_y1_u9_2 = np.mean(ydata_1_u9)
            mean_y2_u1_2 = np.mean(ydata_2_u1)
            mean_y2_u9_2 = np.mean(ydata_2_u9)
                
            for badPoint in bad_point_index_1:
                ydata_1_u1[badPoint][0] = mean_y1_u1[badPoint]
                
            for badPoint in bad_point_index_1:
                ydata_1_u9[badPoint][0] = mean_y1_u9[badPoint]
                
            for badPoint in bad_point_index_2:
                ydata_2_u1[badPoint][0] =   mean_y2_u1[badPoint]
            
            for badPoint in bad_point_index_2:
                ydata_2_u9[badPoint][0] = mean_y2_u9[badPoint]
                              
            for nunNumber in range(0,  len(ydata_1_u1)):
                if  ydata_1_u1[nunNumber][0] == 0:
                    ydata_1_u1[nunNumber][0] = mean_y1_u1_2
                if  ydata_1_u9[nunNumber][0] == 0:
                    ydata_1_u9[nunNumber][0] = mean_y1_u9_2
                if  ydata_2_u1[nunNumber][0] == 0:
                    ydata_2_u1[nunNumber][0] = mean_y2_u1_2
                if  ydata_2_u9[nunNumber][0] == 0:
                    ydata_2_u9[nunNumber][0] = mean_y2_u9_2
           
            xdata_1_f = np.array(xdata_1_f)
            xdata_2_f = np.array(xdata_2_f)
            ydata_1_u1 = np.array(ydata_1_u1)
            ydata_2_u1 = np.array(ydata_2_u1)
            ydata_1_u9 = np.array(ydata_1_u9)
            ydata_2_u9 = np.array(ydata_2_u9)
            
            #data_1 = data_1[outliersMask_1]
            #data_2 = data_2[outliersMask_2]
            
            return (xdata_1_f, xdata_2_f, ydata_1_u1, ydata_2_u1, ydata_1_u9, ydata_2_u9)
        
        else:
            xdata_1_f = data1[:, [0]]
            xdata_2_f = data2[:, [0]]
            ydata_1_u1 = data1[:, [1]]
            ydata_2_u1 = data2[:, [1]]
            ydata_1_u9 = data1[:, [2]]
            ydata_2_u9 = data2[:, [2]]
            
            return (xdata_1_f, xdata_2_f, ydata_1_u1, ydata_2_u1, ydata_1_u9, ydata_2_u9)
        
    def nextPair(self):
        self.plot_start_u1.removePolt()
        self.plot_start_u9.removePolt()
        self.plotFrame.destroy()
        del self.plotFrame
        self.index = self.index + 1
        self.plotingInitialPairs(self.index)
        
    def plotingInitialPairs(self, index):
        self.window.title("After correlation")
        
        pair = self.scanPairs[index]
        self.plotFrame = frame(self.window,(1000,1000), LEFT)
        scanNUmber1 = self.dataFileDir + "/" + pair[0]
        scanNUmber2 = self.dataFileDir + "/" + pair[1]
            
        data_1 = np.fromfile(scanNUmber1, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber1),5))
        data_2 = np.fromfile(scanNUmber2, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber2),5))
        data_1 = np.delete(data_1, (0), axis=0) #izdzes masiva primo elementu
        data_2 = np.delete(data_2, (0), axis=0) #izdzes masiva primo elementu
            
        xdata_1_f, xdata_2_f, ydata_1_u1, ydata_2_u1, ydata_1_u9, ydata_2_u9 = self.__getDataForPolarization__(data_1, data_2, self.filter)
        self.plot_start_u1 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_start_u1.creatPlot(None, 'Frequency Mhz', 'Amplitude', "u1 Polarization")
        self.plot_start_u1.plot(xdata_1_f, ydata_1_u1, 'b', label=pair[0])
        self.plot_start_u1.plot(xdata_1_f, ydata_2_u1, 'r', label=pair[1])
            
        self.plot_start_u9 = Plot(5,5, self.masterFrame, self.plotFrame)
        self.plot_start_u9.creatPlot(None, 'Frequency Mhz', 'Amplitude', "u9 Polarization")
        self.plot_start_u9.plot(xdata_2_f, ydata_1_u9, 'b', label=pair[0])
        self.plot_start_u9.plot(xdata_2_f, ydata_2_u9, 'r', label=pair[1]) 
        
        if index == self.datPairsCount -1:
            self.nextPairButton.destroy()
        
    def __UI__(self):
        if self.index != self.datPairsCount -1:
            self.nextPairButton = Button (self.masterFrame, text="Next pair", command=self.nextPair, activebackground="Blue", background="Blue", font=self.font)
            self.nextPairButton.pack()
            
        self.plotingInitialPairs(self.index)
        
    def _quit(self):
        self.window.destroy()

def main():
    args = parseArguments()
    source = str(args.__dict__["source"])
    date = str(args.__dict__["date"])
    interval = float(args.__dict__["interval"])
    threshold = float(args.__dict__["threshold"])
    filter = str(args.__dict__["filter"])
    paircount = int(args.__dict__["paircount"])
    configFilePath = str(args.__dict__["config"])
    
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    dataFilesPath =  config.get('paths', "dataFilePath")
    badPointRange =  config.getint('parametrs', "badPointRange")
    
    if filter == "True" or filter == "true":
        filtering = True
    else:
        filtering = False
        
    if interval <= 0.0:
        raise Exception("Interval cannot be negative or zero")
    
    if threshold <= 0.0:
        raise Exception("Threshold cannot be negative or zero")
    
    if  paircount% 2 !=0:
        raise Exception("Paircount must Even be " + "got " + str(paircount))
    
    #Create App
    window = tk.Tk()
    window.configure(background='light goldenrod')
    ploting = Analyzer(window, source, date, filtering, threshold, badPointRange, dataFilesPath)
    img = tk.Image("photo", file="viraclogo.png")
    window.call('wm','iconphoto', window._w,img)
        
    ploting.mainloop()
    
    sys.exit(0)

if __name__=="__main__":
    main()