#! /usr/bin/python3
import sys
import os
import argparse
import configparser
from tkinter import *
import tkinter as tk
from tkinter import font
from tkinter import simpledialog
import numpy as np
import scipy.constants
import pandas as pd
from pandas.stats.moments import rolling_mean
from time import strptime
import json

from experimentsLogReader import ExperimentLogReader
from ploting import Plot

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Creates input file for plotting tool. ''', epilog="""PRE PLOTTER.""")
    
    # Positional mandatory arguments
    parser.add_argument("source", help="Experiment source", type=str)
    parser.add_argument("iteration_number", help="iteration number ", type=int)
    parser.add_argument("logFile", help="Experiment log file name", type=str)

    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration Yaml file", type=str, default="config/config.cfg")
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
    f.pack(side = sides, expand=YES)
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

def dopler(ObservedFrequency, velocityReceiver, f0):
    c = scipy.constants.speed_of_light
    #f0 = 6668519200 # Hz 
    velocitySoure = (-((ObservedFrequency/f0)-1)*c + (velocityReceiver * 1000))/1000
    return velocitySoure

def STON(array):
    std = np.std(array) 
    max = np.max(array)
    
    ston = max/(std*3)
    return ston

class Analyzer(Frame):
    def __init__(self, window, source, iteration_number, filter, threshold, badPointRange, dataPath, resultPath, logs):
        Frame.__init__(self)
        self.window = window
        self.source = source
        self.threshold = threshold
        self.filter = filter
        self.badPointRange = badPointRange
        self.dataFilesPath = dataPath
        self.resultPath = resultPath
        self.font = font.Font(family="Times New Roman", size=20, weight=font.BOLD)
        self.font_2 = font.Font(family="Times New Roman", size=10)
        self.masterFrame = frame(self.window,(1000,1000), RIGHT)
        self.index = 0
        self.totalResults_u1 = list()
        self.totalResults_u9 = list()
        self.STON_list_u1 = list()
        self.STON_list_u9 = list()
        self.STON_list_AVG = list()
        self.iteration_number = iteration_number
        self.logs = logs
        self.date = self.logs["header"]["dates"]
        self.dataFileDir = dataPath + self.source + "/" + str(self.iteration_number) + "/"
        self.scanPairs = self.createScanPairs()
        self.datPairsCount = len(self.scanPairs)
        self.f0 = 6668519200
        self.location = self.logs["location"]
        self.calibrationScale = self.calibrationScales[self.location]
        self.expername = self.source + self.date + "_" + self.logs["location"]
    
        self.window.title("Analyze for " + self.source + " " + self.date)
        self.__UI__()
        
        
    def createScanPairs(self):
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
             
            xdata = data1[:, [0]].tolist()
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
           
            xdata = np.array(xdata)
            ydata_1_u1 = np.array(ydata_1_u1)
            ydata_2_u1 = np.array(ydata_2_u1)
            ydata_1_u9 = np.array(ydata_1_u9)
            ydata_2_u9 = np.array(ydata_2_u9)
            
            self.dataPoints = len(xdata)
            
            return (xdata, ydata_1_u1, ydata_2_u1, ydata_1_u9, ydata_2_u9)
        
        else:
            xdata = data1[:, [0]]
            ydata_1_u1 = data1[:, [1]]
            ydata_2_u1 = data2[:, [1]]
            ydata_1_u9 = data1[:, [2]]
            ydata_2_u9 = data2[:, [2]]
            
            self.dataPoints = len(xdata)
            
            return (xdata, ydata_1_u1, ydata_2_u1, ydata_1_u9, ydata_2_u9)
    
    def calibration(self, array_x, data_1, data_2, tsys_1, tsys_2, elevation):
        #from AGN cal sessions (FS /usr2/control/rxg_files/c3.rxg):
        DPFU_max = [0.0442016750, 0.0444381686]
        G_El = [-0.0000740207, 0.0071794464, 0.8243803753]
        Tcal = 3.859
        k = 0.851
        DPFU = np.mean(DPFU_max)*np.polyval(G_El, elevation)
        f_shift = 0.5
        
        P_sig = data_1 # Get Amplitudes
        P_ref = data_2 # Get Amplitudes
        
        Ta_sig = float(tsys_1)*(-P_sig + P_ref)/P_ref #only non-cal phase for dbbc possible...
        Ta_ref = float(tsys_2)*(P_ref - P_sig)/P_sig
        
        f_step = (array_x[self.dataPoints-1]-array_x[0])/(self.dataPoints-1); 
        n_shift = int(f_shift/f_step);
        
        Ta_sig = np.roll(Ta_sig, -n_shift); # pos
        Ta_ref = np.roll(Ta_ref, -n_shift); # neg
        
        #avg shifted spectrums
        Ta = (Ta_sig + Ta_ref)/2 # Creting total spectr
        
        #K->Jy
        Ta = Ta/DPFU/k
  
        return Ta
    
    def nextPair(self):
        self.plot_start_u1.removePolt()
        self.plot_start_u9.removePolt()
        self.plotFrame_start.destroy()
        del self.plotFrame_start
        self.plot_negative_positive_u1.removePolt()
        self.plot_negative_positive_u9.removePolt()
        self.plotFrame_negative_positive.destroy()
        del self.plotFrame_negative_positive
        self.plot_total_u1.removePolt()
        self.plot_total_u9.removePolt()
        self.plotFrame_total.destroy()
        del self.plotFrame_total
        self.index = self.index + 1
        self.plotingPairs(self.index)
        
    def plotingPairs(self, index):
    
        pair = self.scanPairs[index]
        self.plotFrame_start = frame(self.window,(1000, 1000), TOP)
        self.plotFrame_total = frame(self.window,(1000, 1000), BOTTOM)
        self.plotFrame_negative_positive = frame(self.window,(1000, 1000), BOTTOM)
        
        scanNUmber1 = self.dataFileDir + "/" + pair[0]
        scanNUmber2 = self.dataFileDir + "/" + pair[1]
        
        print ("data files ", scanNUmber1, scanNUmber2)
        
        scan_number_1 = pair[0].split(".")[0].split("_")[-1][2:].lstrip("0")
        scan_number_2 = pair[1].split(".")[0].split("_")[-1][2:].lstrip("0")
        
        print ("scan number", scan_number_1, scan_number_2)
        
        scan_1 = self.logs[str(scan_number_2)]
        scan_2 = self.logs[str(scan_number_1)]
        
        # get system temperature
        tsys_u1_1 = scan_1['Systemtemperature'][0]
        tsys_u1_2 = tsys_u1_1
        tsys_u9_1 = scan_1['Systemtemperature'][1]
        tsys_u9_2 = tsys_u9_1
        
        elevation = (float(scan_1["elevation"]) + float(scan_2["elevation"])) /2
        
        print ("tsys", tsys_u1_1, tsys_u9_1)
        
        if float(tsys_u1_1) == 0:
            newT = simpledialog.askfloat("System temperature is zero",  " Expected number between 0 and 300", minvalue = 1, maxvalue = 300)
            tsys_u1_1 = newT
            
        if float(tsys_u9_1) == 0:
            newT = simpledialog.askfloat("System temperature is zero",  " Expected number between 0 and 300", minvalue = 1, maxvalue = 300)
            tsys_u9_1 = newT
                
        data_1 = np.fromfile(scanNUmber1, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber1),9))
        data_2 = np.fromfile(scanNUmber2, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber2),9))
        
        #Delete first row
        data_1 = np.delete(data_1, (0), axis=0) #izdzes masiva primo elementu
        data_2 = np.delete(data_2, (0), axis=0) #izdzes masiva primo elementu
            
        xdata, ydata_1_u1, ydata_2_u1, ydata_1_u9, ydata_2_u9 = self.__getDataForPolarization__(data_1, data_2, self.filter)
           
        self.plot_start_u1 = Plot(4,4, self.masterFrame, self.plotFrame_start)
        self.plot_start_u1.creatPlot(None, 'Frequency Mhz', 'Amplitude', "u1 Polarization")
        self.plot_start_u1.plot(xdata, ydata_1_u1, 'b', label=pair[0])
        self.plot_start_u1.plot(xdata, ydata_2_u1, 'r', label=pair[1])
            
        self.plot_start_u9 = Plot(4,4, self.masterFrame, self.plotFrame_start)
        self.plot_start_u9.creatPlot(None, 'Frequency Mhz', 'Amplitude', "u9 Polarization")
        self.plot_start_u9.plot(xdata, ydata_1_u9, 'b', label=pair[0])
        self.plot_start_u9.plot(xdata, ydata_2_u9, 'r', label=pair[1])
        
        #Calibration  
        data_u1 = self.calibration(xdata, ydata_1_u1, ydata_2_u1, float(tsys_u1_1), float(tsys_u1_2), elevation) 
        data_u9 = self.calibration(xdata, ydata_1_u9, ydata_2_u9, float(tsys_u9_1), float(tsys_u9_2), elevation)
       
        xdata = np.array(xdata)
        
        self.plot_negative_positive_u1 = Plot(4,4, self.masterFrame, self.plotFrame_negative_positive)
        self.plot_negative_positive_u1.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', None)
        self.plot_negative_positive_u1.plot(xdata, data_u1, 'b', label=pair[0] +  "-" + pair[1])
        
        self.plot_negative_positive_u9 = Plot(4,4, self.masterFrame, self.plotFrame_negative_positive)
        self.plot_negative_positive_u9.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', None)
        self.plot_negative_positive_u9.plot(xdata, data_u9, 'b', label=pair[0] +  "-" + pair[1])
        
        self.x = xdata
        f_step = (self.x[self.dataPoints-1]-self.x[0])/(self.dataPoints-1) 
        f_shift = 0.5
        n_shift = int(f_shift/f_step)
        total_u1 = data_u1[(n_shift+1):(self.dataPoints - n_shift - 1)]
        total_u9 = data_u9[(n_shift+1):(self.dataPoints - n_shift - 1)]
        
        self.x = self.x[(n_shift+1):(self.dataPoints - n_shift - 1)] 
        
        self.totalResults_u1.append(total_u1)
        self.totalResults_u9.append(total_u9)
        
        self.plot_total_u1 = Plot(4,4, self.masterFrame, self.plotFrame_total)
        self.plot_total_u1.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', None)
        self.plot_total_u1.plot(self.x, total_u1, 'b')
        
        self.plot_total_u9 = Plot(4,4, self.masterFrame, self.plotFrame_total)
        self.plot_total_u9.creatPlot(None, 'Frequency Mhz', 'Flux density (Jy)', None)
        self.plot_total_u9.plot(self.x, total_u9, 'b')
        
        ston_u1 = STON(total_u1)
        ston_u9 = STON(total_u9)
        stone_AVG = STON(((total_u1 + total_u9)/2))
        
        self.STON_list_u1.append(ston_u1)
        self.STON_list_u9.append(ston_u9)
        self.STON_list_AVG.append(stone_AVG)
        
        if index == self.datPairsCount -1:
            self.nextPairButton.destroy()
            self.totalResultButton = Button (self.masterFrame, text="Move to total results", command=self.plotTotalResults, activebackground="Blue", background="Blue", font=self.font)
            self.totalResultButton.pack()
            
    def plotTotalResults(self):
        self.plot_start_u1.removePolt()
        self.plot_start_u9.removePolt()
        self.plot_negative_positive_u1.removePolt()
        self.plot_negative_positive_u9.removePolt()
        self.plot_total_u1.removePolt()
        self.plot_total_u9.removePolt()
        del self.plot_start_u1
        del self.plot_start_u9
        del self.plot_negative_positive_u1
        del self.plot_negative_positive_u9
        del self.plot_total_u1
        del self.plot_total_u9
        self.plotFrame_start.destroy()
        self.plotFrame_total.destroy()
        self.plotFrame_negative_positive.destroy()
        del self.plotFrame_start
        del self.plotFrame_total
        del self.plotFrame_negative_positive
        self.totalResultButton.destroy()
        del self.totalResultButton
        
        velocitys_avg = np.zeros(self.totalResults_u1[0].shape)
        y_u1_avg = np.zeros(self.totalResults_u1[0].shape)
        y_u9_avg = np.zeros(self.totalResults_u9[0].shape)
        
        for p in range(0,  self.datPairsCount):
            scan_number_1 = self.scanPairs[p][0].split("_")[-1][2:].lstrip("0").split(".")[0]
            scan_number_2 = self.scanPairs[p][1].split("_")[-1][2:].lstrip("0").split(".")[0]
            print ("pairs ", self.scanPairs[p])
            scan_1 = self.logs[str(scan_number_1)]
            scan_2 = self.logs[str(scan_number_2)]
            
            timeStr = scan_1['startTime'].replace(":", " ")
            dateStrList = scan_1['dates'].split()
            dateStrList[1] = strptime(dateStrList[1],'%b').tm_mon
            dateStr = str(dateStrList[2]) + " " + str(dateStrList[1]) + " " + str(dateStrList[0])
            RaStr = " ".join(scan_1["Ra"])
            DecStr = " ".join(scan_1["Dec"])
            FreqStart = (float(scan_1["fs_frequencyfs"])  +   float(scan_2["fs_frequencyfs"]))/2  + float(self.logs["header"]["BBC"])
            print ("FreqStart", FreqStart, "log frecq ", float(scan_1["fs_frequencyfs"]), float(scan_2["fs_frequencyfs"]), "BBC", float(self.logs["header"]["BBC"]))
            dopsetPar = dateStr + " " + timeStr + " " + RaStr + " " + DecStr
            print ("dopsetPar", dopsetPar,  " dateStr ", dateStr + " timeStr " + timeStr + " RaStr " + RaStr + " DecStr" + DecStr)
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
            
            velocitys = dopler((self.x + FreqStart) * (10 ** 6), VelTotal, self.f0)
            y_u1_avg =  y_u1_avg + self.totalResults_u1[p]
            y_u9_avg =  y_u9_avg + self.totalResults_u9[p]
            velocitys_avg = velocitys_avg + velocitys
            
        velocitys_avg =  velocitys_avg/len(self.totalResults_u1)
        y_u1_avg = y_u1_avg/len(self.totalResults_u1)
        y_u9_avg = y_u9_avg/len(self.totalResults_u9)
        
        self.plotFrame_velocity = frame(self.window,(1000, 1000), TOP)
        self.plotFrame_STON = frame(self.window,(1000, 1000), BOTTOM)
        
        # testing data
        #test_file = "/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/dataFiles/cepa_2018-05-07_m193_n02.dat.out"
        #test_data = np.fromfile(test_file, dtype="float64", count=-1, sep=" ") .reshape((file_len(test_file),4))
        #x = test_data[:, [0]].tolist()
        #y = test_data[:, [1]].tolist()
        
        self.plot_velocity_u1 = Plot(5,5, self.masterFrame, self.plotFrame_velocity)
        self.plot_velocity_u1.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "u1 Polarization")
        self.plot_velocity_u1.plot(velocitys_avg, y_u1_avg, 'b')
        #self.plot_velocity_u1.plot(x, y, 'r')
        
        self.plot_velocity_u9 = Plot(5,5, self.masterFrame, self.plotFrame_velocity)
        self.plot_velocity_u9.creatPlot(None, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "u9 Polarization")
        self.plot_velocity_u9.plot(velocitys_avg, y_u9_avg, 'b')
        
        ston_x = np.arange(0, len(self.STON_list_u1))
        self.plot_STON = Plot(5,5, self.masterFrame, self.plotFrame_STON)
        self.plot_STON.creatPlot(None, 'Pair', 'Ratio', "Signal to Noise")
        self.plot_STON.plot(ston_x, self.STON_list_u1, '*r', label="u1 Polarization")
        self.plot_STON.plot(ston_x, self.STON_list_u9, 'og', label="u9 Polarization")
        self.plot_STON.plot(ston_x, self.STON_list_AVG, 'vb', label="AVG Polarization")
        
        totalResults = np.concatenate((velocitys_avg, y_u1_avg, y_u9_avg), axis=1)
        output_file_name = self.dataFilesPath + self.source + self.date.replace(".", "_") + "_" + self.logs["location"] + ".dat"
        output_file_name = output_file_name.replace(" ", "")
        np.savetxt(output_file_name, totalResults)
        
        resultFile = str(self.resultPath) + str(self.source) + ".json"
        
        if os.path.isfile(resultFile):
            pass
        else:
            os.system("touch " + resultFile)
            
            resultFile = open (resultFile, "w")
            resultFile.write("{ \n" + "\n}")
            resultFile.close()
        
        with open(resultFile) as result_data:    
            result = json.load(result_data)
        
        if self.expername not in result:
            result[self.expername] = dict()
        
    def __UI__(self):
        if self.index != self.datPairsCount -1: # cheking if there is not one pair
            self.nextPairButton = Button (self.masterFrame, text="Next pair", command=self.nextPair, activebackground="Blue", background="Blue", font=self.font)
            self.nextPairButton.pack()
            
        self.plotingPairs(self.index)
        
    def _quit(self):
        self.window.destroy()
    
def main():
    args = parseArguments()
    source = str(args.__dict__["source"])
    iteration_number = int(args.__dict__["iteration_number"])
    logFile = str(args.__dict__["logFile"])
    threshold = float(args.__dict__["threshold"])
    filter = str(args.__dict__["filter"])
    singleSourceExperiment = list(args.__dict__["single"])
    configFilePath = str(args.__dict__["config"])
    
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    dataFilesPath =  config.get('paths', "dataFilePath")
    prettyLogsPath =  config.get('paths', "prettyLogsPath")
    logPath = config.get('paths', "logPath")
    resultPath = config.get('paths', "resultFilePath")
    badPointRange =  config.getint('parametrs', "badPointRange")
    logs  = ExperimentLogReader(logPath + logFile, prettyLogsPath, singleSourceExperiment).getLogs()
    
    if filter == "True" or filter == "true":
        filtering = True
    else:
        filtering = False
        
    if threshold <= 0.0:
        raise Exception("Threshold cannot be negative or zero")   
    
    #Create App
    window = tk.Tk()
    window.configure(background='light goldenrod')
    ploting = Analyzer(window, source, iteration_number, filtering, threshold, badPointRange, dataFilesPath, resultPath, logs)
    img = tk.Image("photo", file="viraclogo.png")
    window.call('wm','iconphoto', window._w,img)
        
    ploting.mainloop()
    
    sys.exit(0)

if __name__=="__main__":
    main()
    