#! /usr/bin/python
import sys
import os
import scipy
import scipy.stats
#from scipy.stats import signaltonoise
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas.stats.moments import rolling_mean
import argparse
import configparser
import re

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

def indexies(array, value):
    indexs = list()
    for i in range(0, len(array)-1):
        if array[i] == value:
            indexs.append(i)
    return indexs
   
def createScanPairs(source, date):
    dataFileDir = "dataFiles/" + source + "/" + date
    print (dataFileDir)
    
    dataFiles = list()
    for dataFile in os.listdir(dataFileDir):
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

def is_outlier(points, threshold):
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score < threshold

def PlotScanPairs(scanPairs, source, date, interval, threshold, filter, paircount, dataFilesPath, badPointRange):
    y_u1_results = list()
    y_u9_results = list()
    datPairsCount = len(scanPairs)
    
    filtering = True
    if str(filter) == "True":
        filtering = True
    else:
        filtering = False

    for pair in scanPairs:
        scanNUmber1 = dataFilesPath + source + "/" + date + "/" + pair[0]
        scanNUmber2 = dataFilesPath + source + "/" + date + "/" + pair[1]
        
        data_1 =  np.fromfile(scanNUmber1, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber1),5))
        data_2 =  np.fromfile(scanNUmber2, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber2),5))
        data_1 = np.delete(data_1, (0), axis=0) #izdzes masiva primo elementu
        data_2 = np.delete(data_2, (0), axis=0) #izdzes masiva primo elementu
        
        if filtering == True:
            outliersMask_1 = is_outlier(data_1, threshold)
            outliersMask_2 = is_outlier(data_2, threshold)
            
            bad_point_index_1 = indexies(outliersMask_1, False)
            bad_point_index_2 = indexies(outliersMask_2, False)
             
            xdata_1_f = data_1[:, [0]].tolist()
            xdata_2_f = data_2[:, [0]].tolist()
            ydata_1_u1 = data_1[:, [1]].tolist()
            ydata_2_u1 = data_2[:, [1]].tolist()
            ydata_1_u9 = data_1[:, [2]].tolist()
            ydata_2_u9 = data_2[:, [2]].tolist()
            
            df_y1_u1 = pd.DataFrame(data=ydata_1_u1)
            df_y1_u9 = pd.DataFrame(data=ydata_1_u9)
            df_y2_u1 = pd.DataFrame(data=ydata_2_u1)
            df_y2_u9 = pd.DataFrame(data=ydata_2_u9)
            
            mean_y1_u1 = np.nan_to_num(df_y1_u1.rolling(window=badPointRange, center=True).mean())
            mean_y1_u9 = np.nan_to_num(df_y1_u9.rolling(window=badPointRange, center=True).mean())
            mean_y2_u1 = np.nan_to_num(df_y2_u1.rolling(window=badPointRange, center=True).mean())
            mean_y2_u9 = np.nan_to_num(df_y2_u9.rolling(window=badPointRange, center=True).mean())
            
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
        
        else:
            xdata_1_f = data_1[:, [0]]
            xdata_2_f = data_2[:, [0]]
            ydata_1_u1 = data_1[:, [1]]
            ydata_2_u1 = data_2[:, [1]]
            ydata_1_u9 = data_1[:, [2]]
            ydata_2_u9 = data_2[:, [2]]
            
        plt.figure("polarization u1 after correlation")
        plt.plot(xdata_1_f, ydata_1_u1, label=pair[0] )
        plt.plot(xdata_2_f, ydata_2_u1, label=pair[1]  )
        plt.legend(loc=2)
        plt.grid(True)
        plt.show()
        
        plt.figure("polarization u9 after correlation")
        plt.plot(xdata_1_f, ydata_1_u9, label=pair[0])
        plt.plot(xdata_2_f, ydata_2_u9, label=pair[1])
        plt.legend(loc=2)
        plt.grid(True)
        plt.show()
        
        data_y_u1 = ydata_1_u1 - ydata_2_u1
        data_y_u9 = ydata_1_u9 - ydata_2_u9
        
        plt.figure("polarization u1 first step")
        plt.plot(xdata_1_f, data_y_u1, label=pair[0] + " - " + pair[1])
        plt.legend(loc=2)
        plt.grid(True)
        plt.show()
        
        plt.figure("polarization u9 first step")
        plt.plot(xdata_1_f, data_y_u9, label=pair[0] + " - " + pair[1])
        plt.legend(loc=2)
        plt.grid(True)
        plt.show()
        
        maxFrequency = np.max(xdata_1_f)
        if interval > maxFrequency/4.0:
            raise Exception("Interval cannot be larger than 0.25 of frequency range ", "max interval " + str(maxFrequency/4.0))
        frecquencyRange_1 = (maxFrequency/4.0 - interval, maxFrequency/4.0  + interval) #Negative range
        frecquencyRange_2 = (maxFrequency*(3.0/4.0) - interval, maxFrequency*(3.0/4.0) + interval) #positive range
        
        #Creating index
        index_1_1 = (np.abs(xdata_1_f-frecquencyRange_1[0])).argmin()
        index_1_2 = (np.abs(xdata_1_f-frecquencyRange_1[1])).argmin() 

        index_2_1 = (np.abs(xdata_1_f-frecquencyRange_2[0])).argmin() 
        index_2_2 = (np.abs(xdata_1_f-frecquencyRange_2[1])).argmin()

        #check indexies
        if index_2_2 - index_2_1!= index_1_2 - index_1_1:
            print ("befor correction", index_2_2 - index_2_1 + 1,  index_1_2 - index_1_1 + 1, [index_1_1, index_1_2], [index_2_1, index_2_2])
            if index_2_2 - index_2_1 + 1 > index_1_2:
                index = np.abs(index_2_2 - index_2_1 + 1 - index_1_2) 
                index_1_1 = (np.abs(xdata_1_f-frecquencyRange_1[0])).argmin()
                index_1_2 = (np.abs(xdata_1_f-frecquencyRange_1[1])).argmin() -1
                index_2_1 = (np.abs(xdata_1_f-frecquencyRange_2[0])).argmin()  + index
                index_2_2 = (np.abs(xdata_1_f-frecquencyRange_2[1])).argmin() 
                
            elif index_2_2 - index_2_1 + 1 < index_1_2:
                index = np.abs(index_2_2 - index_2_1 + 1 - index_1_2)
                index_1_1 = (np.abs(xdata_1_f-frecquencyRange_1[0])).argmin()
                index_1_2 = (np.abs(xdata_1_f-frecquencyRange_1[1])).argmin() +1
                index_2_1 = (np.abs(xdata_1_f-frecquencyRange_2[0])).argmin() -index
                index_2_2 = (np.abs(xdata_1_f-frecquencyRange_2[1])).argmin()
                
            print ("after correction", index_2_2 - index_2_1,  index_1_2 - index_1_1, [index_1_1, index_1_2], [index_2_1, index_2_2])
            
        else:
            print ("indexies correct", index_2_2 - index_2_1,  index_1_2 - index_1_1, [index_1_1, index_1_2], [index_2_1, index_2_2])
          
        negativeRange_u1 = data_y_u1[index_1_1:index_1_2]
        positiveveRange_u1 = data_y_u1[index_2_1:index_2_2]
        
        negativeRange_u9 = data_y_u9[index_1_1:index_1_2]
        positiveveRange_u9 = data_y_u9[index_2_1:index_2_2]
                
        result_u1 = (positiveveRange_u1 - negativeRange_u1)/2
        
        print ("total point count ", len( result_u1 ))
        result_u9 = (positiveveRange_u9 - negativeRange_u9)/2
        
        x = np.linspace(0,maxFrequency/2, len(result_u1), dtype="float64").reshape(len(result_u1), 1)
        y_u1_results.append(result_u1)
        y_u9_results.append(result_u9)
        
        #ston = signaltonoise(result_u1)
        #print ("signal vs noise ", ston)
        
        plt.figure("polarization u1 second step")
        plt.plot(x, result_u1)
        plt.grid(True)
        plt.show()
        
        plt.figure("polarization u9 second step")
        plt.plot(x, result_u9)
        plt.grid(True)
        plt.show()
     
    y_u1_avg = np.zeros(y_u1_results[0].shape)
    y_u2_avg = np.zeros(y_u9_results[0].shape)
    dummyData = np.zeros(len(result_u1), dtype="float64").reshape(len(result_u1), 1)
    
    if  paircount == 0:
        i = 0
        # creating average for day
        for result in y_u1_results:
            y_u1_avg = y_u1_avg + result
                
        for result2 in y_u9_results:
            y_u2_avg = y_u2_avg + result2
                
        y_u1_avg = np.array(y_u1_avg/datPairsCount, dtype="float64")
        y_u2_avg = np.array(y_u2_avg/datPairsCount, dtype="float64")
            
        plt.figure("polarization u1 average for experiment")
        plt.plot(x, y_u1_avg)
        plt.grid(True)
        plt.show()
            
        plt.figure("polarization u9 average for experiment")
        plt.plot(x, y_u2_avg)
        plt.grid(True)
        plt.show()
        
        scan_number = int(re.split("([0-9]+)", scanPairs[i][0])[-2])
        totalResults = np.concatenate((x, y_u1_avg, y_u2_avg, dummyData, dummyData), axis=1)
        np.savetxt(dataFilesPath + source + date.replace(".", "_")  + "_n" + str(scan_number) + ".dat", totalResults)
    
    elif paircount == 2:
        
        i = 0
        for output in range(0, len(y_u1_results)):
            scan_number = int(re.split("([0-9]+)", scanPairs[i][0])[-2])
            totalResults = np.concatenate((x, y_u1_results[output], y_u9_results[output], dummyData, dummyData), axis=1)
            np.savetxt(dataFilesPath + source + date.replace(".", "_")  +"_k_" + str(output) + "_n" + str(scan_number) + ".dat", totalResults)
            i = i + 1   
    else:
        pairNumber = 0
        check = 0
        
        i = 0
        for result in range(0, len(y_u1_results)):
            pairNumber = pairNumber + 2
            y_u1_avg = y_u1_avg + y_u1_results[result]
            y_u2_avg = y_u2_avg + y_u9_results[result]
            
            if pairNumber == paircount:
                
                check = check + 1
                y_u1_avg = np.array(y_u1_avg/datPairsCount, dtype="float64")
                y_u2_avg = np.array(y_u2_avg/datPairsCount, dtype="float64")
                    
                plt.figure("polarization u1 average for pair count")
                plt.plot(x, y_u1_avg)
                plt.grid(True)
                plt.show()
                    
                plt.figure("polarization u9 average for pair count")
                plt.plot(x, y_u2_avg)
                plt.grid(True)
                plt.show()
                
                scan_number = int(re.split("([0-9]+)", scanPairs[i][0])[-2])
                totalResults = np.concatenate((x, y_u1_avg, y_u2_avg, dummyData, dummyData), axis=1)
                np.savetxt(dataFilesPath + source + date.replace(".", "_")  +"_k_" + str(result) + "_n" + str(scan_number) + ".dat", totalResults)
                y_u1_avg = np.zeros(y_u1_results[0].shape)
                y_u2_avg = np.zeros(y_u9_results[0].shape)
                
                pairNumber = 0
            i = i + 1       
                
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
    
    if interval <= 0.0:
        raise Exception("Interval cannot be negative or zero")
    
    if threshold <= 0.0:
        raise Exception("Threshold cannot be negative or zero")
    
    if  paircount% 2 !=0:
        raise Exception("Paircount must Even be " + "got " + str(paircount))
    
    scanPairs = createScanPairs(source, date)
    PlotScanPairs(scanPairs, source, date, interval, threshold, filter, paircount, dataFilesPath, badPointRange)
    sys.exit(0)
    
if __name__=="__main__":
    main()
