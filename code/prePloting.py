#! /usr/bin/python
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

def usage():
    print ('Usage: source and date')
    
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
    
def createScanPairs(source, date):
    dataFileDir = "dataFiles/" + source + "/" + date
    
    dataFiles = list()
    for dataFile in os.listdir(dataFileDir):
        dataFiles.append(dataFile)
    
    dataFiles.sort()
            
    scanPairs = list()
    i = 0
    j = 1
    for k in range(0, len(dataFiles) - 3):
        scanPairs.append((dataFiles[i], dataFiles[j]))
        i = i + 2
        j = j + 2
    
    return scanPairs

def PlotScanPairs(scanPairs, source, date):
    y_u1_results = list()
    y_u9_results = list()
    datPairsCount = len(scanPairs)
    
    for pair in scanPairs:
        scanNUmber1 = "dataFiles/" + source + "/" + date + "/" + pair[0]
        scanNUmber2 = "dataFiles/" + source + "/" + date + "/" + pair[1]
        
        data_1 =  np.fromfile(scanNUmber1, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber1),5))
        data_2 =  np.fromfile(scanNUmber2, dtype="float64", count=-1, sep=" ") .reshape((file_len(scanNUmber2),5))
        data_1 = np.delete(data_1, (0), axis=0) #izdzes masiva primo elementu
        data_2 = np.delete(data_2, (0), axis=0) #izdzes masiva primo elementu
        
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
        frecquencyRange_1 = (maxFrequency/4.0 - 0.5, maxFrequency/4.0 + 0.5) #Negative range
        frecquencyRange_2 = (maxFrequency*(3.0/4.0) - 0.5, maxFrequency*(3.0/4.0) + 0.5) #positive range
        
        #Creating index
        index_1_1 = (np.abs(xdata_1_f-frecquencyRange_1[0])).argmin()
        index_1_2 = (np.abs(xdata_1_f-frecquencyRange_1[1])).argmin()
        index_2_1 = (np.abs(xdata_1_f-frecquencyRange_2[0])).argmin()
        index_2_2 = (np.abs(xdata_1_f-frecquencyRange_2[1])).argmin()
        
        negativeRange_u1 = data_y_u1[index_1_1:index_1_2]
        positiveveRange_u1 = data_y_u1[index_2_1:index_2_2]
        
        negativeRange_u9 = data_y_u9[index_1_1:index_1_2]
        positiveveRange_u9 = data_y_u9[index_2_1:index_2_2]
        
        result_u1 = (positiveveRange_u1 - negativeRange_u1)/2
        result_u9 = (positiveveRange_u9 - negativeRange_u9)/2
        
        x = np.linspace(0,maxFrequency/2, len(result_u1), dtype="float64").reshape(len(result_u1), 1)
        y_u1_results.append(result_u1)
        y_u9_results.append(result_u9) 
        
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
    
    dummyData = np.zeros(len(result_u1), dtype="float64").reshape(len(result_u1), 1)
    totalResults = np.concatenate((x, y_u1_avg, y_u2_avg, dummyData, dummyData), axis=1)
    np.savetxt("dataFiles/" + source + date.replace(".", "_")  +"_n1.dat", totalResults)
    
    return totalResults

def main(source, date):
    scanPairs = createScanPairs(source, date)
    cordata  = PlotScanPairs(scanPairs, source, date)
    logFile = source + ".log"
    #maserPloting(cordata, logFile)
    
    sys.exit(0)
    
if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    
    source = sys.argv[1]
    date = sys.argv[2]    
    main(source, date)