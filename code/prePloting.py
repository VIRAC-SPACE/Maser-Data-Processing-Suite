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
        
        plt.figure("polaration u1 after correlation")
        plt.plot(xdata_1_f, ydata_1_u1, label=pair[0])
        plt.plot(xdata_2_f, ydata_2_u1, label=pair[1])
        plt.legend(loc=2)
        plt.show()
        
        plt.figure("polaration u9 after correlation")
        plt.plot(xdata_1_f, ydata_1_u9, label=pair[0])
        plt.plot(xdata_2_f, ydata_2_u9, label=pair[1])
        plt.legend(loc=2)
        plt.show()
        
        
        plt.figure("polaration u1 first step")
        plt.plot(xdata_1_f, ydata_1_u1 - ydata_2_u1, label=pair[0] + " - " + pair[1])
        plt.legend(loc=2)
        plt.show()
        
        plt.figure("polaration u9 first step")
        plt.plot(xdata_1_f, ydata_1_u9 - ydata_2_u9, label=pair[0] + " - " + pair[1])
        plt.legend(loc=2)
        plt.show()
  
def main(source, date):
    scanPairs = createScanPairs(source, date)
    PlotScanPairs(scanPairs, source, date)
    sys.exit(0)
    
if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)
    
    source = sys.argv[1]
    date = sys.argv[2]    
    main(source, date)