#! /usr/bin/python
import sys
import os
import argparse
import configparser
from tkinter import *
import tkinter as tk
import numpy as np

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

class Analyzer():
    def __init__(self, window, source=None, date=None, datafile):
        self.window = window
        self.source = source
        self.date = date
        
        try:
            data = np.fromfile(datafile, dtype="float64", count=-1, sep=" ") .reshape((file_len(datafile),3))
        
        except IOError as e:
            print ("IO Error",  e)
            sys.exit(1)
                
        except:
            print("Unexpected error:", sys.exc_info()[0])
            sys.exit(1)
                 
        else:
            data = np.delete(data, (0), axis=0) #izdzes masiva primo elementu
            self.dataPoints = data.shape[0]
            
            self.xdata = data[:, [0]]
            self.y_u1 = data[:, [1]] 
            self.y_u9 = data[:, [2]]

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
    datafile = list(args.__dict__["datafile"])
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
    ploting = Analyzer(window, source=None, date=None, dataFilesPath + datafile)
    img = tk.Image("photo", file="viraclogo.png")
    window.call('wm','iconphoto', window._w,img)
        
    ploting.mainloop()
        
    sys.exit(0)

if __name__=="__main__":
    main()