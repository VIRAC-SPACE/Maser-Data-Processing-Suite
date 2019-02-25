#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
from parsers._configparser import ConfigParser

def parseArguments():
    parser = argparse.ArgumentParser(description='''Create output file. ''', epilog="""OutputFileCreator.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 1.0')
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-n", "--noGUI", help="Create smoothed and not smothed outputfiles", action='store_true')
    args = parser.parse_args()
    return args

def getArgs(key):
    return str(parseArguments().__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

def getIerattion(file):
    return file.split(".")[0].split("_")[-1]
    
def getFileIterations(source, path):
    return [getIerattion(file) for file in [file for file in os.listdir(getConfigs("paths", path)) if file.startswith(source + "_")]]

def findLogFile(logList, iteration):
    for l in range(0, len(logList)):
        if logList[l].split("/")[-1].split(".")[0].split("_")[-1] == iteration:
            return l 
def main():
    smoohtIndexies = getFileIterations(getArgs("source"), "outputFilePath")
    notSmoohtIndexies = getFileIterations(getArgs("source"), "notSmoohtFilePath")
    hasNoNotSmoohtIndexies = [index for index in smoohtIndexies if index not in notSmoohtIndexies]
    dataFilePathForSource = getConfigs("paths", "dataFilePath") + getArgs("source") + "/"
    dataFileIterations = [iteration for iteration in os.listdir(dataFilePathForSource)]
    logfile_list = [log for log in os.listdir(getConfigs("paths", "logPath")) if log.startswith(getArgs("source"))]
        
    for iteration in hasNoNotSmoohtIndexies:
        if iteration not in dataFileIterations:
            print("Iteration", iteration, "do not exist in data files")
        else:
            frequencyShiftingParametr = getArgs("source") + " " + iteration + " " + str(logfile_list[findLogFile(logfile_list, iteration)])
            print ("\033[1;31;47mExecute ",  "python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr +  " \033[0;29;39m")
            os.system("python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr)
            
    data_files = [d for d in os.listdir(getConfigs("paths", "dataFilePath")) if d.startswith(getArgs("source")) and getIerattion(d) in hasNoNotSmoohtIndexies]
                
    for d in data_files:
        if getArgs("noGUI"):
            print ("\033[1;31;47mExecute ",  "python3  " + "code/totalSpectrumAnalyer_qt5.py " + d  +  " -n \033[0;29;39m") 
            os.system("python3  " + "code/totalSpectrumAnalyer_qt5.py " + d + " -n")
        else:
            print ("\033[1;31;47mExecute ",  "python3  " + "code/totalSpectrumAnalyer_qt5.py " + d  +  " \033[0;29;39m") 
            os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + d)
           
if __name__=="__main__":
    main()
    sys.exit(0)
    