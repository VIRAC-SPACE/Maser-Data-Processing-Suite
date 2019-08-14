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
    return [getIerattion(file) for file in [file for file in os.listdir(getConfigs("paths", path) + "/" + source) if file.startswith(source + "_")]]

def findLogFile(logList, iteration):
    for l in range(0, len(logList)):
        if logList[l].split("/")[-1].split(".")[0].split("_")[-1] == iteration:
            return l 
def main():
    smoohtIndexies = getFileIterations(getArgs("source"), "outputFilePath")
    notSmoohtIndexies = getFileIterations(getArgs("source"), "notSmoohtFilePath")
    hasNoNotSmoohtIndexies = [index for index in smoohtIndexies if index not in notSmoohtIndexies]
    dataFilePathForSourceSDR = getConfigs("paths", "dataFilePath") + "SDR/" + getArgs("source") + "/"
    dataFilePathForSourceDBBC = getConfigs("paths", "dataFilePath") + "DBBC/" + getArgs("source") + "/"
    dataFileIterationsSDR = [iteration for iteration in os.listdir(dataFilePathForSourceSDR)]
    dataFileIterationsDBBC = [iteration for iteration in os.listdir(dataFilePathForSourceDBBC)]

    logfile_listSDR = [log for log in os.listdir(getConfigs("paths", "logPath") + "/SDR") if log.startswith(getArgs("source"))]
    logfile_listDBBC = [log for log in os.listdir(getConfigs("paths", "logPath") + "/DBBC") if log.startswith(getArgs("source"))]

    for iteration in hasNoNotSmoohtIndexies:
        if iteration in dataFileIterationsSDR:
            frequencyShiftingParametr =  getArgs("source") + " " + getArgs("line") + " " + iteration + " " + str(logfile_listSDR[findLogFile(logfile_listSDR, iteration)])
            print("\033[1;31;47mExecute ","python3  " + "code/SDR_fs.py " + frequencyShiftingParametr + " \033[0;29;39m")
            os.system("python3 " + "code/SDR_fs.py " + frequencyShiftingParametr)
        else:
            print("Iteration", iteration, "do not exist in SDR data files")

        if iteration in dataFileIterationsDBBC:
            frequencyShiftingParametr = getArgs("source") + " " + iteration + " " + str(logfile_listDBBC[findLogFile(logfile_listDBBC, iteration)])
            print("\033[1;31;47mExecute ","python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr + " \033[0;29;39m")
            os.system("python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr)
        else:
            print("Iteration", iteration, "do not exist in DBBC data files")

    data_filesSDR = [d for d in os.listdir(getConfigs("paths", "dataFilePath") + "/SDR") if d.startswith(getArgs("source") + "_") and getIerattion(d) in hasNoNotSmoohtIndexies]
    data_filesDBBC = [d for d in os.listdir(getConfigs("paths", "dataFilePath") + "/DBBC") if d.startswith(getArgs("source") + "_") and getIerattion(d) in hasNoNotSmoohtIndexies]
                
    for d in data_filesSDR:
        if getArgs("noGUI"):
            print ("\033[1;31;47mExecute ",  "python3  " + "code/totalSpectrumAnalyer_qt5.py " + d  +  " -n \033[0;29;39m") 
            os.system("python3  " + "code/totalSpectrumAnalyer_qt5.py " + d + " -n")
        else:
            print ("\033[1;31;47mExecute ",  "python3  " + "code/totalSpectrumAnalyer_qt5.py " + d  +  " \033[0;29;39m") 
            os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + d)

    for d in data_filesDBBC:
        if getArgs("noGUI"):
            print("\033[1;31;47mExecute ", "python3  " + "code/totalSpectrumAnalyer_qt5.py " + d + " -n \033[0;29;39m")
            os.system("python3  " + "code/totalSpectrumAnalyer_qt5.py " + d + " -n")
        else:
            print("\033[1;31;47mExecute ", "python3  " + "code/totalSpectrumAnalyer_qt5.py " + d + " \033[0;29;39m")
            os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + d)
           
if __name__=="__main__":
    main()
    sys.exit(0)
    