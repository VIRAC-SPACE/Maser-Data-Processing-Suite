#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import configparser
import json

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''automatically call frequencyShiftingAnalyzer and totalSpectrumAnalyer. ''',
    epilog="""Main program.""")

    # Positional mandatory arguments
    parser.add_argument("source", help="Source Name", type=str)
    
    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 3.0')

    # Parse arguments
    args = parser.parse_args()

    return args

def findLogFile(logList, iteration):
    for l in range(0, len(logList)):
        if logList[l].endswith(iteration + ".log"):
            return l
        
def main():
    # Parse the arguments
    args = parseArguments()
    sourceName = str(args.__dict__["source"])
    configFilePath = str(args.__dict__["config"])
    
    #Creating config parametrs
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    dataFilesPath = config.get('paths', "dataFilePath")
    resultPath = config.get('paths', "resultFilePath")
    logPath = config.get('paths', "logPath")
    
    path = dataFilesPath + sourceName + "/"
    iterations = list()
     
    for iteration in os.listdir(path):
        iterations.append(iteration)
            
    iterations.sort(key=int, reverse=False)
     
    logfile_list = list()
    
    for log in os.listdir(logPath):
        if log.startswith(sourceName):
            logfile_list.append(log)
            
    resultFileName = sourceName + ".json"
        
    if os.path.isfile(resultPath + resultFileName):
        pass
    else:
        os.system("touch " + resultPath +  resultFileName)
            
        resultFile = open (resultPath +  resultFileName, "w")
        resultFile.write("{ \n" + "\n}")
        resultFile.close()
        
    with open(resultPath + resultFileName) as result_data:    
        result = json.load(result_data)
    
    processed_iteration = list()
      
    for experiment in result:
        if experiment[-1]  in iterations:
            processed_iteration.append(experiment[-1])
    
    for i in range(0, len(iterations)):
        
        if i not in processed_iteration:
            frequencyShiftingParametr = sourceName + " " + iterations[i] + " " + str(logfile_list[findLogFile(logfile_list, iterations[i])])
            print ("Execute ",  "python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr)
            os.system("python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr)
           
if __name__=="__main__":
    main()
    