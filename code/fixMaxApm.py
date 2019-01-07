#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os
import json
import argparse
import configparser
import numpy as np

from help import *

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Fix Amplitudes in output file. ''',
    epilog="""Fix.""")
    
    # Positional mandatory arguments
    parser.add_argument("source", help="source name", type=str, default="")

    # Positional mandatory arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')

    # Parse arguments
    args = parser.parse_args()

    return args
  
def getLocalMax(outputfile, outputFilePath, source_velocities, index_range_for_local_maxima):
    data = np.fromfile(outputFilePath + outputfile, dtype="float64", count=-1, sep=" ") .reshape((file_len(outputFilePath + outputfile),4))
    x = data[:, [0]]
    y_u1 = data[:, [1]] 
    y_u9 = data[:, [2]]
    y_avg = data[:, [2]]
    
    indexies_for_source_velocities = [0] * len(source_velocities)
    for index in range (0, len(source_velocities)):
            indexies_for_source_velocities[index] =  (np.abs(x-float(source_velocities[index]))).argmin()
                      
    max_amplitude_list_u1 = list()
    max_amplitude_list_u9 = list()
    max_amplitude_list_uavg = list()

    for index in indexies_for_source_velocities:
            max_amplitude_list_tmp_u1 = list()
            max_amplitude_list_tmp_u9 = list()
            max_amplitude_list_tmp_uavg = list()
           
            for i in range (index - index_range_for_local_maxima, index + index_range_for_local_maxima):
                max_amplitude_list_tmp_u1.append(y_u1[i])
                max_amplitude_list_tmp_u9.append(y_u9[i])
                max_amplitude_list_tmp_uavg.append(y_avg[i])
                  
            max_amplitude_list_u1.append(max_amplitude_list_tmp_u1)
            max_amplitude_list_u9.append(max_amplitude_list_tmp_u9)
            max_amplitude_list_uavg.append(max_amplitude_list_tmp_uavg)
        
    max_apmlitudes_u1 = [np.max(value) for value  in max_amplitude_list_u1]
    max_apmlitudes_u9 = [np.max(value) for value  in max_amplitude_list_u9]
    max_apmlitudes_uavg = [np.max(value) for value  in max_amplitude_list_uavg]
        
    for max in range(0,len(max_apmlitudes_u1)):
        max_apmlitudes_u1[max] = [source_velocities[max], max_apmlitudes_u1[max]]
        max_apmlitudes_u9[max] = [source_velocities[max], max_apmlitudes_u9[max]]
        max_apmlitudes_uavg[max] = [source_velocities[max], max_apmlitudes_uavg[max]]
   
    return (max_apmlitudes_u1, max_apmlitudes_u9, max_apmlitudes_uavg)


def chanegResultAmplitudes(outputfiles, resutlFile, outputFilePath, source_velocities, index_range_for_local_maxima): 
  
    with open(resutlFile) as result:
        result_json = json.load(result)
        
    for outputfile  in outputfiles:
        itarration = getIerattion(outputfile)
        max_apmlitudes_u1, max_apmlitudes_u9, max_apmlitudes_uavg = getLocalMax(outputfile, outputFilePath, source_velocities, index_range_for_local_maxima)
        
        for key in result_json.keys():
            if key.endswith("_" + itarration):
                result_json[key]["polarizationU1"] =  max_apmlitudes_u1
                result_json[key]["polarizationU9"] = max_apmlitudes_u9
                result_json[key]["polarizationAVG"] = max_apmlitudes_uavg
            
    resultFile = open (resutlFile, "w")
    resultFile.write(json.dumps(result_json, indent=2))
    resultFile.close()

def getIerattion(outputFile):
    itarration = outputFile.split(".")[0].split("_")[-1]
    return itarration
    
def getAllOutputfiles(outputFilePath, source):
    outputFiles = list()
    
    for file in os.listdir(outputFilePath):
        if file.startswith(source):
            outputFiles.append(file)
            
    return outputFiles

def main():
    args = parseArguments()
    source= str(args.__dict__["source"])
    configFilePath = str(args.__dict__["config"])
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    outputFilePath =  config.get('paths', "outputFilePath")
    resultFilePath =  config.get('paths', "resultFilePath")
    source_velocities = config.get('velocities', source).replace(" ", "").split(",")
    index_range_for_local_maxima = int(config.get('parameters', "index_range_for_local_maxima"))
    resutlFile = resultFilePath + source + ".json"
    outputfiles = getAllOutputfiles(outputFilePath, source)
    chanegResultAmplitudes(outputfiles, resutlFile, outputFilePath, source_velocities, index_range_for_local_maxima)
         
if __name__=="__main__":
    main()
