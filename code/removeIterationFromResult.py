#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import json

from parsers._configparser import ConfigParser

def parseArguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("iteration", help="Iteration", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()

    return args

def getResult(resultFilePath, source):
    with open(resultFilePath + source + ".json") as result_data:    
            result = json.load(result_data)
            
    return result

def deleteIteration(result, iteration):
    tmpKey = None
    for key in result.keys():
        if key.endswith("_" + str(iteration)):
            tmpKey = key
        
    if tmpKey:
        del result[tmpKey]
              
    return result

def removeIterationFromResultJson(resultFilePath, source, updateResult):
    os.system("cp  " + resultFilePath + source + ".json " + resultFilePath + "/oldresult")
    resultFile = open (resultFilePath + source + ".json", "w")
    resultFile.write(json.dumps(updateResult, indent=2))
    resultFile.close()
    
def main():
    # Parse the arguments
    args = parseArguments()
    source = str(args.__dict__["source"])
    iteration_number = int(args.__dict__["iteration"])
    configFilePath = str(args.__dict__["config"])
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    resultFilePath = config.getConfig('paths', "resultFilePath")
    
    result = getResult(resultFilePath, source)
    updateResult = deleteIteration(result, iteration_number)
    removeIterationFromResultJson(resultFilePath, source, updateResult)
       
    sys.exit(0)

if __name__ == "__main__":
    main()