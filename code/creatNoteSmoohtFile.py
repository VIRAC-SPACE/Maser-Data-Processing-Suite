#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
import pickle

from parsers._configparser import ConfigParser
from help import *

class Result():
    def __init__(self, matrix, specie):
        self.matrix = matrix
        self.specie = specie
        
    def getMatrix(self):
        return self.matrix
    
    def getSpecie(self):
        return self.specie

def parseArguments():
    parser = argparse.ArgumentParser(description='''Fix Amplitudes in output file. ''', epilog="""Fix.""")
    parser.add_argument("source", help="source name", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')
    args = parser.parse_args()
    return args

def getAllOutputfiles(outputFilePath, source):
    outputFiles = list()
    
    for file in os.listdir(outputFilePath):
        if file.startswith(source):
            outputFiles.append(file)
            
    return outputFiles

def getAllDatafiles(dataFilePath, source):
    dataFiles = list()
    
    for file in os.listdir(dataFilePath):
        if file.startswith(source + "_"):
            dataFiles.append(file)
            
    return dataFiles

def createData(datafile, outputfile):
    output_Data = np.fromfile(outputfile, dtype="float64", count=-1, sep=" ").reshape((file_len(outputfile),4))
    output_x_Data = output_Data[:, [0]]
    output_x_first = output_x_Data[0][0]
    output_x_last = output_x_Data[-1][0]
    
    result = pickle.load(open(datafile, "rb"))
    x = result.getMatrix()[:, [0]]
    print(x)
    indexFirst = np.where(x==output_x_first)
    indexLast = np.where(x==output_x_last)
    x = x[0] #[indexFirst:indexLast]
    
    
    print(indexFirst, indexLast)
    
    #indexies = 

def main():
    args = parseArguments()
    source= str(args.__dict__["source"])
    configFilePath = str(args.__dict__["config"])
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    outputFilePath= config.getConfig('paths', "outputFilePath")
    dataFilePath = config.getConfig('paths', "dataFilePath")
    
    outPutFiles = getAllOutputfiles(outputFilePath, source)
    dataFiles = getAllDatafiles(dataFilePath, source)
    
    createData(dataFilePath + dataFiles[0], outputFilePath + outPutFiles[0])
    
    sys.exit(0)
    
if __name__=="__main__":
    main()