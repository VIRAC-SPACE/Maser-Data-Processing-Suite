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
    
    for file in os.listdir(outputFilePath + source):
        if file.startswith(source):
            outputFiles.append(file)
            
    return outputFiles

def getAllDatafiles(dataFilePath, source):
    dataFiles = list()

    for file in os.listdir(dataFilePath):
        if file.startswith(source + "_"):
            dataFiles.append(file)
            
    return dataFiles

def getIteration(file):
    return file.split(".")[0].split("_")[-1]

def findOutputFile(dataFileIteration, outputFileList):
    for file in outputFileList:
        if file.endswith("_" +dataFileIteration + ".dat"):
            return file
        
def correctDataRead(input):
    tmpList = list()
    
    for i in input:
        tmpList.append(i[0]) 
        
    return np.array(tmpList)
    
def createData(datafile, outputfile, cuts, outputFilePath, source):
    output_Data = np.fromfile(outputfile, dtype="float64", count=-1, sep=" ").reshape((file_len(outputfile),4))    
    output_x_Data = output_Data[:, [0]]

    output_x_first = output_x_Data[0][0]
    output_x_last = output_x_Data[-1][0]

    result = pickle.load(open(datafile, "rb"))
    x = result.getMatrix()[:, [0]]
    x = correctDataRead(x)
    x = np.flip(x,0)

    y_u1 = result.getMatrix()[:, [1]]
    y_u1 = correctDataRead(y_u1)
    y_u1 = np.flip(y_u1,0)
    
    y_u9 = result.getMatrix()[:, [2]]
    y_u9 = correctDataRead(y_u9)
    y_u9 = np.flip(y_u9,0)

    indexFirst = findNearestIndex(x, output_x_first)
    indexLast = findNearestIndex(x, output_x_last)

    x = x[indexFirst:indexLast]
    y_u1 = y_u1[indexFirst:indexLast]
    y_u9 = y_u9[indexFirst:indexLast]
    
    cutsIndex = list()
    cutsIndex.append(0)

    for cut in cuts:
        cutsIndex.append((np.abs(x-float(cut[0]))).argmin()) 
        cutsIndex.append((np.abs(x-float(cut[1]))).argmin())
            
    cutsIndex.append(len(x))

    polyArray_x = list()
    polyArray_u1 = list()
    polyArray_u9 = list()
        
    i = 0
    j = 1
        
    while i != len(cutsIndex):
        polyArray_x.append(x[cutsIndex[i] : cutsIndex[j]])
        polyArray_u1.append(y_u1[cutsIndex[i] : cutsIndex[j]])
        polyArray_u9.append(y_u9[cutsIndex[i] : cutsIndex[j]])
        i = i + 2 
        j = j + 2
            
    poly_x = list()
    poly_u1 = list()
    poly_u9 = list()
        
    for p in polyArray_x:
        for p1 in p:
            poly_x.append(p1)
                
    for p in polyArray_u1:
        for p1 in p:
            poly_u1.append(p1)
        
    for p in polyArray_u9:
        for p1 in p:
            poly_u9.append(p1)
            
    polyx = np.array(poly_x)
    polyu1 = np.array(poly_u1)
    polyu9 = np.array(poly_u9)
    
    polynomialOrder = 3    
    z_u1 = np.polyfit(polyx, polyu1, polynomialOrder)
    p_u1 = np.poly1d(z_u1)
        
    z_u9 = np.polyfit(polyx, polyu9, polynomialOrder)
    p_u9 = np.poly1d(z_u9)
    
    localMax_Array_u1 = y_u1 - p_u1(x)
    localMax_Array_u9 = y_u9 - p_u9(x)
    
    z1_NotSmoohtData = localMax_Array_u1
    z2_NotSmoohtData = localMax_Array_u9
    avg_y_NotSmoohtData = (z1_NotSmoohtData + z2_NotSmoohtData) / 2
    
    location = datafile.split("/")[-1].split(".")[0].split("_")[-2]
    date = "_".join([datafile.split("/")[-1].split(".")[0].split("_")[1], datafile.split("/")[-1].split(".")[0].split("_")[2], datafile.split("/")[-1].split(".")[0].split("_")[3]])
    time = datafile.split("/")[-1].split(".")[0].split("_")[-3]
    iteration_number = datafile.split("/")[-1].split(".")[0].split("_")[-1]
    totalResults = [x,  z1_NotSmoohtData,  z2_NotSmoohtData,  avg_y_NotSmoohtData]
    output_file_name = outputFilePath + "/NotSmooht/"  + source + "_" + time.replace(":", "_") + "_" + date.replace(" ", "_") + "_" + location + "_" + str(iteration_number) + ".dat"
    output_file_name = output_file_name.replace(" ", "")
    np.savetxt(output_file_name, np.transpose(totalResults))
           
def main():
    args = parseArguments()
    source= str(args.__dict__["source"])
    configFilePath = str(args.__dict__["config"])
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    cuts = config.getConfig('cuts', source).split(";")
    cuts = [c.split(",") for c in cuts]

    outputFilePath = config.getConfig('paths', "outputFilePath")
    dataFilePathSDR = config.getConfig('paths', "dataFilePath") + "SDR/"
    dataFilePathDBBC = config.getConfig('paths', "dataFilePath") + "DBBC/"

    outPutFiles = getAllOutputfiles(outputFilePath, source)
    dataFilesSDR = getAllDatafiles(dataFilePathSDR, source)
    dataFilesDBBC = getAllDatafiles(dataFilePathDBBC, source)

    for dataFile in dataFilesSDR:
        outputFile = findOutputFile(getIteration(dataFile), outPutFiles)

        if "NoneType" not in str(type(outputFile)):
            createData(dataFilePathSDR + dataFile, outputFilePath + source + "/" + outputFile, cuts, outputFilePath, source)

    for dataFile in dataFilesDBBC:
        outputFile = findOutputFile(getIteration(dataFile), outPutFiles)
        if "NoneType" not in str(type(outputFile)):
            createData(dataFilePathDBBC + dataFile, outputFilePath + source + "/" + outputFile, cuts, outputFilePath, source)
    
    sys.exit(0)
    
if __name__=="__main__":
    main()
    sys.exit(0)
