#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from scipy.integrate import trapz
import matplotlib.pyplot as plt
import json
from parsers._configparser import ConfigParser
from help import *

def parseArguments():
    parser = argparse.ArgumentParser(description='''Velocity Density. ''', epilog="""Velocity Density.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')
    parser.add_argument("source", help="Experiment source", type=str, default="")
    args = parser.parse_args()
    return args

def getArgs(key):
    args = parseArguments()
    return str(args.__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

def getOputFiles(source):
    outputDir = getConfigs('paths', 'outputFilePath')
    outputFileList = list()
    for outputFile in os.listdir(outputDir):
        if outputFile.startswith(source):
            outputFileList.append(outputFile)
    return outputFileList

def getData(outputFiles):
    outputDir = getConfigs('paths', 'outputFilePath')
    velocitiesList = list()
    amplitudeList = list()
    itarationList = list()
    for outputFile in outputFiles:
        file = outputDir + outputFile
        data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),4))
        velocity = data[:, [0]]
        ampvid = data[:, [3]]
        velocitiesList.append(velocity)
        amplitudeList.append(ampvid)
        iteration = outputFile.split(".")[0].split("_")[-1]
        itarationList.append(iteration)
    return (velocitiesList, amplitudeList, itarationList)

def addAreaToResultFile(source, iteration, areas):
    resultDir = getConfigs('paths', 'resultFilePath')
    with open(resultDir + source + ".json") as resultData:
        results = json.load(resultData)
    
    for experiment in results:
        if experiment.endswith("_" + str(iteration)):
            results[experiment]["areas"] = areas
            
    resultFile = open (resultDir +  source + ".json", "w")
    resultFile.write(json.dumps(results, indent=2))
    resultFile.close()
    
def gauss(x, A, μ, σ):
    return A / (σ * np.sqrt(2 * np.pi)) * np.exp(-(x-μ)**2 / (2*σ**2))
  
def computeDensity(velocity, ampvid, source_velocities, source, iteration):
    velocity = np.array([v[0] for v in velocity])
    ampvid = np.array([a[0] for a in ampvid])
    dx = np.abs(np.abs(velocity[-1]) - np.abs(velocity[0]))/velocity.size
    areas = list()
    #colors = ['b', 'g', 'c', 'm', 'y', 'k', 'w']
    #fig = plt.figure()
    #plt.plot(velocity, ampvid, "r")
    i = 0
    for vel in source_velocities:
        index = (np.abs(velocity - float(vel))).argmin()
        g = gauss(velocity, ampvid[index], float(vel), 0.1)
        #plt.plot(velocity, g, "--" + colors[i], label="--" + colors[i])
        area = trapz(g, dx=dx)
        areas.append(area)
        addAreaToResultFile(source, iteration, areas)
        i = i + 1
        
    #plt.legend()
    #plt.show()
    #print ("area", areas)

def main():
    args = parseArguments()
    source = str(args.__dict__["source"])
    velocitiesList, amplitudeList, itarationList = getData(getOputFiles(source))
    source_velocities = getConfigs('velocities', source).split(",")
    
    for i in range(0, len(velocitiesList)):
        computeDensity(velocitiesList[i], amplitudeList[i], source_velocities, source, itarationList[i])
    sys.exit(0)
    
if __name__ == "__main__":
    main()