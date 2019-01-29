#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from astropy.modeling import models, fitting
from scipy.integrate import simps
from parsers._configparser import ConfigParser
from help import *

def parseArguments():
    parser = argparse.ArgumentParser(description='''Velocity Density. ''', epilog="""Velocity Density.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')
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
    for outputFile in outputFiles:
        file = outputDir + outputFile
        data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),4))
        velocity = data[:, [0]]
        ampvid = data[:, [3]]
        velocitiesList.append(velocity)
        amplitudeList.append(ampvid)
    return (velocitiesList, amplitudeList)
    
def computeDensity(velocity, ampvid, vel):
    index = (np.abs(velocity - vel)).argmin()
    g_init = models.Gaussian1D(amplitude=ampvid[index], mean=np.mean(ampvid), stddev=np.std(ampvid))
    #g_init = models.Gaussian1D(amplitude=ampvid[index], mean=1, stddev=0)
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, velocity, ampvid)
    dx = np.abs(np.abs(velocity[-1]) - np.abs(velocity[0]))/velocity.size
    area = np.sum(simps(g(ampvid), dx=dx))
    print ("area", area)

def main():
    source = "cepa"
    velocitiesList, amplitudeList = getData(getOputFiles(source))
    source_velocities = getConfigs('velocities', source).split(",")
    
    for vel in source_velocities:
        for i in range(0, len(velocitiesList)):
            computeDensity(velocitiesList[i], amplitudeList[i], float(vel))
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()