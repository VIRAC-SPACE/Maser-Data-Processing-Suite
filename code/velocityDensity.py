#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from lmfit.models import ExponentialModel, GaussianModel
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
    
def computeDensity(velocity, ampvid):
    exp_mod = ExponentialModel(prefix='exp_')
    pars = exp_mod.guess(ampvid, x=velocity)
    gauss1 = GaussianModel(prefix='g1_')
    pars.update(gauss1.make_params())
    
    pars['g1_center'].set(105, min=75, max=125)
    pars['g1_sigma'].set(15, min=3)
    pars['g1_amplitude'].set(2000, min=10)
    print (pars)

def main():
    source = "cepa"
    velocitiesList, amplitudeList = getData(getOputFiles(source))
    source_velocities = getConfigs('velocities', source).split(",")
    
    for i in range(0, len(velocitiesList)):
        computeDensity(velocitiesList[i], amplitudeList[i])
    
    sys.exit(0)
    
if __name__ == "__main__":
    main()