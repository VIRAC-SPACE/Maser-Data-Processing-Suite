#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from scipy.signal import deconvolve
from astropy.convolution import Gaussian1DKernel

from parsers._configparser import ConfigParser
from help import *

def parseArguments():
    parser = argparse.ArgumentParser(description='''Creates input file for plotting too. ''', epilog="""Deconvolution.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    return args

def getArgs(key):
    return str(parseArguments().__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

def getIteration(fileName):
    return fileName.split(".")[0].split("_")[-1]

def findOutputFile(dataFileIteration, outputFileList):
    for file in outputFileList:
        if file.endswith("_" +dataFileIteration + ".dat"):
            return file
def main():
    smoothPath = getConfigs("paths", "outputFilePath") + getArgs("source")
    notSmoothPath = getConfigs("paths", "notSmoohtFilePath") + getArgs("source")

    smoothFiles = [file for file in os.listdir(smoothPath) if file.startswith(getArgs("source") + "_") ]
    smoothIterations = [getIteration(file) for file in os.listdir(smoothPath) if file.startswith(getArgs("source") + "_") ]
    notSmoothIterations = [getIteration(file) for file in os.listdir(notSmoothPath)]

    for iteration in smoothIterations:
        if iteration not in notSmoothIterations:
            smoothFile = smoothPath +"/"+ findOutputFile(iteration, smoothFiles)
            notSmoothFile = notSmoothPath +"/"+ findOutputFile(iteration, smoothFiles) + ".new"

            print("Deconvolving", smoothFile)
            data = np.fromfile(smoothFile, dtype="float64", count=-1, sep=" ").reshape((file_len(smoothFile), 4))
            f = correctNumpyReadData(data[:, [0]])
            amp1 = correctNumpyReadData(data[:, [1]])
            amp9 = correctNumpyReadData(data[:, [2]])
            amp = correctNumpyReadData(data[:, [3]])

            gauss = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
            orginalAmp1 = np.nan_to_num(deconvolve(amp1, gauss)[0].astype('float64'))
            orginalAmp9 = np.nan_to_num(deconvolve(amp9, gauss)[0].astype('float64'))
            orginalAmp = np.nan_to_num(deconvolve(amp, gauss)[0].astype('float64'))

            pad = (0, f.size - orginalAmp1.size)
            results = [f, np.pad( orginalAmp1, pad, 'constant'), np.pad(orginalAmp9, pad, 'constant'), np.pad(orginalAmp, pad, 'constant')]
            np.savetxt(notSmoothFile, np.transpose(results))

    sys.exit(0)

if __name__ == "__main__":
    main()