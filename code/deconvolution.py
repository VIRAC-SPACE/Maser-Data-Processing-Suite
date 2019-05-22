#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from scipy.signal import deconvolve
from astropy.convolution import Gaussian1DKernel
from functools import reduce

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

            ampTets = np.zeros(len(f))
            gaussTest = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
            b = np.fft.fft(gaussTest)
            i = 0
            j = 19

            ampConvolvList = list()
            while j<=len(f):
                ampTetsTmp = np.abs(np.fft.ifft(np.fft.fft(amp[i:j])/b))
                ampTetsTmp = np.pad(ampTetsTmp, (i, len(f) - j), 'constant', constant_values=(0, 0))
                ampConvolvList.append(ampTetsTmp)
                i += 1
                j += 1

            ampTets = reduce(lambda a,b : a+b,ampConvolvList)
            ampTets = ampTets / len(ampConvolvList)

            gauss = np.exp(-((np.linspace(0, 50) - 25.) / float(12)) ** 2)
            orginalAmp1 = np.nan_to_num(deconvolve(amp1, gauss)[0].astype('float64'))
            orginalAmp9 = np.nan_to_num(deconvolve(amp9, gauss)[0].astype('float64'))
            orginalAmp = np.nan_to_num(deconvolve(amp, gauss)[0].astype('float64'))

            n = len(f) - len(gauss) + 1
            s = int((len(f) - n) / 2)
            deconv_res = np.zeros(len(f))
            deconv_res[s:len(f) - s - 1] = orginalAmp
            orginalAmp = deconv_res

            pad = (0, f.size - orginalAmp1.size)
            results = [f, np.pad( orginalAmp1, pad, 'constant'), np.pad(orginalAmp9, pad, 'constant'), ampTets]
            np.savetxt(notSmoothFile, np.transpose(results))

    sys.exit(0)

if __name__ == "__main__":
    main()