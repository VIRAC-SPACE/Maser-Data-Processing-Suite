#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os,sys
import argparse
import numpy as np
from numpy import trapz
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter
from functools import reduce
import json

from parsers._configparser import ConfigParser
from help import *

def parseArguments():
    parser = argparse.ArgumentParser(description='''Velocity Density. ''', epilog="""Velocity Density.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')
    parser.add_argument("-o", "--output", help="Plot fits for outputfile ", default="")
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

def iterationFromDataFile(dataFileName):
    return dataFileName.split(".")[0].split("_")[-1]

def getLocalMaximum(xarray, yarray, source):
    dataPoints = len(xarray)
    cuts = getConfigs("cuts", source).split(";")
    cuts = [c.split(",") for c in cuts]
    cutsIndex = list()
    cutsIndex.append(0)

    for cut in cuts:
        cutsIndex.append((np.abs(xarray - float(cut[0]))).argmin())
        cutsIndex.append((np.abs(xarray - float(cut[1]))).argmin())

    cutsIndex.append(dataPoints)

    polyArray_x = list()
    polyArray_y = list()

    while i != len(cutsIndex):
        polyArray_x.append(xarray[cutsIndex[i]: cutsIndex[j]])
        polyArray_y.append(yarray[cutsIndex[i]: cutsIndex[j]])
        i = i + 2
        j = j + 2

    poly_x = list()
    poly_y = list()

    for p in polyArray_x:
        for p1 in p:
            poly_x.append(p1)

    for p in polyArray_y:
        for p1 in p:
            poly_y.append(p1)

    polyx = np.array(poly_x)
    polyy = np.array(poly_y)

    z = np.polyfit(polyx, polyy, 3)
    p = np.poly1d(z)
    localMax_Array = yarray - p(xarray)

def computeGauss(file):
    gaussLines = getConfigs("gauss_lines", getArgs("source")).replace(" ", "").split(",")
    dataFile = getConfigs("paths", "notSmoohtFilePath") + getArgs("source") + "/" + file  #getArgs("output")

    data = np.fromfile(dataFile, dtype="float64", count=-1, sep=" ").reshape((file_len(dataFile), 4))
    velocity = correctNumpyReadData(data[:, [0]])
    ampvid = correctNumpyReadData(data[:, [3]])
    indexies = [(np.abs(velocity - float(line))).argmin() for line in gaussLines]
    mons = [max(ampvid[index - 5:index + 5]) for index in indexies]
    gaussian = [models.Gaussian1D(mons[index], gaussLines[index], 0.05, bounds={'stddev': (None, 0.15)}) for index in range(0, len(mons))]

    def sum(a, b):
        return a + b

    gg_init = reduce(sum, gaussian)
    fitter = fitting.SLSQPLSQFitter()
    fit = LevMarLSQFitter()
    gg_fit = fit(gg_init, velocity, ampvid)
    sts = [models.Gaussian1D(gg_fit[index].amplitude, gg_fit[index].mean, gg_fit[index].stddev) for index in range(0, len(gaussLines))]
    gaussianAreas = []

    for st in sts:
        gaussianAreas.append(trapz(st(velocity), velocity))

    return (gaussianAreas, sts, gg_fit, velocity, ampvid, gaussLines)

def addAreasToResultFiles(file, gaussianAreas):
    with open(getConfigs("paths", "resultFilePath") + getArgs("source") + ".json", "r") as resultFile:
        results = json.load(resultFile)

    for experiment in results:
        if int(results[experiment]["Iteration_number"]) == int(iterationFromDataFile(file)):
            results[experiment]["areas"] = gaussianAreas

    with open(getConfigs("paths", "resultFilePath") + getArgs("source") + ".json", "w") as resultFile:
        resultFile.write(json.dumps(results, indent=2))

if __name__ == "__main__":
    if getArgs("output") != "":
        gaussianAreas, sts, gg_fit, velocity, ampvid, gaussLines = computeGauss(getArgs("output"))

        plt.figure("Gausian fits")
        plt.plot(velocity, ampvid, 'C0+', label="data")
        plt.plot(velocity, gg_fit(velocity), "C1-", label="total fit", linewidth=4)

        colors = []

        for index in range(2, len(gaussLines) + 2):
            if index %2:
                index -= 1
            else:
                index += 1
            colors.append("C" + str(index))

        c = 0
        for st in sts:
            plt.plot(velocity, st(velocity), colors[c], label="ST" + str(c))
            c += 1

        plt.plot(velocity, (gg_fit(velocity) - ampvid) * 1 - 35)
        plt.plot(velocity, velocity * 0 - 35)
        plt.legend()
        plt.show()

        print("Fit value",  np.mean(ampvid - gg_fit(velocity)))

        addAreasToResultFiles(getArgs("output"), gaussianAreas)

    else:
        for resultFile in os.listdir(getConfigs("paths", "notSmoohtFilePath") + getArgs("source")):
            print("processing file", resultFile)
            addAreasToResultFiles(resultFile, computeGauss(resultFile)[0])

    sys.exit(0)
