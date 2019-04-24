#! /usr/bin/python3
# -*- coding: utf-8 -*-

import argparse
import numpy as np
from numpy import trapz
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter
from functools import reduce

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

if __name__ == "__main__":
    gaussLines = getConfigs("gauss_lines", getArgs("source")).replace(" ", "").split(",")
    dataFile = getConfigs("paths", "notSmoohtFilePath") + getArgs("source") + "/" + getArgs("output")

    data = np.fromfile(dataFile, dtype="float64", count=-1, sep=" ").reshape((file_len(dataFile), 4))
    velocity = correctNumpyReadData(data[:, [0]])
    ampvid = correctNumpyReadData(data[:, [3]])
    indexies = [(np.abs(velocity - float(line))).argmin() for line in gaussLines]
    mons = [max(ampvid[index-5:index+5]) for index in indexies]
    gaussian = [models.Gaussian1D(mons[index], gaussLines[index], 0.05, bounds={'stddev': (None, 0.15)}) for index in range(0, len(mons))]

    def sum(a,b):
        return a + b
    gg_init = reduce(sum, gaussian)

    fitter = fitting.SLSQPLSQFitter()
    fit = LevMarLSQFitter()
    gg_fit = fit(gg_init, velocity, ampvid)
    print(gg_fit.parameters)
    sts = [models.Gaussian1D(gg_fit[index].amplitude, gg_fit[index].mean, gg_fit[index].stddev) for index in range(0, len(gaussLines))]

    #trapz(sts[0], 5)

    plt.plot(velocity, ampvid, 'b+', label="data")
    plt.plot(velocity, gg_fit(velocity), "-", label="total fit", linewidth=4,)

    colors = []

    for index in range(0, len(gaussLines)):
        if index %2:
            index -= 1
        else:
            index += 1
        colors.append("C" + str(index))

    c = 0
    for st in sts:
        #print(st(velocity)[st(velocity)>0], st(velocity)[st(velocity)>0].ndim)
        #print(trapz(st(velocity)[st(velocity)>0], 5))
        print("Area of st" + str(c), trapz(st(velocity), velocity))
        plt.plot(velocity, st(velocity), colors[c], label="ST" + str(c))
        c += 1

    plt.plot(velocity, (gg_fit(velocity) - ampvid) * 1 - 35)
    plt.plot(velocity, velocity * 0 - 35)
    plt.legend()
    plt.show()

    print("Fit value",  np.mean(ampvid - gg_fit(velocity)))
