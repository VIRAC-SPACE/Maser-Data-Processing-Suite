#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
from scipy.signal import correlate, resample
import matplotlib.pyplot as plt

from help import *
from parsers._configparser import ConfigParser

def parseArguments():
    parser = argparse.ArgumentParser(description='''Cross - Correlatet two maser componets. ''', epilog="""CrossCorr.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("file", help="Experiment source", type=str, default="")
    parser.add_argument("a", help="componet A ", type=int)
    parser.add_argument("b", help="componet B", type=str)
    parser.add_argument("--index1", help="index 1", type=str, default="")
    parser.add_argument("--index2", help="index 2", type=str, default="")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args

def getArgs(key):
    return str(parseArguments().__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

def getData(file):
    a = int(getArgs("a"))
    b = int(getArgs("b"))
    compunetCount = len(getConfigs("velocities", getArgs("source")).replace(" ", "").split(","))
    data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),compunetCount + 1))
    x = correctNumpyReadData(data[:, [a]])
    y = correctNumpyReadData(data[:, [b]])
    time = correctNumpyReadData(data[:, [0]])
    return (x, y, time)

def main():
    index1 = getArgs("index1")
    index2 = getArgs("index2")
    file = getConfigs("paths", "monitoringFilePath") + "/" + getArgs("file")

    if index1 == "":
        index1 = 0

    if index2 == "":
        index2 = len(getData(file)[0])

    index1 = int(index1)
    index2 = int(index2)

    x = getData(file)[0][index1:index2]
    y = getData(file)[1][index1:index2]
    time = getData(file)[2][index1:index2]

    print("Resampling\n")

    xResample = resample(x, len(x), t=time)
    yResample = resample(y, len(y), t=time)

    plt.subplot(1, 2, 1)
    plt.plot(xResample[1], xResample[0])
    plt.xlabel('Time')
    plt.ylabel('Resampled components A' )
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(yResample[1], yResample[0])
    plt.xlabel('Time')
    plt.ylabel('Resampled components B')
    plt.grid(True)

    plt.show()

    print("corr coef befor resampling", np.corrcoef(x, y), "\n\n")
    print ("corr coef after resampling", np.corrcoef(xResample[0], yResample[0]), "\n\n")

    print("Direct")
    crossCorr = np.correlate(xResample[0], yResample[0], "full")
    a = len(xResample[0])
    b = len(crossCorr)
    ratio = a/b
    print("Input point count", a)
    print("Output point count", b)
    print("Ratio between input data length and cross-correlation is", ratio)
    points = np.linspace(0, b, b)
    time = points * ratio
    print("Delay is", time[np.argmax(crossCorr)], "time units")

    plt.subplot(1,2,1)
    plt.plot(time, crossCorr)
    plt.xlabel('Time')
    plt.ylabel('Cross-correlation between components ' + getArgs("a") + " " + getArgs("b"))
    plt.grid(True)
    plt.title("Direct")

    print("\n\nFFT")
    crossCorr = correlate(xResample[0], yResample[0], "full", "fft")
    a = len(x)
    b = len(crossCorr)
    ratio = a / b
    print("Input point count", a)
    print("Output point count", b)
    print("Ratio between input data length and cross-correlation is", ratio)
    points = np.linspace(0, b, b)
    time = points * ratio
    print("Delay is", time[np.argmax(crossCorr)], "time units")

    plt.subplot(1, 2, 2)
    plt.plot(time, crossCorr)
    plt.xlabel('Time')
    plt.ylabel('Cross-correlation between components ' + getArgs("a") + " " + getArgs("b"))
    plt.grid(True)
    plt.title("FFT")
    plt.show()

    sys.exit(0)
    
if __name__ == "__main__":
    main()
