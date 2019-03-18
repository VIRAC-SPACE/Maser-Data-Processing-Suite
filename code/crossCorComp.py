#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
import matplotlib.pyplot as plt

from help import *
from parsers._configparser import ConfigParser

def parseArguments():
    parser = argparse.ArgumentParser(description='''Cross - Correlatet two maser componets. ''', epilog="""CrossCorr.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("a", help="componet A ", type=int)
    parser.add_argument("b", help="componet B", type=str)
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

def getData():
    monitoringFile = "monitoring/" + getArgs("source") + ".txt"
    a = int(getArgs("a"))
    b = int(getArgs("b"))
    compunetCount = len(getConfigs("velocities", getArgs("source")).replace(" ", "").split(","))
    data = np.fromfile(monitoringFile, dtype="float64", count=-1, sep=" ").reshape((file_len(monitoringFile),compunetCount + 1))    
    x = correctNumpyReadData(data[:, [a]])
    y = correctNumpyReadData(data[:, [b]])
    z = correctNumpyReadData(data[:, [1]])
    return (x, y, z)

def main():
    x = getData()[0] / getData()[2]
    y = getData()[1] / getData()[2]
    print ("corr coef", np.corrcoef(x,y))
    crossCorr = np.correlate(x, y, "full")
    plt.plot(crossCorr)
    plt.grid(True)
    plt.show()
    sys.exit(0)
    
if __name__ == "__main__":
    main()
    