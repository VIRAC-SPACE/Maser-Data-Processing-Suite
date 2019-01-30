#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import numpy as np
from scipy.integrate import trapz
from scipy import optimize
import matplotlib.pyplot as plt
import json
from multiprocessing import Process
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

def getData2(outputFile):
    outputDir = getConfigs('paths', 'outputFilePath')
    file = outputDir + outputFile
    data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),4))
    velocity = data[:, [0]]
    ampvid = data[:, [3]]
    return (velocity, ampvid)

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
    return float(A) / (float(σ) * np.sqrt(2 * np.pi)) * np.exp(-(x-float(μ))**2 / (2*float(σ)**2))

def createGaussFit(velocity, vel, amplitude, sigma):
    g = gauss(velocity, amplitude, float(vel), sigma)
    return g
  
def computeDensity(velocity, ampvid, source_velocities, source, iteration):
    velocity = np.array([v[0] for v in velocity])
    ampvid = np.array([a[0] for a in ampvid])
    dx = np.abs(np.abs(velocity[-1]) - np.abs(velocity[0]))/velocity.size
    areas = list()
    i = 0
    def cost(parameters):
        a, b, c = parameters
        return np.sum(np.power(gauss(velocity, a, b, c) - ampvid, 2)) / velocity.size
    for vel in source_velocities:
        index = (np.abs(velocity - float(vel))).argmin()
        initial_guess = [float(vel), ampvid[index], 0.38]
        result = optimize.minimize(cost, initial_guess)
        g = createGaussFit(velocity, result.x[0], result.x[1], result.x[2])
        area = trapz(g, dx=dx)
        areas.append(area)
        addAreaToResultFile(source, iteration, areas)
        i = i + 1
            
def plotDensity(velocity, amplitude, source_velocities):
    colors = ['b', 'g', 'c', 'm', 'y', 'k', 'w']
    fig = plt.figure("Gauss fit")
    ax = fig.add_subplot(111)
    ax.set_xlabel('Velocity (km sec$^{-1}$)')
    ax.set_ylabel('Flux density (Jy)')
    plt.plot(velocity, amplitude, "r", label="spectre")
    i = 0
    def cost(parameters):
        a, b, c = parameters
        return np.sum(np.power(gauss(velocity, a, b, c) - amplitude, 2)) / velocity.size
    for vel in source_velocities:
        index = (np.abs(velocity - float(vel))).argmin()
        initial_guess = [float(vel), amplitude[index][0], 0.38]
        result = optimize.minimize(cost, initial_guess)
        g = createGaussFit(velocity, result.x[0], result.x[1], result.x[2])
        plt.plot(velocity, g, "--" + colors[i], label=float(vel))
        i = i + 1
        
    plt.grid()    
    plt.legend()
    plt.show()

def main():
    args = parseArguments()
    output = str(args.__dict__["output"])
    source = str(args.__dict__["source"])
    source_velocities = getConfigs('velocities', source).split(",")
    
    if output == "":
        velocitiesList, amplitudeList, itarationList = getData(getOputFiles(source))
        
        for i in range(0, len(velocitiesList)):
            p = Process(target=computeDensity, args=(velocitiesList[i], amplitudeList[i], source_velocities, source, itarationList[i],))
            p.start()
            p.join()
    else:
        velocity, amplitude = getData2(output)
        plotDensity(velocity, amplitude, source_velocities)
        
    sys.exit(0)
    
if __name__ == "__main__":
    p = Process(target=main)
    p.start()
    #p.join()
    
    
    