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
import threading
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

def addAreaToResultFileProcess(source, iteration, areas):
    resultDir = getConfigs('paths', 'resultFilePath')
    with open(resultDir + source + ".json") as resultData:
        results = json.load(resultData)
         
    for experiment in results:
        if experiment.endswith("_" + str(iteration)):
            results[experiment]["areas"] = areas
                 
    resultFile = open (resultDir +  source + ".json", "w")
    resultFile.write(json.dumps(results, indent=2))
    resultFile.close()

def addAreaToResultFile(source, iteration, areas):
    p = Process(target=addAreaToResultFileProcess, args=(source, iteration, areas))
    p.start()
    p.join()

class DataFileProcessing(threading.Thread):
    def __init__(self, source):
        __slots__ = ['source', 'outputDir']
        super().__init__()
        threading.Thread.__init__(self) 
        self.source = source
        self.outputDir = getConfigs('paths', 'outputFilePath')
        self.velocitiesList = list()
        self.amplitudeList = list()
        self.itarationList = list()
        
    def getOputFiles(self):
        outputFileList = [outputFile for  outputFile in   os.listdir(self.outputDir) if outputFile.startswith(self.source)]
        return outputFileList
    
    def getDataForSingleDataFile(self, outputFile):
        file = self.outputDir + outputFile
        data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),4))
        velocity = data[:, [0]]
        ampvid = data[:, [3]]
        return (velocity, ampvid)
    
    def getDataForAllDataFile(self):
        for outputFile in self.getOputFiles():
            velocity, ampvid = self.getDataForSingleDataFile(outputFile)
            self.velocitiesList.append(velocity)
            self.amplitudeList.append(ampvid)
            iteration = outputFile.split(".")[0].split("_")[-1]
            self.itarationList.append(iteration)
            
    def getData(self):
        return (self.velocitiesList, self.amplitudeList, self.itarationList)
        
    def run(self):
        self.getDataForAllDataFile()

class VelocityDensity(threading.Thread):
    def __init__(self, velocity, amplitude, source, itaration):
        __slots__ = ['velocity', 'amplitude', 'source_velocities', 'source', 'itaration', 'areasList']
        super().__init__()
        threading.Thread.__init__(self) 
        self.velocity = np.array([v[0] for v in velocity])
        self.amplitude = np.array([a[0] for a in amplitude])
        self.source = source
        self.source_velocities = getConfigs('velocities', self.source).split(",")
        self.itaration = itaration
        self.areas = list()
        
    def gauss(self, A, μ, σ):
        return float(A) / (float(σ) * np.sqrt(2 * np.pi)) * np.exp(-(self.velocity-float(μ))**2 / (2*float(σ)**2))
    
    def createGaussFit(self, vel, A, sigma):
        return  self.gauss(A, float(vel), sigma)
    
    def cost(self,parameters):
        a, b, c = parameters
        return np.sum(np.power(self.gauss(a, b, c) - self.amplitude, 2)) / self.velocity.size
    
    def computeDensity(self):
        dx = np.abs(np.abs(self.velocity[-1]) - np.abs(self.velocity[0]))/self.velocity.size
        
        i = 0        
        for vel in self.source_velocities:
            index = (np.abs(self.velocity - float(vel))).argmin()
            initial_guess = [float(vel), self.amplitude[index], 0.38]
            result = optimize.minimize(self.cost, initial_guess)
            g = self.createGaussFit(result.x[0], result.x[1], result.x[2])
            area = trapz(g, dx=dx)
            self.areas.append(area)
            i = i + 1
                    
    def write(self):
        addAreaToResultFile(self.source, self.itaration, self.areas)
        
    def run(self):
        self.computeDensity()
                
class VelocityDensityPolter(VelocityDensity):
    def __init__(self, velocity, amplitude, source):
        __slots__ = ['velocity', 'amplitude', 'source_velocities', 'source']
        self.velocity = np.array([v[0] for v in velocity])
        self.amplitude = np.array([a[0] for a in amplitude])
        self.source = source
        self.source_velocities = getConfigs('velocities', self.source).split(",")
          
    def plotDensity(self):
        colors = ['b', 'g', 'c', 'm', 'y', 'k', 'w']
        fig = plt.figure("Gauss fit")
        ax = fig.add_subplot(111)
        ax.set_xlabel('Velocity (km sec$^{-1}$)')
        ax.set_ylabel('Flux density (Jy)')
        plt.plot(self.velocity, self.amplitude, "r", label="spectre")
        i = 0
        
        for vel in self.source_velocities:
            index = (np.abs(self.velocity - float(vel))).argmin()
            initial_guess = [float(vel), self.amplitude[index], 0.38]
            result = optimize.minimize(self.cost, initial_guess)
            g = self.createGaussFit(result.x[0], result.x[1], result.x[2])
            plt.plot(self.velocity, g, "--" + colors[i], label=float(vel))
            i = i + 1
            
        plt.grid()    
        plt.legend()
        plt.show()
        
def main():
    args = parseArguments()
    output = str(args.__dict__["output"])
    source = str(args.__dict__["source"])
    
    dataFileProcessing = DataFileProcessing(source)
    
    if output == "":
        dataFileProcessing.start()
        dataFileProcessing.join()
        velocitiesList, amplitudeList, itarationList = dataFileProcessing.getData()
        
        for i in range(0, len(velocitiesList)):
            velocityDensity = VelocityDensity(velocitiesList[i], amplitudeList[i], source, itarationList[i])
            velocityDensity.start()
            velocityDensity.join()
            velocityDensity.write()
                    
    else:
        velocity, amplitude = dataFileProcessing.getDataForSingleDataFile(output)
        velocityDensity = VelocityDensityPolter(velocity, amplitude, source)
        velocityDensity.plotDensity()
         
    sys.exit(0)
    
if __name__ == "__main__":
    p = Process(target=main)
    p.start()
    #p.join()
    