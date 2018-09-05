import json
import os
import pickle

class Result():
    def __init__(self, matrix, specie):
        self.matrix = matrix
        self.specie = specie
        
    def getMatrix(self):
        return self.matrix
    
    def getSpecie(self):
        return self.specie
    
    
def getLocalMax(data):
    pass


resultDir = "/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/results/" 
outputDir = "/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/output/" 

for resultFile in os.listdir(outputDir):
    result = pickle.load(open(resultFile, "rb"))
    specie = result.getSpecie()
    data = result.getMatrix()
    
    resultFile = resultDir + specie + ".json"
    
    localMaxima =  getLocalMax(data)