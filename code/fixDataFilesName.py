import os

dataFileDir = "/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/dataFiles/s255/20.03.2018"
for dataFile in os.listdir(dataFileDir):
    print dataFile, dataFile[5:len(dataFile)]
    new = dataFile[5:len(dataFile)]
    os.rename(dataFileDir + "/" + dataFile, dataFileDir+"/"+new)