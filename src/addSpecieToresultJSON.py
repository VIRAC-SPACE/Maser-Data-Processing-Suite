import json
import os

key = "specie"
value = "ch3oh"

inputDir = "/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/results/" 

for resultFile in os.listdir(inputDir):
    with open(inputDir + resultFile) as data_file:
        data = json.load(data_file)
     
    for k in data:
        data[k][key] = value
        
    data_file = open(resultFile, "w")
    data_file.write(json.dumps(data, indent=2))