import json

key = "specie"
value = "ch3oh"

inputFile = "/home/janis/Documents/workspace-sts/DataProcessingForMaserObservation/results/" + "cepa.json"


with open(inputFile) as data_file:
     data = json.load(data_file)
     

for k in data:
    data[k][key] = value

data_file = open(inputFile, "w")
data_file.write(json.dumps(data, indent=2))