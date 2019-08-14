import json
import os

key = "type"
value = "DBBC"
inputDir = "results/"

for resultFile in os.listdir(inputDir):
    with open(inputDir + resultFile, "r") as data_file:
        data = json.load(data_file)

    for k in data:
        data[k][key] = value

    data_file = open(inputDir + resultFile, "w")
    data_file.write(json.dumps(data, indent=2))
    data_file.close()