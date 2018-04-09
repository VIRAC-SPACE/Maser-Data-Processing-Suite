#! /usr/bin/python
import os
import sys

try:
    import json
except:
    import simplejson as json
    pass

from experimentsLogReader import ExperimentLogReader

def usage():
    print ('Usage: ' + sys.argv[0] + ' Source name')

if __name__=="__main__":
    
    if len(sys.argv) < 1:
        usage()
        sys.exit(1)
    
    source_name = sys.argv[1]
     
    #paths 
    logFileDir = "logs/"
    dataFileDir = "dataFiles/"
    prettyLogDir = "prettyLogs/"
    resultDir = "results/"
    resultFileName = source_name + ".json"
    
    if os.path.isfile(resultDir + resultFileName):
        pass
    else:
        os.system("touch " + resultDir +  resultFileName)
            
        resultFile = open (resultDir +  resultFileName, "w")
        resultFile.write("{ \n" + "\n}")
        resultFile.close()
    
    #open result file    
    with open(resultDir + resultFileName) as result_data:    
        result = json.load(result_data)
    
    #check if experiment has result
    for logFileName in os.listdir(logFileDir):
        experName = logFileName.split(".")[0][:-2]
        
        if experName in result:
            print ("Experiment " + experName + " already is processed")
        
        else:
            scan_numbers = None
            try:
                scan_numbers = ExperimentLogReader("logs/" + logFileName, prettyLogDir, []).getScansForSource(source_name) # find all scansg
            except:
                print("Unexpected error:", sys.exc_info()[0])
                print ("Got Logreader Error")
                
            if scan_numbers != None:
                print (scan_numbers)
                for scan in scan_numbers:
                    dataFile = experName + "_n" + scan + ".dat"
                    
                    print ("Log file is", logFileDir + logFileName, "Data file is", dataFile)
                    os.system("python3  " +  "code/plotAmplitudeFrequencies.py " + logFileName + " " + dataFile)
                  
    sys.exit(0)