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
         
    with open(resultDir + resultFileName) as result_data:    
        result = json.load(result_data)
    
    for logFileName in os.listdir(logFileDir):
        if logFileName.split(".")[0][:-2] in result:
            pass
        else:
            result[logFileName.split(".")[0][:-2]] = dict()
            scan_numbers = ExperimentLogReader(logFileDir + logFileName, prettyLogDir).getScansForSource(source_name)
            
            for scan in scan_numbers:
                result[logFileName.split(".")[0][:-2]][scan] =  dict()
                dataFile = logFileName.split(".")[0][:-2] + "_n" + scan + ".dat"
                
                print "Log file is", logFileDir + logFileName, "Data file is", dataFile
                os.system("python2  " +  "code/plotAmplitudeFrequencies.py " + logFileName + " " + dataFileDir + dataFile)
                
    print(json.dumps(result, indent=4))    
    sys.exit(0)