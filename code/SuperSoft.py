import os
import sys

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
      
    logFileList = list()
    dataFileList = list()
    
    for filename in os.listdir(logFileDir):
        logs  = ExperimentLogReader(logFileDir + filename, prettyLogDir).getLgs() 
        scanNumbers = ExperimentLogReader(logFileDir + filename, prettyLogDir).getAllScansNumbers()
        
        for i in range(0, len(scanNumbers)):
            scan = logs[scanNumbers[i]]
            sourceNameInLogs = scan["sourceName"]
            
            if sourceNameInLogs ==  source_name:
                datafile = dataFileDir + filename.split(".")[0][:-2] + "_n" + scanNumbers[i] + ".dat"
                print "Log file is", logFileDir + filename, "Data file is", datafile
                os.system("python2  " +  "code/plotAmplitudeFrequencies.py " + filename + " " + datafile)
        
    sys.exit(0)