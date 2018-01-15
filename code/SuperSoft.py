import os
import sys
from subprocess import call

from experimentsLogReader import ExperimentLogReader

def usage():
    print ('Usage: ' + sys.argv[0] + ' Source name')

if __name__=="__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    
    source_name = sys.argv[1]
    
    logFileDir = "logs/"
    dataFileDir = "dataFiles/"
    prettyLogDir = "prettyLogs/"
      
    logFileList = list()
    dataFileList = list()
    
    for filename in os.listdir(logFileDir):
        logs  = ExperimentLogReader(logFileDir + filename).getLgs() 
        scanNumbers = ExperimentLogReader(logFileDir + filename).getAllScansNumbers()
        
        for i in range(0, len(scanNumbers)):
            scan = logs[scanNumbers[i]]
            sourceNameInLogs = scan["sourceName"]
            
            if sourceNameInLogs ==  source_name:
                datafile = dataFileDir + filename.split(".")[0][:-2] + "_n" + scanNumbers[i] + ".dat"
                print logFileDir + filename, datafile
                os.system("python2  " +  "code/plotAmplitudeFrequencies.py " + logFileDir + filename + " " + datafile)
                #call(["python2 ", "code/plotAmplitudeFrequencies.py" ,logFileDir + filename, datafile])
        
        
        #sourceFromLogFile = logs[]
    sys.exit(0)