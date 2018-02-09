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
    
    for logFileName in os.listdir(logFileDir):
        scan_numbers = ExperimentLogReader(logFileDir + logFileName, prettyLogDir).getScansForSource(source_name)
        
        for scan in scan_numbers:
            dataFile = logFileName.split(".")[0][:-2] + "_n" + scan + ".dat"
            
            print "Log file is", logFileDir + logFileName, "Data file is", dataFile
            os.system("python2  " +  "code/plotAmplitudeFrequencies.py " + logFileName + " " + dataFileDir + dataFile)
        
        
       
                
        
    sys.exit(0)