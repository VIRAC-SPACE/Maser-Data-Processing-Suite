import sys, os
import numpy as np
import matplotlib.pyplot  as plt

try:
    import json
except:
    import simplejson as json
    pass #paths 

if __name__=="__main__":
    
    if len(sys.argv) < 1:
        #usage()
        sys.exit(1)
    
    source_name = sys.argv[1]  
    logFileDir = "logs/"
    dataFileDir = "dataFiles/"
    prettyLogDir = "prettyLogs/"
    resultDir = "results/"
    resultFileName = source_name + ".json"
    
    with open(resultDir + resultFileName) as result_data:    
            result = json.load(result_data)
          
    for experiment in result:
        for scan in result[experiment]:
            
            amplitudes_for_u1 =  result[experiment][scan]["amplitude_for_polarizationU1"]
            #amplitudes_for_u9 =  result[experiment][scan]["amplitude_for_polarizationU9"]
            #amplitudes_for_avg = result[experiment][scan]["amplitude_for_polarizationAVG"]
            print amplitudes_for_u1
    
    '''        
    amplitudes = np.array([amplitudes_for_u1, amplitudes_for_u9, amplitudes_for_avg])
    
    row, col = amplitudes.shape
    figure = plt.Figure()
    for c in range (0, col):
        amplitude = amplitudes[:, [c]]
        plt.plot(amplitude, "ro")
        plt.show()
    '''      