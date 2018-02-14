import sys, os
import numpy as np
import matplotlib.pyplot  as plt
from matplotlib.dates import date2num
import pandas as pd
from datetime import datetime

try:
    import json
except:
    import simplejson as json
    pass #paths 

if __name__=="__main__":
    
    if len(sys.argv) < 1:
        #usage()
        sys.exit(1)
    
    months = {"Jan":"1", "Feb":"2", "Mar":"3", "Apr":"4", "May":"5", "Jun":"6", "Jul":"7", "Aug":"8", "Sep":"9", "Oct":"10", "Nov":"11", "Dec":"12"}
    
    source_name = sys.argv[1]  
    logFileDir = "logs/"
    dataFileDir = "dataFiles/"
    prettyLogDir = "prettyLogs/"
    resultDir = "results/"
    resultFileName = source_name + ".json"
    
    with open(resultDir + resultFileName) as result_data:    
            result = json.load(result_data)
    
    amplitude_for_u1_monitoring = dict()     
    for experiment in result:
        for scan in result[experiment]:
            scanData = result[experiment][scan]
            amplitudes_for_u1 = scanData["polarizationU1"]
            amplitude_for_u1 = amplitudes_for_u1[0][1]
            time = scanData["startTime"]
            date = scanData["Date"]
            dates = date.split(" ")
            monthsNumber = months[dates[1]]
            dates[1] = monthsNumber
            date = " ".join(dates)
            key_u1 = date + " " + time
            dateNumber = datetime.strptime(key_u1, '%d %m %Y %H:%M:%S')
            amplitude_for_u1_monitoring[dateNumber] = amplitude_for_u1
            #amplitudes_for_u9 =  result[experiment][scan]["amplitude_for_polarizationU9"]
            #amplitudes_for_avg = result[experiment][scan]["amplitude_for_polarizationAVG"]
    
    dataStrings = sorted(amplitude_for_u1_monitoring)
    values = list() 
    for data in dataStrings:
        values.append(amplitude_for_u1_monitoring[data])
    x = [date2num(date) for date in  dataStrings]
    y = values

    fig = plt.figure()

    graph = fig.add_subplot(111)

    # Plot the data as a red line with round markers
    graph.plot(x,y,'ro')

    # Set the xtick locations to correspond to just the dates you entered.
    graph.set_xticks(x)

    # Set the xtick labels to correspond to just the dates you entered.
    #graph.set_xticklabels([date.strftime("%Y-%m-%d") for (date, value) in amplitude_for_u1_monitoring])

    plt.show()