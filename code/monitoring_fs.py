#! /usr/bin/python
import sys
import matplotlib.pyplot  as plt
from matplotlib.dates import date2num
import mplcursors
from datetime import datetime
import json
import argparse
import configparser
from operator import itemgetter

from result import  Result

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''',
    epilog="""Monitor.""")

    # Positional mandatory arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("source", help="Source to monitor", type=str)

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')

    # Parse arguments
    args = parser.parse_args()

    return args

def main():
    # Parse the arguments
    args = parseArguments()
    configFilePath = str(args.__dict__["config"])
    source = str(args.__dict__["source"])
     
    #Creating config parametrs
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    resultDir = config.get('paths', "resultFilePath")
    source_velocities = config.get('velocities', source).split(",")
    velocities_range = 0.06
    result_list = list()
    velocitie_dict = {"u1":dict(), "u9":dict(), "avg":dict()}
    iteration_list = list()
    date_list = list()
    
    for velocitie in velocitie_dict:
        for vel in source_velocities:
            velocitie_dict[velocitie][vel] = list()
                            
    resultFileName = source + ".json"
    
    months = {"Jan":"1", "Feb":"2", "Mar":"3", "Apr":"4", "May":"5", "Jun":"6", "Jul":"7", "Aug":"8", "Sep":"9", "Oct":"10", "Nov":"11", "Dec":"12"}
    
    with open(resultDir + resultFileName) as result_data:    
            results = json.load(result_data)
    
    labels = list()
    for experiment in results:
            scanData = results[experiment]
            date = scanData["Date"]
            location = scanData["location"]
            amplitudes_for_u1 = scanData["polarizationU1"] # Got poitns for all experiments for polarization u1
            amplitudes_for_u9 = scanData["polarizationU9"] # Got poitns for all experiments for polarization u9
            amplitudes_for_uAVG = scanData["polarizationAVG"] # Got poitns for all experiments for polarization uAVG
            iter_number = scanData["Iteration_number"]
            
            label = "Station is " + location + "\n" + "Date is " + " ".join(date.split("_")) + "\n " + iter_number
            labels.append(label)
            
            dates = date.split("_")
            monthsNumber = dates[1]
            dates[1] = months[monthsNumber]
            date = " ".join(dates)
            key_u1 = date #+ " " + time
            #dateNumber = datetime.strptime(key_u1.strip(),  '%d %m %Y')
            
            result = Result(location, date, amplitudes_for_u1, amplitudes_for_u9, amplitudes_for_uAVG, iter_number)
            
            result_list.append(dict(result))
            
    result_list = sorted(result_list, key=itemgetter('iteration_number'), reverse=False)
    
    for erperiment in result_list:
        u1 = erperiment["polarizationU1"]
        u9 = erperiment["polarizationU9"]
        avg = erperiment["polarizationUAVG"]
        iteration_list.append(erperiment["iteration_number"])
        date_list.append(erperiment["date"])
        
        for i in u1:
            for vel in source_velocities:
                if float(vel) - velocities_range <= i[0] <= float(vel) + velocities_range:
                    velocitie_dict["u1"][vel].append(i[1]) 
        
        for j in u9:
            for vel in source_velocities:
                if float(vel) - velocities_range <= j[0] <= float(vel) + velocities_range:
                    velocitie_dict["u9"][vel].append(j[1]) 
                    
        for k in avg:
            for vel in source_velocities:
                if float(vel) - velocities_range <= k[0] <= float(vel) + velocities_range:
                    velocitie_dict["avg"][vel].append(k[1]) 
   
    x = list()
    for a in range(0, len(date_list)):
        x.append("Date " + date_list[a] + " Iteration number " + iteration_list[a])
    
    Symbols =  ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    
    fig = plt.figure()
    #fig, ax = plt.subplots()
    graph = fig.add_subplot(111)
    for i in range(0, len(source_velocities)-1):
        graph.plot(x, velocitie_dict["u1"][source_velocities[i]], Symbols[i]+"r", label="polarization U1 " + "Velocity " + source_velocities[i])
        graph.plot(x, velocitie_dict["u9"][source_velocities[i]], Symbols[i]+"g", label="polarization U9 " + "Velocity " + source_velocities[i])
        graph.plot(x, velocitie_dict["avg"][source_velocities[i]], Symbols[i]+"b", label="polarization AVG " + "Velocity " + source_velocities[i])
    
    plt.legend()
    
    cursor =  mplcursors.cursor(hover = True, highlight=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(labels[sel.target.index]))
   
    plt.show()
    
    sys.exit(0)
    
if __name__=="__main__":
    main()
    