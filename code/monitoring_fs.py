#! /usr/bin/python
import sys
import matplotlib.pyplot  as plt
from matplotlib.dates import date2num
import mplcursors
from datetime import datetime
import json
import argparse
import configparser

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
    
    resultFileName = source + ".json"
    
    months = {"Jan":"1", "Feb":"2", "Mar":"3", "Apr":"4", "May":"5", "Jun":"6", "Jul":"7", "Aug":"8", "Sep":"9", "Oct":"10", "Nov":"11", "Dec":"12"}
    
    with open(resultDir + resultFileName) as result_data:    
            result = json.load(result_data)
    
    velocitys_U1 = list()
    velocitys_U9 = list() 
    velocitys_AVG = list() 
    dateList = list()
    labels = list()
    for experiment in result:
            scanData = result[experiment]
            date = scanData["Date"]
            location = scanData["location"]
            amplitudes_for_u1 = scanData["polarizationU1"] # Got poitns for all experiments for polarization u1
            amplitudes_for_u9 = scanData["polarizationU9"] # Got poitns for all experiments for polarization u9
            amplitudes_for_uAVG = scanData["polarizationAVG"] # Got poitns for all experiments for polarization uAVG
            
            label = "Station is " + location + "\n" + "Date is " + " ".join(date.split("_"))  
            labels.append(label)
            
            v_list_u1 = list()
            for v in amplitudes_for_u1:
                v_list_u1.append(v[1])
            velocitys_U1.append(v_list_u1)
            
            v_list_u9 = list()
            for v in amplitudes_for_u9:
                v_list_u9.append(v[1])
            velocitys_U9.append(v_list_u9)
            
            v_list_AVG = list()
            for v in amplitudes_for_uAVG:
                v_list_AVG.append(v[1])
            velocitys_AVG.append(v_list_AVG)
            
            dates = date.split("_")
            monthsNumber = dates[1]
            dates[1] = months[monthsNumber]
            date = " ".join(dates)
            key_u1 = date #+ " " + time
            dateNumber = datetime.strptime(key_u1.strip(),  '%d %m %Y')
            dateList.append(dateNumber)
    
    y_u1 = list()
    velocityCount = len(velocitys_U1[0])
    dummy = 0
    while dummy != velocityCount -1:
        velo = list()
        for vel in velocitys_U1:
            velo.append(vel[dummy])
        y_u1.append(velo)
        dummy = dummy + 1
        
    y_u9 = list()
    velocityCount = len(velocitys_U9[0])
    dummy = 0
    while dummy != velocityCount -1:
        velo = list()
        for vel in velocitys_U9:
            velo.append(vel[dummy])
        y_u9.append(velo)
        dummy = dummy + 1
    
    y_avg = list()
    velocityCount = len(velocitys_AVG[0])    
    dummy = 0
    while dummy != velocityCount -1:
        velo = list()
        for vel in velocitys_AVG:
            velo.append(vel[dummy])
        y_avg.append(velo)
        dummy = dummy + 1
    
    x = [date2num(date) for date in  dateList]
    
    Symbols =  ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    
    fig = plt.figure()
    #fig, ax = plt.subplots()
    graph = fig.add_subplot(111)
    for i in range(0, len(y_u1)):
        graph.plot(x, y_u1[i], Symbols[i]+"r", label="polarization U1 " + "Velocity " + source_velocities[i])
        graph.plot(x, y_u9[i], Symbols[i]+"g", label="polarization U9 " + "Velocity " + source_velocities[i])
        graph.plot(x, y_avg[i], Symbols[i]+"b", label="polarization AVG " + "Velocity " + source_velocities[i])
    graph.set_xticks(x)
    graph.set_xticklabels([date.strftime("%d %m %Y") for date in  dateList])
    plt.legend()
    
    cursor =  mplcursors.cursor(hover = True, highlight=True)
    cursor.connect("add", lambda sel: sel.annotation.set_text(labels[sel.target.index]))
   
    plt.show()
    
    sys.exit(0)

if __name__=="__main__":
    main()
    