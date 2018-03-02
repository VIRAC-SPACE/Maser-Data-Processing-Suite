#! /usr/bin/python

import os
import sys
import time
import datetime
import argparse
from decimal import Decimal
import collections

import yaml

#nolasit /gps-fmout/

os.environ['TZ'] = 'UTC'
time.tzset()

timeFormat = "%Yy%jd%Hh%Mm%Ss"

logFile = dict()

def file_size(fname):
        statinfo = os.stat(fname)
        return statinfo.st_size
    
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Reads log file and create pretty logs. ''',
    epilog="""LOGGREADER.""")

    # Positional mandatory arguments
    parser.add_argument("logFile", help="Experiment log file name", type=str)

    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration Yaml file", type=str, default="config/logConfig.yaml")

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 1.0')

    # Parse arguments
    args = parser.parse_args()

    return args
    
def usage():
    print ('Usage:   log file')
    
class ExperimentLogReader():
    def __init__(self, logs, prettyLogs):
        self.logs = logs
        self.prettyLogs = prettyLogs
        self.scan_names = list()
        self.sources = list()
        self.dates = ""
        self.timeStarts = list()
        self.timeStops = list()
        self.DurationsMin = list()
        self.DurationsSec = list()
        self.RAs = list()
        self.DECs = list()
        self.Epochs = list()
        self.Systemtemperatures = list()
        self.SystemtemperaturesForScan = list()
        self.FreqBBC1s = list()
        self.FreqBBC2s = list()
        self.FreqStart = list()
        self.FreqStop = list()
        self.loas = list()
        self.locs = list()
        self.scanNameString = list()
        self.sourceName = list()
        self.clocks = list()
        self.scanLines = dict()
        self.tcount = 0
        self.Location = ""
        self.bbc1count = 0
        self.bbc2count = 0
        self.scan_count = 0
        self.loacount = 0
        self.loccount = 0
        self.sourcecount = 0
        self.datecount = 0
        self.clockcount = 0
        
        self.logfile = open(self.logs, "r")
        self.datafile = open(self.prettyLogs + self.logs.split(".")[0].split("/")[1] + "log.dat", "w")
        
        append = False
        key = 0
        for x in range(0, file_size(self.logs)):
            logLine = self.logfile.readline()
            
            if "location" in logLine:
                self.Location = logLine.split(";")[1].split(",")[1].strip()
                year = logLine.split(".")[0]
                dayNumber = logLine.split(".")[1]
                date = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(dayNumber) - 1)
                day = date.day
                monthNr = date.month
                month = datetime.date(1900, int(monthNr) , 1).strftime('%B')[0:3]
                
                self.dates = str(day).zfill(2) + " " + month + " " + str(year)
            
            elif "scan_name=no" in logLine:
                append = True 
                self.scan_count = self.scan_count + 1
                self.scan_name = logLine.split(":")[3].split(",")[0].split("=")[1][2:].lstrip("0")
                key = int(self.scan_name)
                self.scan_names.append(self.scan_name)
                self.scanLines[key] = list()
                
            if len(logLine) != 0 and append:
                self.scanLines[key].append(logLine)
                    
        for key in self.scanLines:
            if  self.sourcecount != 0:
                self.sourcecount = self.sourcecount -1  
            for line in self.scanLines[key]:
                if "source="  in line:
                    self.sourcecount =  self.sourcecount + 1
                    print "key", key, "self.sourcecount", self.sourcecount 
                    
                    #print self.scan_count
                    if self.sourcecount != key:
                        print True
                        self.RAs.append(['0','0','0'])
                        self.DECs.append(['0','0','0'])
                        self.sources.append("None")
                        self.sourceName.append("None")
                        self.Epochs.append('0')
                        self.sourcecount =  self.sourcecount + 1

                    if line.endswith(",\n"):
                        logLineSplit = line[:-2].split(",")
                    else:
                        logLineSplit = line[:-1].split(",")
                
                    self.source = logLineSplit[0].split("=")[1]+","+logLineSplit[-3]+","+logLineSplit[-2]
                    
                    self.sources.append(self.source)
                    self.sourceName.append(logLineSplit[0].split("=")[1])
               
                    self.ra = list()
                    self.dec = list()
                    
                    self.RA = self.source.split(",")[1]
                    self.DEC = self.source.split(",")[2]
                    self.Epoch = logLineSplit[-1]
                    self.Epochs.append(self.Epoch)
                    
                    self.ra.append(self.RA[0:2])
                    self.ra.append(self.RA[2:4])
                    self.ra.append(self.RA[4:len(self.RA)])
                    
                    if self.DEC[0] == "-":
                        self.dec.append(self.DEC[0:3])
                        self.dec.append(self.DEC[3:5])
                        self.dec.append(self.DEC[5:len(self.DEC)])
                    else:    
                        self.dec.append(self.DEC[0:2])
                        self.dec.append(self.DEC[2:4])
                        self.dec.append(self.DEC[4:len(self.DEC)])   
                    
                    self.RAs.append(self.ra)
                    self.DECs.append(self.dec)
                    
                    self.ra = list()
                    self.dec = list()
                    
        print self.sources
        '''    
            elif "disk_record=on" in logLine:
                timeStart =  logLine.split(".")[2]
                self.timeStarts.append(timeStart)  
            
            elif "disk_record=of" in logLine:
                timeStop =  logLine.split(".")[2]
                self.timeStops.append(timeStop)
                
            elif "/tsys/" in logLine:
                self.tcount =  self.tcount + 1
                t  = logLine.split("/")[2].split(",")[1]
                self.SystemtemperaturesForScan.append(t)
                    
                if  self.tcount == 2:
                    self.Systemtemperatures.append(self.SystemtemperaturesForScan)
                    self.tcount = 0
                    self.SystemtemperaturesForScan = list()
                        
            elif "/bbc01=" in logLine:
                self.bbc1count =  self.bbc1count + 1
                
                if self.scan_count != self.bbc1count:
                    self.FreqBBC1s.append(0)
                    self.bbc1count =  self.bbc1count + 1
                
                FreqBBC1 =  logLine.split("=")[1].split(",")[0]
                self.FreqBBC1s.append(FreqBBC1)
            
            elif "/bbc02=" in logLine:
                self.bbc2count =  self.bbc2count + 1
                
                if self.scan_count != self.bbc2count:
                    self.FreqBBC2s.append(0)
                    self.bbc2count =  self.bbc2count + 1
                    
                FreqBBC2 =  logLine.split("=")[1].split(",")[0]
                self.FreqBBC2s.append(FreqBBC2)
            
            elif "lo=loa" in logLine:
                self.loacount = self.loacount + 1
                
                if self.scan_count != self.loacount:
                    self.loas.append(0)
                    self.loacount =  self.loacount + 1
                loa =  logLine.split("=")[1].split(",")[1]
                self.loas.append(loa)
            
            elif "lo=loc" in logLine:
                self.loccount = self.loccount + 1
                
                if self.scan_count != self.loccount:
                    self.locs.append(0)
                    self.loccount =  self.loccount + 1
                
                loc =  logLine.split("=")[1].split(",")[1]
                self.locs.append(loc)
            
            elif "/gps-fmout/" in logLine:
                self.clockcount = self.clockcount + 1
                
                if self.scan_count != self.clockcount:
                    self.clocks.append(0)
                    self.clockcount = self.clockcount + 1
                
                self.clocks.append(Decimal(logLine.split("/")[2]) * 10000)
                
            #elif  "/bbc01/" in logLine:
                #self.bbc1count =  self.bbc1count + 1
                #if self.bbc1count == 1:
                    #bbc1 = logLine.split("/")[2].split(",")[0]
                    #self.FreqBBC1s.append(bbc1)
                #elif self.bbc1count == 6:
                    #self.bbc1count = 0
            
            #elif  "/bbc02/" in logLine:
                #print(logLine)
                #self.bbc2count =  self.bbc2count + 1
                #if self.bbc2count == 1:
                    #bbc2 = logLine.split("/")[2].split(",")[0]
                    #self.FreqBBC2s.append(bbc2)
                #elif self.bbc2count == 6:
                    #self.bbc2count = 0
        
        #if len(self.scan_names) > len(self.FreqBBC1s):
        while len(self.scan_names) > len(self.FreqBBC1s):
            self.FreqBBC1s.append(0)
            
        #if len(self.scan_names) > len(self.FreqBBC2s):
        while len(self.scan_names) > len(self.FreqBBC2s):
            self.FreqBBC2s.append(0)
            
        #if len(self.scan_names) > len(self.loas):
        while len(self.scan_names) > len(self.loas):
            self.loas.append(0)
        
        #if len(self.scan_names) > len(self.locs):
        while len(self.scan_names) > len(self.locs):
            self.locs.append(0)

        while len(self.scan_names) > len(self.timeStops):
            self.timeStops.append('0')

        while len(self.scan_names) > len(self.timeStarts):
            self.timeStarts.append('0')
                        
        for y in range(0, len(self.scan_names)):
            try:
                DurMin =datetime.datetime.strptime(self.timeStops[y], "%H:%M:%S") -  datetime.datetime.strptime(self.timeStarts[y], "%H:%M:%S")
                self.DurationsMin.append(DurMin.seconds)
                self.DurationsSec.append(DurMin.seconds/60)
            except:
                self.DurationsMin.append(0)
                self.DurationsSec.append(0)
                continue
            
        for f in range(0, len(self.scan_names)):
            self.FreqStart.append(float(self.FreqBBC1s[f]) + float(self.loas[f]))
            self.FreqStop.append(float(self.FreqBBC2s[f])  + float(self.locs[f]))
        
        #print(len(self.FreqStart), " ", len(self.scan_names), " ", len(self.FreqBBC1s), " ", len(self.FreqBBC2s))
        '''
                  
    def writeOutput(self):
        self.datafile.write("Start;Header;")
        self.datafile.write("\n")
        self.datafile.write("Station;" + self.Location)
        self.datafile.write("\n")
    
        self.datafile.write("End;Header;----------------------;")
        self.datafile.write("\n")
    
        for scan in range (0, len(self.scan_names)):
            self.datafile.write("Scan;" + self.scan_names[scan].zfill(2) + ";")
            self.datafile.write("\n")
        
            self.datafile.write("Source;" + self.sources[scan] + ";")
            self.datafile.write("\n")
            
            self.datafile.write("Date;" + self.dates + ";")
            self.datafile.write("\n")
            
            self.datafile.write("TimeStart;" + self.timeStarts[scan] + ";" + "UT;")
            self.datafile.write("\n")
            self.datafile.write("TimeStop;" + self.timeStops[scan] + ";" + "UT;")
            self.datafile.write("\n")
            
            self.datafile.write("Duration;" + str(self.DurationsMin[scan]) + ";sec;" + str(self.DurationsSec[scan]) + ";min")
            self.datafile.write("\n")
            
            self.datafile.write("RA;" + " ".join(self.RAs[scan]) + ";")
            self.datafile.write("\n")
            
            self.datafile.write("DEC;" + " ".join(self.DECs[scan]) + ";")
            self.datafile.write("\n")
    
            self.datafile.write("Epoch;" + self.Epochs[scan] + ";")
            self.datafile.write("\n")
            
            try:
                self.datafile.write("Systemtemperature1;1u;" + self.Systemtemperatures[scan][0] + ";")
                self.datafile.write("\n")
                self.datafile.write("Systemtemperature2;9u;" + self.Systemtemperatures[scan][1] + ";")
                self.datafile.write("\n")
            except:
                pass
            
            self.datafile.write("FreqBBC1;" + str(self.FreqBBC1s[scan]) + ";MHz")
            self.datafile.write("\n")
            
            self.datafile.write("FreqBBC2;" + str(self.FreqBBC2s[scan]) + ";MHz")
            self.datafile.write("\n")
            
            self.datafile.write("FreqStart;" + str(self.FreqStart[scan]) + ";MHz")
            self.datafile.write("\n")
            
            self.datafile.write("FreqStop;" + str(self.FreqStop[scan]) + ";MHz")
            self.datafile.write("\n")
            
            self.datafile.write("ClockOffset;" + str(self.clocks[scan]))
            self.datafile.write("\n")
            
            self.datafile.write("End;" + self.scan_names[scan].zfill(2) + ";----------------------;")
            self.datafile.write("\n")
            
        print ("Created file " + "prettyLogs/" + self.logs.split(".")[0].split("/")[1] + "log.dat")

    def getLgs(self):
        self.writeOutput()
        logs = dict()
        
        logs["location"] = self.Location
        
        for i in range(0, len(self.scan_names)):

            logs[self.scan_names[i]] = {"Systemtemperature":self.Systemtemperatures[i], "Ra":self.RAs[i] , "Dec":self.DECs[i], "dates":self.dates, "startTime":self.timeStarts[i], "FreqStart": self.FreqStart[i], "sourceName":self.sources[i], "source":self.sourceName[i], "stopTime": self.timeStops[i], "clockOffset": self.clocks[i]}

        return logs
    
    def getAllScansNumbers(self):
        return  self.scan_names
    
    def getScansForSource(self, sourceName):
        indices_source = [i for i, x in enumerate(self.sourceName) if x == sourceName]
        scanNamesForSource = list()
        
        for j in range(0,len(indices_source)):
            scanNamesForSource.append(self.scan_names[indices_source[j]])
        
        return scanNamesForSource
        
    def __del__(self):
        self.logfile.close()
        self.datafile.close()
     
if __name__=="__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    
    # Parse the arguments
    args = parseArguments()
    
    logFileName = str(args.__dict__["logFile"])
    configFilePath = str(args.__dict__["config"])
    
    #Creating config parametrs
    logPath =  ""
    prettyLogsPath = ""
    config = dict()
    for key, value in yaml.load(open(configFilePath)).iteritems():
        config[key] = value

    logPath =  config["logPath"]
    prettyLogsPath = config["prettyLogsPath"]
       
    experimentLogReader = ExperimentLogReader(logPath + logFileName, prettyLogsPath)
    #experimentLogReader.writeOutput()
    experimentLogReader.__del__()
    sys.exit(0)
        
