#! /usr/bin/python
import os
import sys
import time
import datetime
import argparse
import platform

import yaml
from scan import Scan
    
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Reads log file and create pretty logs. ''',
    epilog="""LOGGREADER.""")

    # Positional mandatory arguments
    parser.add_argument("logFile", help="Experiment log file name", type=str)

    # Optional arguments
    parser.add_argument("-c", "--config", help="Configuration Yaml file", type=str, default="config/logConfig.yaml")
    parser.add_argument("-s", "--source", help="Set RA, DEC, Epoch, Source name", nargs="*", type=str, default=[])
    # option -s example 

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 3.0')

    # Parse arguments
    args = parser.parse_args()

    return args
    
def usage():
    print ('Usage:   log file')
    
class ExperimentLogReader():
    def __init__(self, logs, prettyLogs, singleSourceName=None):
        self.logs = logs
        self.prettyLogs = prettyLogs
        self.singleSourceName = singleSourceName
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
        self.FreqBBC1s = list()
        self.FreqBBC2s = list()
        self.FreqStart = list()
        self.FreqStop = list()
        self.loas = list()
        self.locs = list()
        self.scanNameString = list()
        self.sourceName = list()
        self.clocks = list()
        self.scanList = list()
        self.headerLines = list()
        self.scanLines = dict()
        self.Location = ""
        
        if len(self.singleSourceName) != 0:
            self.single = True
        else:
            self.single = False
    
        try:
            self.logfile = open(self.logs, "r")
            self.datafile = open(self.prettyLogs + self.logs.split(".")[0].split("/")[1] + "log.dat", "w")
            
        except(IOError):
            print "IO Error"
            
        except:
            print("Unexpected error:", sys.exc_info()[0])
            raise
            
        else:
            append = False
            key = 0
            previousScan = 0
            
            for line in self.logfile.readlines():
            
                if "location" in line:
                    self.Location = line.split(";")[1].split(",")[1].strip()
                    year = line.split(".")[0]
                    dayNumber = line.split(".")[1]
                    date = datetime.datetime(int(year), 1, 1) + datetime.timedelta(int(dayNumber) - 1)
                    day = date.day
                    monthNr = date.month
                    month = datetime.date(1900, int(monthNr) , 1).strftime('%B')[0:3]
                    
                    self.dates = str(day).zfill(2) + " " + month + " " + str(year)
                
                elif "scan_name=no" in line:
                    append = True 
                    self.scan_name = line.split(":")[3].split(",")[0].split("=")[1][2:].lstrip("0")
                    key = int(self.scan_name)
                    
                    if self.scan_name in self.scan_names:
                        raise Exception("Two scans with same name " + self.scan_name)
                        
                    self.scan_names.append(self.scan_name)
                    self.scanLines[key] = list()
                    
                    if  previousScan !=0 and previousScan +1 != key:
                        print "Skipped scan ", previousScan + 1
                        
                    previousScan = key
                
                #Testing if line is not in header and it is not empty   
                if len(line) != 0 and append:
                    self.scanLines[key].append(line)
                
                #Testing if line is  in header and it is not empty    
                elif len(line) != 0 and append==False:
                    self.headerLines.append(line)
                    
            header = Scan(self.headerLines)
            header.getParametrs()
            header_source, header_sourceName, header_epoch, header_ra, header_dec, header_timeStart, header_timeStop, header_SystemtemperaturesForScan, header_freqBBC1, header_freqBBC2, header_loa, header_loc, header_clock = header.returnParametrs()
            
            for scan in self.scanLines:
                scanData = Scan(self.scanLines[scan])
                self.scanList.append(scanData)
                scanData.setScanNumber(scan)
                scanData.getParametrs()
                source, sourceName, epoch, ra, dec, timeStart, timeStop, SystemtemperaturesForScan, freqBBC1, freqBBC2, loa, loc, clock = scanData.returnParametrs()
                
                if self.single:
                    source =  self.singleSourceName[0] + "," + self.singleSourceName[1] + "," + self.singleSourceName[2] 
                    ra = list()
                    dec = list()
                    Ra = self.singleSourceName[1]
                    Dec = self.singleSourceName[2]
                    epoch =  self.singleSourceName[3]
                    
                    ra.append(Ra[0:2])
                    ra.append(Ra[2:4])
                    ra.append(Ra[4:len(Ra)])
                    
                    if Dec[0] == "-":
                        dec.append(Dec[0:3])
                        dec.append(Dec[3:5])
                        dec.append(Dec[5:len(self.Dec)])
                    else:    
                        dec.append(Dec[0:2])
                        dec.append(Dec[2:4])
                        dec.append(Dec[4:len(Dec)])
                    
                    freqBBC2 = header_freqBBC2
                    loc = header_loc
                    
                self.sources.append(source)
                self.sourceName.append(sourceName)
                self.Epochs.append(epoch)
                self.RAs.append(ra)
                self.DECs.append(dec)
                self.timeStarts.append(timeStart)
                self.timeStops.append(timeStop)
                self.Systemtemperatures.append(SystemtemperaturesForScan)
                self.FreqBBC1s.append(freqBBC1)
                self.FreqBBC2s.append(freqBBC2)
                self.loas.append(loa)
                self.locs.append(loc)
                self.clocks.append(clock)
                            
            for y in range(0, len(self.scan_names)):
                try:
                    DurMin = datetime.datetime.strptime(self.timeStops[y], "%H:%M:%S") -  datetime.datetime.strptime(self.timeStarts[y], "%H:%M:%S")
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
            self.logfile.close()
           
    def writeOutput(self):
        self.datafile.write("Start;Header;")
        self.datafile.write("\n")
        self.datafile.write("Station;" + self.Location)
        self.datafile.write("\n")
        
        if any(scan.getmanualyChangedSystemTemU1() or scan.getmanualyChangedSystemTemU9() or scan.getmanualyChangedBBC1() or scan.getmanualyChangedBBC2() for scan in self.scanList):
            self.datafile.write("Manual Changes !!!" )
            self.datafile.write("\n")
        
        for scan in self.scanList:
            if scan.getmanualyChangedSystemTemU1():
                self.datafile.write("For scan Number " + str(scan.getScanNumber()) + " Manually Changed System Temperature U1")
                self.datafile.write("\n")
                 
            if scan.getmanualyChangedSystemTemU9():
                self.datafile.write("For scan Number " + str(scan.getScanNumber()) + " Manually Changed System Temperature U9")
                self.datafile.write("\n")
                
            if scan.getmanualyChangedBBC1():
                self.datafile.write("For scan Number " + str(scan.getScanNumber()) + " Manually Changed BBC1")
                self.datafile.write("\n")
                
            if scan.getmanualyChangedBBC2():
                self.datafile.write("For scan Number " + str(scan.getScanNumber()) + " Manually Changed BBC2")
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
    
            self.datafile.write("Epoch;" + str(self.Epochs[scan]) + ";")
            self.datafile.write("\n")
            
            self.datafile.write("Systemtemperature1;1u;" + str(self.Systemtemperatures[scan][0]) + ";")
            self.datafile.write("\n")
            self.datafile.write("Systemtemperature2;9u;" + str(self.Systemtemperatures[scan][1]) + ";")
            self.datafile.write("\n")
            
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
            
        self.datafile.close()  
        print ("Created file " + "prettyLogs/" + self.logs.split(".")[0].split("/")[1] + "log.dat")

    def getLogs(self):
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
        del self.scan_names
        del self.sources
        del self.timeStarts
        del self.timeStops
        del self.DurationsMin
        del self.DurationsSec
        del self.RAs
        del self.DECs
        del self.Epochs
        del self.Systemtemperatures
        del self.FreqBBC1s
        del self.FreqBBC2s
        del self.FreqStart
        del self.FreqStop
        del self.loas
        del self.locs
        del self.scanNameString
        del self.sourceName
        del self.clocks
        del self.scanList
        del self.scanLines
           
def main():     
    if platform.system() == "Linux":
        os.environ['TZ'] = 'UTC'
        time.tzset()
    
    elif platform.system() == "Windows":
        os.environ['TZ'] = 'UTC'
        pass
     
    # Parse the arguments
    args = parseArguments()
    
    logFileName = str(args.__dict__["logFile"])
    configFilePath = str(args.__dict__["config"])
    singleSourceExperiment = list(args.__dict__["source"])
    
    #Creating config parametrs
    config = dict()
    for key, value in yaml.load(open(configFilePath)).iteritems():
        config[key] = value

    logPath =  config["logPath"]
    prettyLogsPath = config["prettyLogsPath"]
    
    experimentLogReader = ExperimentLogReader(logPath + logFileName, prettyLogsPath, singleSourceExperiment)
    experimentLogReader.writeOutput()
    sys.exit(0)
    
if __name__=="__main__":
    main()
    
    
