#! /usr/bin/python

import os
import sys
import re
import time
import datetime

os.environ['TZ'] = 'UTC'
time.tzset()

timeFormat = "%Yy%jd%Hh%Mm%Ss"

logFile = dict()

def file_size(fname):
        statinfo = os.stat(fname)
        return statinfo.st_size
    
def usage():
    print ('Usage: '+sys.argv[0]+' log file')
    
class ExperimentLogReader():
    def __init__(self, logs):
        self.logs = logs
        self.scan_names = list()
        self.sources = list()
        self.dates = list()
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
        self.tcount = 0
        self.Location = ""
        
        self.logfile = open(self.logs, "r")
        self.datafile = open("m80ir.log".split(".")[0] + "sch.dat", "w")
        
        for x in range(0, file_size(self.logs)):
            logLine = self.logfile.readline()
        
            if "location" in logLine:
                self.Location = logLine.split(";")[1].split(",")[1].strip()
            
            if "scan_name" in logLine:
                self.scan_name = logLine.split(":")[3].split(",")[0].split("=")[1][2:].lstrip("0")
                self.scan_names.append(self.scan_name)
                
            if "source="  in logLine:
                self.source = logLine.split(":")[3].split(",")[0].split("=")[1] + "," + logLine.split(":")[3].split(",")[1]  + "," + logLine.split(":")[3].split(",")[2]
                self.sources.append(self.source)
                
                self.ra = list()
                self.dec = list()
                
                self.RA = logLine.split("=")[1].split(",")[1]
                self.DEC = logLine.split("=")[1].split(",")[2]
                self.Epoch = logLine.split("=")[1].split(",")[3]
                
                self.Epochs.append(self.Epoch)
                
                self.ra.append(self.RA[0:2])
                self.ra.append(self.RA[2:4])
                self.ra.append(self.RA[4:len(self.RA)])
                
                self.dec.append(self.DEC[0:2])
                self.dec.append(self.DEC[2:4])
                self.dec.append(self.DEC[4:len(self.RA)])
                
                self.RAs.append(self.ra)
                self.DECs.append(self.dec)
                
                self.ra = list()
                self.dec = list()
                
            if "scan_check" in logLine and "no" in logLine:
                datestring = logLine.split(",")[4].split(".")[0] + "s"
                tupletime = time.strptime(datestring, "%Yy%jd%Hh%Mm%Ss");
                year = tupletime.__getattribute__("tm_year")
                day = tupletime.__getattribute__("tm_mday")
                month = datetime.date(1900, int(tupletime.__getattribute__("tm_mon")) , 1).strftime('%B')[0:3]
                self.dates.append(str(day) + " " + month + " " + str(year))
                
            if "disk_record=on" in logLine:
                timeStart =  logLine.split(".")[2]
                self.timeStarts.append(timeStart)  
            
            if "disk_record=of" in logLine:
                timeStop =  logLine.split(".")[2]
                self.timeStops.append(timeStop)
                
            if "/tsys/" in logLine:
                self.tcount =  self.tcount + 1
                t = Systemtemperature = logLine.split("/")[2].split(",")[1]
                self.SystemtemperaturesForScan.append(t)
                    
                if  self.tcount == 2:
                    self.Systemtemperatures.append(self.SystemtemperaturesForScan)
                    self.tcount = 0
                    self.SystemtemperaturesForScan = list()
                        
            if "/bbc01=" in logLine:
                FreqBBC1 =  logLine.split("=")[1].split(",")[0]
                self.FreqBBC1s.append(FreqBBC1)
            
            if "/bbc02=" in logLine:
                FreqBBC2 =  logLine.split("=")[1].split(",")[0]
                self.FreqBBC2s.append(FreqBBC2)
            
            if "lo=loa" in logLine:
                loa =  logLine.split("=")[1].split(",")[1]
                self.loas.append(loa)
            
            if "lo=loc" in logLine:
                loc =  logLine.split("=")[1].split(",")[1]
                self.locs.append(loc)                          
    
        for y in range(0, len(self.scan_names)):
            DurMin =datetime.datetime.strptime(self.timeStops[y], "%H:%M:%S") -  datetime.datetime.strptime(self.timeStarts[y], "%H:%M:%S")
            self.DurationsMin.append(DurMin.seconds)
            self.DurationsSec.append(DurMin.seconds/60)
            
        for f in range(0, len(self.scan_names)):
            self.FreqStart.append(float(self.FreqBBC1s[f]) + float(self.loas[f]))
            self.FreqStop.append(float(self.FreqBBC2s[f])  + float(self.locs[f]))
            
    def writeOutput(self):
        self.datafile.write("Start;Header;")
        self.datafile.write("\n")
        self.datafile.write("Station;" + self.Location)
        self.datafile.write("\n")
    
        self.datafile.write("End;Header;----------------------;")
        self.datafile.write("\n")
    
        for scan in range (0, len(self.scan_names)):
            self.datafile.write("Scan;" + self.scan_names[scan] + ";")
            self.datafile.write("\n")
        
            self.datafile.write("Source;" + self.sources[scan] + ";")
            self.datafile.write("\n")
            
            self.datafile.write("Date;" + self.dates[scan] + ";")
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
            
            self.datafile.write("Systemtemperature1;1u;" + self.Systemtemperatures[scan][0] + ";")
            self.datafile.write("\n")
            self.datafile.write("Systemtemperature2;9u;" + self.Systemtemperatures[scan][1] + ";")
            self.datafile.write("\n")
            
            self.datafile.write("FreqBBC1;" + self.FreqBBC1s[scan] + ";Mhz")
            self.datafile.write("\n")
            
            self.datafile.write("FreqBBC2;" + self.FreqBBC2s[scan] + ";Mhz")
            self.datafile.write("\n")
            
            self.datafile.write("FreqStart;" + str(self.FreqStart[scan]) + ";Mhz")
            self.datafile.write("\n")
            
            self.datafile.write("FreqStop;" + str(self.FreqStop[scan]) + ";Mhz")
            self.datafile.write("\n")
            
            self.datafile.write("End;" + self.scan_names[scan] + ";----------------------;")
            self.datafile.write("\n")

    def getLgs(self):
        logs = dict()
        
        logs["location"] = self.Location
        
        for i in range(0, len(self.scan_names)):
            logs[self.scan_names[i]] = {"Systemtemperature":self.Systemtemperatures[i]}
            
        return logs
        
    def __del__(self):
        self.logfile.close()
        self.datafile.close()
     
if __name__=="__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
        
    experimentLogReader  = ExperimentLogReader(sys.argv[1])
    experimentLogReader.writeOutput()
    sys.exit(0)
        