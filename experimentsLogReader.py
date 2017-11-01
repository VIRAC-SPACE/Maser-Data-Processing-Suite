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
                self.Location = self.logLine.split(";")[1].split(",")[1].strip()
            
            if "scan_name" in logLine:
                self.scan_name = self.logLine.split(":")[3].split(",")[0].split("=")[1][2:].lstrip("0")
                self.scan_names.append(self.scan_name)
                
            if "source="  in logLine:
                self.source = self.logLine.split(":")[3].split(",")[0].split("=")[1] + "," + self.logLine.split(":")[3].split(",")[1]  + "," + self.logLine.split(":")[3].split(",")[2]
                self.sources.append(self.source)
                
                self.ra = list()
                self.dec = list()
                
                self.RA = self.logLine.split("=")[1].split(",")[1]
                self.DEC = self.logLine.split("=")[1].split(",")[2]
                self.Epoch = self.logLine.split("=")[1].split(",")[3]
                
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
                
            if "scan_check" in self.logLine and "no" in logLine:
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
        self.selfdatafile.write("\n")
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
    
if __name__=="__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    '''
    
    scan_names = list()
    sources = list()
    dates = list()
    timeStarts = list()
    timeStops = list()
    DurationsMin = list()
    DurationsSec = list()
    RAs = list()
    DECs = list()
    Epochs = list()
    Systemtemperatures = list()
    SystemtemperaturesForScan = list()
    FreqBBC1s = list()
    FreqBBC2s = list()
    FreqStart = list()
    FreqStop = list()
    loas = list()
    locs = list()
    tcount = 0
    Location = ""
    
    logfile = open(sys.argv[1], "r")
    datafile = open("m80ir.log".split(".")[0] + "sch.dat", "w") 
    
    for x in range(0, file_size(sys.argv[1])):
        logLine = logfile.readline()
        
        if "location" in logLine:
            Location = logLine.split(";")[1].split(",")[1].strip()
        
        if "scan_name" in logLine:
            scan_name = logLine.split(":")[3].split(",")[0].split("=")[1][2:].lstrip("0")
            scan_names.append(scan_name)
            
        if "source="  in logLine:
            source = logLine.split(":")[3].split(",")[0].split("=")[1] + "," + logLine.split(":")[3].split(",")[1]  + "," + logLine.split(":")[3].split(",")[2]
            sources.append(source)
            
            ra = list()
            dec = list()
            
            RA = logLine.split("=")[1].split(",")[1]
            DEC = logLine.split("=")[1].split(",")[2]
            Epoch = logLine.split("=")[1].split(",")[3]
            
            Epochs.append(Epoch)
            
            ra.append(RA[0:2])
            ra.append(RA[2:4])
            ra.append(RA[4:len(RA)])
            
            dec.append(DEC[0:2])
            dec.append(DEC[2:4])
            dec.append(DEC[4:len(RA)])
            
            RAs.append(ra)
            DECs.append(dec)
            
            ra = list()
            dec = list()
            
        if "scan_check" in logLine and "no" in logLine:
            datestring = logLine.split(",")[4].split(".")[0] + "s"
            tupletime = time.strptime(datestring, "%Yy%jd%Hh%Mm%Ss");
            year = tupletime.__getattribute__("tm_year")
            day = tupletime.__getattribute__("tm_mday")
            month = datetime.date(1900, int(tupletime.__getattribute__("tm_mon")) , 1).strftime('%B')[0:3]
            dates.append(str(day) + " " + month + " " + str(year))
            
        if "disk_record=on" in logLine:
            timeStart =  logLine.split(".")[2]
            timeStarts.append(timeStart)  
        
        if "disk_record=of" in logLine:
            timeStop =  logLine.split(".")[2]
            timeStops.append(timeStop)
            
        if "/tsys/" in logLine:
                tcount = tcount + 1
                t = Systemtemperature = logLine.split("/")[2].split(",")[1]
                SystemtemperaturesForScan.append(t)
                
                if tcount == 2:
                    Systemtemperatures.append(SystemtemperaturesForScan)
                    tcount = 0
                    SystemtemperaturesForScan = list()
                    
        if "/bbc01=" in logLine:
            FreqBBC1 =  logLine.split("=")[1].split(",")[0]
            FreqBBC1s.append(FreqBBC1)
        
        if "/bbc02=" in logLine:
            FreqBBC2 =  logLine.split("=")[1].split(",")[0]
            FreqBBC2s.append(FreqBBC2)
        
        if "lo=loa" in logLine:
            loa =  logLine.split("=")[1].split(",")[1]
            loas.append(loa)
        
        if "lo=loc" in logLine:
            loc =  logLine.split("=")[1].split(",")[1]
            locs.append(loc)                          
    
    for y in range(0, len(scan_names)):
        DurMin =datetime.datetime.strptime(timeStops[y], "%H:%M:%S") -  datetime.datetime.strptime(timeStarts[y], "%H:%M:%S")
        DurationsMin.append(DurMin.seconds)
        DurationsSec.append(DurMin.seconds/60)
        
    for f in range(0, len(scan_names)):
        FreqStart.append(float(FreqBBC1s[f]) + float(loas[f]))
        FreqStop.append(float(FreqBBC2s[f])  + float(locs[f]))
    
    datafile.write("Start;Header;")
    datafile.write("\n")
    datafile.write("Station;" + Location)
    datafile.write("\n")
    
    datafile.write("End;Header;----------------------;")
    datafile.write("\n")
    
    for scan in range (0, len(scan_names)):
        datafile.write("Scan;" + scan_names[scan] + ";")
        datafile.write("\n")
        
        datafile.write("Source;" + sources[scan] + ";")
        datafile.write("\n")
        
        datafile.write("Date;" + dates[scan] + ";")
        datafile.write("\n")
        
        datafile.write("TimeStart;" + timeStarts[scan] + ";" + "UT;")
        datafile.write("\n")
        datafile.write("TimeStop;" + timeStops[scan] + ";" + "UT;")
        datafile.write("\n")
        
        datafile.write("Duration;" + str(DurationsMin[scan]) + ";sec;" + str(DurationsSec[scan]) + ";min")
        datafile.write("\n")
        
        datafile.write("RA;" + " ".join(RAs[scan]) + ";")
        datafile.write("\n")
        
        datafile.write("DEC;" + " ".join(DECs[scan]) + ";")
        datafile.write("\n")

        datafile.write("Epoch;" + Epochs[scan] + ";")
        datafile.write("\n")
        
        datafile.write("Systemtemperature1;1u;" + Systemtemperatures[scan][0] + ";")
        datafile.write("\n")
        datafile.write("Systemtemperature2;9u;" + Systemtemperatures[scan][1] + ";")
        datafile.write("\n")
        
        datafile.write("FreqBBC1;" + FreqBBC1s[scan] + ";Mhz")
        datafile.write("\n")
        
        datafile.write("FreqBBC2;" + FreqBBC2s[scan] + ";Mhz")
        datafile.write("\n")
        
        datafile.write("FreqStart;" + str(FreqStart[scan]) + ";Mhz")
        datafile.write("\n")
        
        datafile.write("FreqStop;" + str(FreqStop[scan]) + ";Mhz")
        datafile.write("\n")
        
        datafile.write("End;" + scan_names[scan] + ";----------------------;")
        datafile.write("\n")
         
        
        
    for i in range(0, len(scan_names)):
        logFile[scan_names[i]] = {"Systemtemperature":Systemtemperatures[i]}

    print(logFile)
    logfile.close()
    datafile.close()
    '''
        
    ExperimentLogReader experimentLogReader = ExperimentLogReader(sys.argv[1])
    sys.exit(0)
        