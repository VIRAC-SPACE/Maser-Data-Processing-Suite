#! /usr/bin/python

import os
import sys
import re
import time
import datetime


#nolasit /gps-fmout/

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
        self.bbc1count = 0
        self.bbc2count = 0
        self.scan_count = 0
        self.loacount = 0
        self.loccount = 0
        self.sourcecount = 0
        self.datecount = 0
        
        self.logfile = open(self.logs, "r")
        self.datafile = open("prettyLogs/" + self.logs.split(".")[0].split("/")[1] + "sch.dat", "w")
        
        #print self.datafile
        
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
                
                self.dates.append(str(day).zfill(2) + " " + month + " " + str(year))
            
            elif "scan_name=no" in logLine:
                self.scan_count = self.scan_count + 1
                self.scan_name = logLine.split(":")[3].split(",")[0].split("=")[1][2:].lstrip("0")
                self.scan_names.append(self.scan_name)
                
            elif "source="  in logLine:
                self.sourcecount =  self.sourcecount + 1
                
                if self.scan_count != self.sourcecount:
                    self.RAs.append(['0','0','0'])
                    self.DECs.append(['0','0','0'])
                    self.sources.append("None")
                    self.Epochs.append('0')
                    self.sourcecount =  self.sourcecount + 1
                #self.source = logLine.split(":")[3].split(",")[0].split("=")[1] + "," + logLine.split(":")[3].split(",")[1]  + "," + logLine.split(":")[3].split(",")[2]
                if logLine.endswith(",\n"):
                    logLineSplit = logLine[:-2].split(",")
                else:
                    logLineSplit = logLine[:-1].split(",")
                self.source = logLineSplit[0].split("=")[1]+","+logLineSplit[-3]+","+logLineSplit[-2]
                self.sources.append(self.source)
                #print logLine
                #print self.source
                
                self.ra = list()
                self.dec = list()
                
                #self.RA = logLine.split("=")[1].split(",")[1]
                #self.DEC = logLine.split("=")[1].split(",")[2]
                #self.Epoch = logLine.split("=")[1].split(",")[3]
                self.RA = self.source.split(",")[1]
                self.DEC = self.source.split(",")[2]
                self.Epoch = logLineSplit[-1]
                
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
                #print(len(self.scan_names), " ", len(self.timeStops), " ", len(self.timeStarts))
                #print self.scan_names
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
            
            self.datafile.write("Date;" + self.dates[0] + ";")
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
            
            self.datafile.write("End;" + self.scan_names[scan].zfill(2) + ";----------------------;")
            self.datafile.write("\n")

    def getLgs(self):
        self.writeOutput()
        logs = dict()
        
        logs["location"] = self.Location
        
        for i in range(0, len(self.scan_names)):
            logs[self.scan_names[i]] = {"Systemtemperature":self.Systemtemperatures[i], "Ra":self.RAs[i] , "Dec":self.DECs[i], "dates":self.dates[0], "startTime":self.timeStarts[i], "FreqStart": self.FreqStart[i], "sourceName":self.sources[i]}

        return logs
    
    def getAllScansNumbers(self):
        return  self.scan_names
        
    def __del__(self):
        self.logfile.close()
        self.datafile.close()
     
if __name__=="__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
        
    experimentLogReader = ExperimentLogReader("logs/" + sys.argv[1])
    experimentLogReader.writeOutput()
    experimentLogReader.__del__()
    sys.exit(0)
        