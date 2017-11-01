#! /usr/bin/python

import os
import sys
import re
import time
import datetime

os.environ['TZ'] = 'UTC'
time.tzset()

timeFormat = "%Yy%jd%Hh%Mm%Ss"

def file_size(fname):
        statinfo = os.stat(fname)
        return statinfo.st_size
    
def usage():
    print ('Usage: '+sys.argv[0]+' log file')
    
def FileCheck(fn):
    try:
        open(fn, "r")
        return 1
    except IOError:
        print "Error: File does not appear to exist."
        return 0

if __name__=="__main__":
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)
    
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
     
    logfile.close()
    datafile.close()
    sys.exit(0)
        