# Methanol Maser Spectral Line Script parameterisfrom schedf ir file
# by K.Berzins 2017
# Requires python2
# Ver.1.3

def dummyLineRead (fileUnit,nr):
   for i in range(nr):
       line=fileUnit.readline()
       #i+=1
       #print i,": ",line
   return


import sys, os
if len(sys.argv[:]) > 1:
    irSchedFile = sys.argv[1]
else:
    irSchedFile="m39sch.ir"

#------------------------------------------------------------------------------------
print
print "Sched parameter file: \t",irSchedFile
if irSchedFile[-2:] != "ir":
     print "WARNING: The sched generated *.ir iparameter file was expected!" 

irSchedFileOUT = irSchedFile[:-2]+"dat"
print "OUT  parameter  file: \t",irSchedFileOUT

fpar = open(irSchedFile,'r')
fout = open(irSchedFileOUT,'w')

paramsOK = False
lineCount = 0
scanNr = 0

while paramsOK != True:
    line=fpar.readline()
    #print line
    #print line[:8]
    lineCount +=1
    if "Station:" in line:
        station = line[9:-2].strip()
        print "Station: ",station
    if "Experiment:" in line:
        experiment = line[12:].strip()
        print "Experiment: ",experiment
    if "Exp. Code:" in line:
        exp_code = line[11:].strip()
        print "Station: ",exp_code
    if line[:8] == "Start UT" or lineCount > 40: 
        paramsOK = True
    #line=fpar.readline()

print >> fout , "Start;Header;"
print >> fout , "Station;"+station+";"
print >> fout , "Experiment;"+experiment+";"
print >> fout , "Exp_Code;"+exp_code+";"
print >> fout , "End;Header;----------------------;"

dummyLineRead(fpar,3)

dateLine=fpar.readline()
dateLine=dateLine.strip()
datums = dateLine[9:20]
#print "DateLine:",dateLine
dummyLineRead(fpar,1)
line=fpar.readline() # Next scan frequencies

while "SETUP FILE INFORMATION" not in line: 

	if "Next scan frequencies" in line:
		scanNr +=1
		scanNrStr = str(scanNr).zfill(2)
		#print scanNrStr+": "+line.strip()
		lineSplit = line.split(':')
		dataLine = ' '.join(lineSplit[1].split())
		dataSplit = dataLine.split(' ')
		Freq1 = dataSplit[0] # MHz
		Freq2 = dataSplit[2] # Mhz
		#print Freq1, Freq2

		line=fpar.readline() # Next BBC frequencies
		#print line
		lineSplit = line.split(':')
		dataLine = ' '.join(lineSplit[1].split())
		dataSplit = dataLine.split(' ')
		FreqBBC1 = dataSplit[0] # MHz
		FreqBBC2 = dataSplit[2] # Mhz
		#print FreqBBC1, FreqBBC2

		line=fpar.readline() #
		#print "@@"+line
		if line.strip():
		    #print line.strip()
		    dummyLineRead(fpar,1)

		line=fpar.readline()	# START
		#print line.strip()
		Time1 = line[:8] 	# START-TIME UT (hh mm ss)
		Time1Str = Time1.replace(" ", "-")	# (hh-mm-ss)
		Name  = line[10:22].strip()		# SOURCE NAME
		EL0  = line[32:37].strip()		# ELEVATION on Start 
		AZ0  = line[38:43].strip()		# AZIMUTH on Start 
		#print Name, Time1Str

		line=fpar.readline()	# STOP
		#print line.strip()
		Time2 = line[:8]	# STOP-TIME UT (hh mm ss)
		Time2Str = Time2.replace(" ", "-")	# (hh-mm-ss)
		EL9  = line[32:37].strip()		# ELEVATION on End 
		AZ9  = line[38:43].strip()		# AZIMUTH on End 
		#print Name, Time1Str

		T1 = Time1Str.split('-')
		T2 = Time2Str.split('-')

		T1s = int(T1[0]) * 3600 + int(T1[1]) * 60 + int(T1[2])
		T2s = int(T2[0]) * 3600 + int(T2[1]) * 60 + int(T2[2])
		if T2s > T1s:
			duration = T2s - T1s 	# Exposition Duration in seconds
		else:
		        duration = T2s - T1s + 24*3600	# Asume Date changed during the scan
		# print duration



		#--------
		fpar2 = open(irSchedFile,'r')
		line2 = fpar2.readline()	# skip
		while "POSITIONS OF SOURCES USED IN RECORDING SCANS" not in line2:
			line2 = fpar2.readline()	# skip

		dummyLineRead(fpar2,4)
		while Name not in line2:
			line2 = fpar2.readline()	# Some Source data

		#print "$$",line2			# The needed source
		RA=line2[41:59].strip()
		line2 = fpar2.readline()		# The needed source Line Nr 2
		DEC=line2[41:59].strip()
		dummyLineRead(fpar2,2)
		line2 = fpar2.readline()		# The needed source Line Nr 2 with LSR velocities
		#line2 = ' '.join(line2.split())
                velLSR= line2.split()[0]		# Velocity at Local Standard of Rest (asume the same for all chanels)
                #print velLSR

		fpar2.close()
		#--------


	        print >> fout , "Scan;"+scanNrStr+";"
	        print >> fout , "Source;"+Name+";"
	        print >> fout , "Date;"+datums+";"
	        print >> fout , "Time;"+Time1.replace(" ", ":")+";UT;"
	        print >> fout , "TimeStop;"+Time2.replace(" ", ":")+";UT;"
	        print >> fout , "Duration;"+str(duration)+";sec;"+str(round(float(duration)/60,1))+";min"
	        print >> fout , "RA;"+RA+";"
	        print >> fout , "DEC;"+DEC+";"
	        print >> fout , "Epoch;2000;"
	        print >> fout , "ElevationStart;"+EL0+";deg;"
	        print >> fout , "AzimuthStart;"+AZ0+";deg;"
	        print >> fout , "ElevationStop;"+EL9+";deg;"
	        print >> fout , "AzimuthStop;"+AZ9+";deg;"
	        print >> fout , "FreqStart;"+Freq1+";MHz;"
	        print >> fout , "FreqBBC1;"+FreqBBC1+";MHz;"
	        print >> fout , "FreqStop;"+Freq2+";MHz;"
	        print >> fout , "FreqBBC2;"+FreqBBC2+";MHz;"
	        print >> fout , "VelocityLSR;"+velLSR+";km/s;"
	        print >> fout , "End;"+scanNrStr+";----------------------;"

	else:
		line=fpar.readline() # Try next line

fpar.close()
fout.close()

print "Data extraction: Done"
print
