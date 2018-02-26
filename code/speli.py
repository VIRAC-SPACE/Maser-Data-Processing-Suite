# Methanol Maser Spectral Line Scrip
# by K.Berzins & A.Aberfelds 2016
# Requires python2, dopsetpy.f
# Ver.3.31 alpha

#def FWHM(X,Y):
#    half_max = max(Y) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
#    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    #plot(X,d) #if you are interested
    #find the left and right most indexes
#    left_idx = find(d > 0)[0]
#    right_idx = find(d < 0)[-1]
#    return X[right_idx] - X[left_idx] #return the difference (full width) 

from scipy.interpolate import splrep, sproot, splev
import scipy.constants
#import scipy.interpolate

#Relative unit -> K -Jy transformation calculation----------------------------------
def calibration1(Station,Tsys):
    if "IRBENE16" in stationStr:
        scale1=26
    elif "IRBENE" in stationStr:
        scale1=12
    scale1u=scale1*Tsys_1u
    return scale1u

def calibration2(Station,Tsys):
    if "IRBENE16" in stationStr:
        scale1=26
    elif "IRBENE" in stationStr:
        scale1=12
    scale9u=scale1*Tsys_9u
    return scale9u      
    
#----------------------------------------------------------------------------------    
  
def FWHM(x, y, k=10):
    """
    Determine full-with-half-maximum of a peaked set of points, x and y.

    Assumes that there is only one peak present in the datasset.  The function
    uses a spline interpolation of order k.
    """

    class MultiplePeaks(Exception): pass
    class NoPeaksFound(Exception): pass

    half_max = max(y)/2.0
    print half_max
    nknots = k
    knots = np.arange(x[1],x[len(x)-1],(x[len(x)-1]-x[1])/np.double(nknots))
    s = splrep(x, y - half_max, t=knots)
    roots = sproot(s)

    if len(roots) > 2:
        raise MultiplePeaks("The dataset appears to have multiple peaks, and "
                "thus the FWHM can't be determined.")
    if len(roots) < 2:
        raise NoPeaksFound("No proper peaks were found in the data set; likely "
                "the dataset is flat (e.g. all zeros).")
    else:
        return abs(roots[1] - roots[0])
        
        
#radial velocity calculation
def freq2vel (ObservedFrequency,VelReceiver):
    c = scipy.constants.speed_of_light
    f0 = 6668519200 # Hz 
    VelSource = (c*(f0-ObservedFrequency)/f0)/1000+VelReceiver
    return VelSource
 
def rangex (X,Y1,Y2,X0,deltaX):
    NN = len(X)
    print "NN=",NN
    i1 = 0
    i2 = NN-1
    XX0 = []
    YY1 = []
    YY2 = []
    for i in range (0,NN-1):
        if X[i] >= (X0 - deltaX) and X[i] <= (X0 +deltaX) :
           XX0.append(X[i])
           YY1.append(Y1[i])
           YY2.append(Y2[i])
    return XX0,YY1,YY2
#________________________________________________________________
#Allovs to select points from plot
class MyClick(object):
    def __init__(self):
        # initialize mode
        self._mode = 1 
        self.peak = []
        self.pk = [] # new - list where to store artists/pts already picked
        
    def press(self, event):
        """record last key pressed"""
        sys.stdout.flush()
        if event.key == '1':   
            self._mode = 1 
            print('Peak Mode')

        elif event.key == '4':
            self._mode = 4
            print('Delete Mode')

	elif event.key == '5' or event.key=='q' or event.key=='Q':
	    self._mode = 0
	    print('LOG_EVENT: Exit editing mode.')
	    plt.close()

    def pick(self,event):
        mouseevent = event.mouseevent
        artist = event.artist
        x, y = np.array(artist.get_xdata()), np.array(artist.get_ydata()) #used numpy array for argsort feature
        ind = event.ind

        #remove = [artist for artist in self.peak if artist.contains()]

        if self._mode == 1:
            # peak
            peakind = y[ind].argsort()[-1:] #argsort to find peak
            self.peak.append(zip(x[ind[peakind]], y[ind[peakind]])) #save peak
            pt, = ax.plot(x[ind[peakind]],y[ind[peakind]],'ro',picker=5) #new
            self.pk.append(ind[peakind])  #new - store into a list
            print('peak: ', zip(x[ind[peakind]], y[ind[peakind]]))
            pass

        elif self._mode == 4:
            # delete
            for i in self.pk:   # this loop can probably be done in 1 line
                for j in ind:
                  if i == artist: #see if point in list matches current point
                    i.remove()  #remove the point
                        
        fig.canvas.draw()
#_______________________________________________________________________
def onkeyclick (event):   #def for regon selection and grafth manigment
    global NrCut,valuesFromPlot,NN,aa1,aa2
    global bb1,bb2,lineCut,Nline
    if event.xdata==None:
        print ('Did not read button, try moving the mouse first')
        return
        #else if event.key=='a' or event.key=='A':
    if event.key=='0':
        x,y = -10000,1
        aa1=x
        print "Selected from: ", aa1, "\tX,Y= ",x,y
    if event.key=='O':
        x,y = +10000,1
        aa2=x
        print "Selected from: ", aa1, "\tX,Y= ",x,y
    if event.key=='o':
        x,y = +10000,1
        aa2=x
        print "Selected from: ", aa1, "\tX,Y= ",x,y
    if event.key=='1':
        x,y = event.xdata,event.ydata
        aa1=x
        print "Selected from: ", aa1, "\tX,Y= ",x,y
    if event.key=='2':
        x,y = event.xdata,event.ydata
        aa2=x
        print "Selected   to: ", aa2, "\tX,Y= ",x,y
    if event.key=='3':
        x,y = event.xdata,event.ydata
        bb1=x
        print "Spectral Line selected from: ", bb1, "\tX,Y= ",x,y
    if event.key=='4':
        x,y = event.xdata,event.ydata
        bb2=x
        print "Spectral Line selected   to: ", bb2, "\tX,Y= ",x,y
    if event.key=='8':
        if bb1 > bb2:
            bb,bb2 = bb2,bb1
        lineCut.append((bb1,bb2))
        Nline +=1
        print "Selection LINE SAVED",Nline,lineCut[Nline-1]
        print "N!B! This version is limited to only ONE LINE"
    if event.key=='w':
        print "Width = ",bb2 - bb1, "km/sec"
    if event.key=='9':
        if aa1 > aa2:
            aa1,aa2 = aa2,aa1
        valuesFromPlot.append((aa1,aa2))
        NN +=1
        print "Selection SAVED",NN,valuesFromPlot[NN-1]
    if event.key=='5' or event.key=='q' or event.key=='Q':
        plt.close('all')
    return
#_________________________________________________________

def dummyLineRead (fileUnit,nr):
   for i in range(nr):
       line=fileUnit.readline()
   return

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#------------------------------------------------------------------------------------
from scipy import signal

import sys, os
import re
if len(sys.argv[:]) > 1:
    filename = sys.argv[1]
else:
    filename='dataFiles/m82_n19.dat' # Filename in format <ExperimentCode>_<Scan_Number>.dat

print
print "Data file: \t",filename

baseName = filename[:-8]
start='m*_n'
end='.dat'
scanNrStr = re.search('%s(.*)%s' % (start, end), filename).group(1)
print "Scan Nr: "+scanNrStr
irSchedPrms = baseName+"sch.dat"
#irLOGPrms = "prettyLogs/" + baseName+"irsch.dat" #open log reader made file
irLOGPrms = "prettyLogs/" + "m82irsch.dat"
#------------------------------------------------------------------------------------
print "Experment LOG file: \t",irLOGPrms

logpar = open(irLOGPrms,'r')

line1 = logpar.readline() # Start;Header;
while "End;Header" not in line1:
    line1 = logpar.readline()
    Header = line1.split(';')
    vards = Header[0]
    if vards == "Station":
        stationStr = Header[1]
    
print "Station: ",stationStr
print

while "Scan;"+scanNrStr not in line1:
  line1 = logpar.readline()
  if line1 is None:
    print "Data not found ERROR"
    break
# reading info from **irsch.dat file
nHeader=1
vards="Scan"
while vards != "End":
    line1 = logpar.readline()
    Header = line1.split(';')
    vards = Header[0]
    if vards == "Date":
        dateStr = Header[1]
    elif vards == "TimeStart":
        laiks = Header[1]
    elif vards == "TimeStop":
        TimeStop = Header[1]    
    elif vards == "RA":
        RaStr = Header[1]
    elif vards == "DEC":
        DecStr = Header[1]
    elif vards == "Source":
        Source = Header[1]
        Source = Source.split(",")
        Source = Source[0]
    elif vards == "FreqStart":
        FreqStart = float(Header[1])
    elif vards == "FreqStop":
        Freq2 = Header[1]
    elif vards == "FreqBBC1":
        FreqBBC1 = Header[1]
    elif vards == "Systemtemperature1":
        Tsys_1u = float(Header[2])
    elif vards == "Systemtemperature2":
        Tsys_9u = float(Header[2])
    elif vards == "AzimuthStart":
        AZ0 = Header[1]
    elif vards == "AzimuthStop":
        AZ9 = Header[1]
    elif vards == "VelocityLSR":
        VelLSR = float(Header[1])
    nHeader +=1
    
print
print "Source: ",Source
print "RA,Dec: ",RaStr+", "+DecStr
print("FreqStart", FreqStart )  
print("Tsys_1u", Tsys_1u )
print("Tsys_9u", Tsys_9u )
#------------------------------------------------------------------------------------
global bb1,bb2,lineCut,Nline
bb1 = bb2 = Nline= 0
lineCut = []

global aa1,aa2
global NrCut,valuesFromPlot
aa1=0
aa2=0
NN=0
valuesFromPlot = []
#------------------------------------------------------------------------------------
# Read Data
f = open(filename, 'r')

print "DateTime =", dateStr

import calendar, datetime
datums = dateStr[:2]
menesis = list(calendar.month_abbr).index(dateStr[3:7].strip())
gads = dateStr[-4:]
dateStr = str(gads) + " " + str(menesis) + " " + str(datums) 
try:
    timeStr = laiks.replace(":", " ")
    #print dateStr, timeStr #, timeH, timeM, timeS
except:
    print
    print "Warning: Can not get Time from data file!" 
    print
#caling dopsetpy - calculating reciver LSR velocity
dopsetPar= dateStr + " " + timeStr + " " + RaStr + " " +DecStr
print "dopsetPar "  + dopsetPar 
os.system("./dopsetpy_v1.5 "+dopsetPar)
with open('lsrShift.dat') as openfileobject:
    for line in openfileobject:
        Header = line.split(';')
        vards = Header[0]
        if vards == "Date":
            dateStr = Header[1]
        elif vards == "Time":
            laiks = Header[1]
        elif vards == "RA":
            RaStr = Header[1]
        elif vards == "DEC":
            DecStr = Header[1]
        elif vards == "Source":
            Source = Header[1]
        elif vards == "LSRshift":
            lsrShift = Header[1]
        elif vards == "MJD":
            mjd = Header[1]
            print "MJD: \t",mjd
        elif vards == "Vobs":
            Vobs = Header[1]
            print "Vobs: \t",Vobs
        elif vards == "AtFreq":
            AtFreq = Header[1]
            print "At Freq: \t",AtFreq
        elif vards == "FreqShift":
            FreqShift = Header[1]
            print "FreqShift: \t",FreqShift
        elif vards == "VelTotal":
            VelTotal = float(Header[1])
            print "VelTotal: \t",VelTotal
        nHeader +=1


Vobs = float(Vobs)
lsrCorr = float(lsrShift)*1.e6 # for MHz

print 
frequency = []
velocity = []
amplitude = []
amplit2 = []	# 2nd polarization

k=0
next(f) # skiping first line
for line in f.xreadlines():
    line = ' '.join(line.split())
    var = line.split(' ')
    var[1] = float(var[1])
    var[0] = float(var[0])
     
    frequency.append(var[0])
    velocity.append(freq2vel((var[0]+FreqStart)*10**6,VelTotal))# freqvency to velocity transformation
    
    if float(var[1]) > 10.: ###################### clean
    	#print var[1] 
    	var[1] = 0.
    if float(var[2]) > 10.: ###################### clean
    	#print var[2] 
    	var[2] = 0.
    amplitude.append(float(var[1]))
    amplit2.append(float(var[2]))
    k +=1

f.close()

#------------------------------------------------------------------------------------
pointCountOriginal = k
#print 4*f0/c

#print frequency[k-1], amplitude[k-1]

print "Number of datapoints: \t\t", k


#------------------------------------------------------------------------------------
import ast

cutFolder = "cuts/"
cutFile = filename
cutFile +=".cut"
cutFile = cutFolder+cutFile
cutFileIr = os.path.isfile(cutFile) # True if file exists
cutFileAll = cutFolder + filename.split("_")[0] + ".cutall"
if os.path.isfile(cutFileAll): # True if file exists
    # if defined cut for the source the same for all scans
    # N!B! The cut boundaries may change with time !!!
    # Delete or raname the file, eg. to .off, if you do not want to use it!
    cutFile = cutFileAll
    cutFileIr = True
if cutFileIr:
    fcut = open(cutFile)
    NN=int(fcut.readline().strip())
    cutline = fcut.readline()
#    with open(cutFile) as openfileobject:
#        for line in openfileobject:
#            cutline = line
    valuesFromPlot = ast.literal_eval(cutline)
    NN = len(valuesFromPlot)

import matplotlib; from matplotlib import pyplot as plt
import numpy as np

scale1u=calibration1(stationStr,Tsys_1u)
scale9u=calibration2(stationStr,Tsys_9u)
print 'scale1u', scale1u
print 'scale9u', scale9u
a = amplitude
a0 = np.average(a)
print "Level offset draft: \t\t",a0 
amplitude0 = (amplitude-a0) * scale1u 
amplit20 = (amplit2-a0) * scale9u

print min(velocity),max(velocity)
#_______________________________________________________________________
#Taking off Chebeshef filter bondery regons
velocity=velocity[900:3200]
amplitude0=amplitude0[900:3200]
amplit20=amplit20[900:3200]
#------------------------------------------------------------------------------------------------
# Select bad points!
fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])
ax = plt.axes()
#plt.plot(velocity,amplitude0,'go')
plt.plot(velocity,amplitude0,'bo',picker=5)
#plt.plot(velocity,amplit20,'ro')

plt.xlabel ('Velocity (km sec$^{-1}$)')
plt.ylabel ('Flux density (Jy)')
plt.title (filename[:-4]+"  "+Source)
plt.title ("Select bad data points from amplitude0")

amticks=1. # Multiplicatin factor for Y max range in plot to fit all graph
asmax = max(amplitude0)
#print asmax,asmax//amticks*amticks

am0 = min(amplitude0)
am9 = (asmax//amticks+1)*amticks # Y range
plt.ylim([am0,am9]) # Y range
#use point selec moule
browser = MyClick()
cid = fig.canvas.mpl_connect('pick_event', browser.pick)
cid = fig.canvas.mpl_connect('key_press_event', browser.press)

plt.show()
#----------------------------------------------------------------------------------
#bad point remuve from data arrays
selectedPoints = browser.pk
print selectedPoints
am11=[]
for i in selectedPoints:
  print i
  am11.append(i)
amplitude0=np.delete(amplitude0, am11)
velocity=np.delete(velocity, am11)
amplit20=np.delete(amplit20, am11)
#------------------------------------------------------------------------------------------------
# Select bad points!
fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])
ax = plt.axes()
#plt.plot(velocity,amplitude0,'go')
plt.plot(velocity,amplit20,'bo',picker=5)
#plt.plot(velocity,amplit20,'ro')

plt.xlabel ('Velocity (km sec$^{-1}$)')
plt.ylabel ('Flux density (Jy)')
plt.title (filename[:-4]+"  "+Source)
plt.title ("Select bad data points from amplit20")

amticks=1. # Multiplicatin factor for Y max range in plot to fit all graph
asmax = max(amplit20)
#print asmax,asmax//amticks*amticks

am0 = min(amplit20)
am9 = (asmax//amticks+1)*amticks # Y range
plt.ylim([am0,am9]) # Y range
#use point selec moule
browser = MyClick()
cid = fig.canvas.mpl_connect('pick_event', browser.pick)
cid = fig.canvas.mpl_connect('key_press_event', browser.press)

plt.show()
#----------------------------------------------------------------------------------    
#bad point remuve from data arrays
selectedPoints = browser.pk
print selectedPoints
am00=[]
for i in selectedPoints:
  print i
  am00.append(i)
  
print am00
amplitude0=np.delete(amplitude0, am00)
velocity=np.delete(velocity, am00)
amplit20=np.delete(amplit20, am00)  
#------------------------------------------------------------------------------------------------
# Produce initial plot raw data with noise with level offset
fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])

plt.plot(velocity,amplitude0,'go')
plt.plot(velocity,amplitude0,'r-')
plt.plot(velocity,amplit20,'ro')

plt.xlabel ('Velocity (km sec$^{-1}$)')
plt.ylabel ('Flux density (Jy)')
plt.title (filename[:-4]+"  "+Source)

print "MAX amplitude Jy draft: \t",max(amplitude0)
print
if cutFileIr:
    print "Bgr cut selections read from file: ",cutFile
    #print valuesFromPlot
else:
    print "Select the line velocity regions on plot for background calculation!"
    print "Use keys: 1,2 - select the region; 9 - save the region; 5 - exit"


amticks=1. # Multiplicatin factor for Y max range in plot to fit all graph
asmax = max(amplitude0)
#print asmax,asmax//amticks*amticks

am0 = min(amplitude0)
am9 = (asmax//amticks+1)*amticks # Y range
plt.ylim([am0,am9]) # Y range

cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
plt.show()
#----------------------------------------------------------------------------------
print "Number of Selected Cut Regions:",NN
valuesFromPlot.sort()
#print valuesFromPlot

if NN > 0 and not cutFileIr: # Do not rewrite if cut or cutall file exists
    if  os.path.isdir(cutFolder):  
        fout = open(cutFile,'w')
        print >> fout , NN
        print >> fout , valuesFromPlot
        fout.close()
        print "Cut Data saved in: \t\t\t",cutFile
    else:
        print
        print "Warning: Create the folder: ",cutFolder
        print "Warning: The Cut Data are NOT saved!"
    print

velBgr=[]
amplBgr=[]
amplBgr2 = []
nBgr = 0
velTmp = velocity[:]
amplTmp = amplitude0[:]
amplTmp2 = amplit20[:]

#print "Size= ",len(velTmp)
pointCount = len(velocity)

# Array velocity is sorted at this point!
for i in range (NN):
    #print "Start ",i
    nBgr = 0
    del velBgr[:]
    amplBgr = []
    del amplBgr[:]
    del amplBgr2[:]
    velBgr = []
    for j in range (pointCount):
    	#print "ij",i,j,nBgr
        #print valuesFromPlot[i][0],valuesFromPlot[i][1]
        #print velTmp[j]
        if velTmp[j] > valuesFromPlot[i][0] and velTmp[j] < valuesFromPlot[i][1]:
            continue
        velBgr.append(velTmp[j])
        amplBgr.append(amplTmp[j])
        amplBgr2.append(amplTmp2[j])
        nBgr +=1
    velTmp=velBgr[:]
    amplTmp=amplBgr[:]
    amplTmp2=amplBgr2[:]
    pointCount = nBgr
    #print i,nBgr

if nBgr == 0:
    velBgr = velTmp # Assume no signal all is background :) 
    amplBgr = amplTmp # Assume no signal all is background :) 
    amplBgr2 = amplTmp2 # Assume no signal all is background :) 

fitBgr =np.poly1d(np.polyfit(velBgr, amplBgr,11))
fitBgr2 =np.poly1d(np.polyfit(velBgr, amplBgr2,11))

velBgrFit = np.arange (min(velBgr),max(velBgr),1)
amplBgrFit = fitBgr(velBgrFit)
amplBgrFit2 = fitBgr2(velBgrFit)



import pylab

fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])
plt.plot(velBgr,amplBgr,'go')
plt.plot(velBgrFit,amplBgrFit,'r-')
plt.plot(velBgr,amplBgr2,'co')
plt.plot(velBgrFit,amplBgrFit2,'r-')
plt.xlabel ('Velocity (km sec$^{-1}$)')
plt.ylabel ('Flux density (Jy)')
plt.title ("Background_"+filename[:-4])
cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
figFolder="figs/"
if os.path.isdir(figFolder):  
    pylab.savefig(figFolder+filename[0:-4]+'_bg.jpg')
x1 = 0.05 * (max(velBgr) - min(velBgr)) + min(velBgr)
y1 = 0.05 * (max(amplBgr) - min(amplBgr))*(-1) + max(amplBgr)
#x1 = -30.
#y1 = 1.
#print x1,y1
#plt.text(x1, y1, r'$3\sigma=$'+str(sigma3)+' Jy')
plt.show()
# fixing data to 0 level
amplCalibrated = amplitude0-fitBgr(velocity)
amplCalibrated2 = amplit20-fitBgr2(velocity)

#---------------------------------------------------------------
# Calculate Background Noise
sigma_1u=round(3*np.std(amplCalibrated[0:500]))
sigma_9u=round(3*np.std(amplCalibrated2[0:500]))
print("3simgma L polarization: ",sigma_1u)
print("3simgma R polarization: ",sigma_9u)
#---------------------------------------------------------------

# Calculate line parameters
#------------------------------------------------------
lineCut.sort()
velLine = []
amplLine= []

#print "Size= ",len(velTmp)
hwx = hwy = []
pointCount = pointCountOriginal
if Nline > 0:
    print "Line parameters:"
    for i in range (Nline):
        nl = 0
        del velLine[:]
        velLine = []
        del amplBgr[:]
        amplLine = []
        for j in range (pointCount):
            if velocity[j] > lineCut[i][0] and velocity[j] < lineCut[i][1]:
                velLine.append(velocity[j])
                amplLine.append(amplCalibrated[j])
                nl +=1
        hwp = 0.4 # Make it shorter 1=100%, ie. -hwp/2 & -hwp/2
        hw0 = (max(velLine) - min(velLine)) *(1-hwp/2)
        hwx = [min(velLine)+hw0,max(velLine)-hw0] 
        hwy = [max(amplLine)/2,max(amplLine)/2] 
        print i,lineCut[i][0],lineCut[i][1],"\t Line max = ",max(amplLine),max(amplLine)/2
    
#------------------------------------------------------

fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])

plt.plot(velocity,amplCalibrated,'k-')
plt.plot(velocity,amplCalibrated2,'k:')
plt.xlabel ('Velocity (km sec$^{-1}$)', fontsize=20)
plt.ylabel ('Flux density (Jy)', fontsize=20)
plt.title (filename[:-4]+"  "+Source)

if Nline > 0:
    plt.plot(hwx,hwy,'b-')

amticks=1. # Multiplicatin factor for Y max range in plot to fit all graph
#amticks=350.
asmax = max(amplCalibrated)
am0 = min(amplCalibrated)
am9 = (asmax//amticks+1)*amticks # Y range
plt.ylim([am0,am9]) # Y range
#plt.xlim([55,65])
cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
print
print "MAX amplitude Jy ",asmax

if os.path.isdir(figFolder):  
    pylab.savefig(figFolder+filename[0:-4]+'.eps')
    pylab.savefig(figFolder+filename[0:-4]+'.jpg')
else:
    print "Please create folder to save plot images: ",figFolder
    print "Warning: Plot images are NOT saved!"
plt.show()


#fileOut1 = filename[0:-4]
fileOut1 = filename
fileOut1 +=".out"
print "Calibrated data out: \t\t\t",fileOut1
fout = open(fileOut1,'w')
for i in range (len(velocity)):
    #print >> fout , velocity[i],(amplitude0[i]-fitBgr(velocity[i]))
    poll = amplitude0[i]-fitBgr(velocity[i])
    polk = amplit20[i]-fitBgr2(velocity[i]) # TODO Change Bgr2
    print >> fout , velocity[i],poll,polk,(poll+polk)/2
fout.close()
print


#------------------------------------------------------
# Plot R and L polarization graph
fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])

#
velrange = 20. # half-range
XX0,YY1,YY2 = rangex(velocity,amplCalibrated,amplCalibrated2,np.mean(velocity),velrange)

velocity = XX0
amplCalibrated = YY1
amplCalibrated2 = YY2

plt.plot(velocity,amplCalibrated,'k:')
plt.plot(velocity,amplCalibrated2,'k-')
plt.xlabel ('Velocity (km sec$^{-1}$)', fontsize=20)
plt.ylabel ('Flux density (Jy)', fontsize=20 )
plt.title (filename[:-4]+"  "+Source)

if Nline > 0:
    plt.plot(hwx,hwy,'b-')

amticks=5. # Multiplicatin factor for Y max range in plot to fit all graph
#amticks=350.
asmax1 = max(amplCalibrated)
asmax2 = max(amplCalibrated2)
asmax = max(asmax1,asmax2)
am0 = min(amplCalibrated)
am9 = (asmax//amticks+1)*amticks # Y range
plt.ylim([am0,am9]) # Y range
#plt.xlim([55,65])
cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
print
print "MAX amplitude R polarization Jy ",asmax1
print "MAX amplitude L polarization Jy ",asmax2

if os.path.isdir(figFolder):  
    pylab.savefig(figFolder+filename[0:-4]+'.eps')
    pylab.savefig(figFolder+filename[0:-4]+'.jpg')
else:
    print "Please create folder to save plot images: ",figFolder
    print "Warning: Plot images are NOT saved!"
plt.show()

#------------------------------------------------------
# Plot average <R,L> polarization graph
fig = plt.figure()
fig.add_subplot(111) #,xticks=[],yticks=[])

amplAverage = (np.array(amplCalibrated) + np.array(amplCalibrated2))/2 # <R,L>

#Fiting from Astro.py-----------------------------------------------------------
import astropy
from astropy.convolution import Gaussian1DKernel, convolve, Trapezoid1DKernel, MexicanHat1DKernel
from astropy.modeling import models, fitting
from astropy.modeling.models import MexicanHat1D, Gaussian1D

# Create kernel
#trap = Trapezoid1DKernel(4, slope=1, mode='oversample',factor=2) // works good
#trap = Trapezoid1DKernel(4, slope=1, mode='oversample',factor=2) #works good
gaus = Gaussian1DKernel(stddev=2, x_size=19, mode='center', factor=100) #givs best convol. 
#mex = MexicanHat1DKernel(0.7, x_size=1, mode='oversample', factor=500000)
#han = np.hanning(12) 
#han=han/han.sum()   #normalaizing
#hann = np.convolve(han, amplAverage, mode='SAME')
# Convolve data
z = convolve(amplAverage, gaus)
#-------------------------------------------------------------------------------

plt.plot(velocity,amplAverage,'k:')
plt.plot(velocity,z,'r-') #ploting convolved data
#plt.plot(velocity,hann,'g-')
plt.xlabel ('Velocity (km sec$^{-1}$)', fontsize=20)
plt.ylabel ('Flux density (Jy)', fontsize=20)
plt.title (filename[:-4]+"  "+Source+" <R,L> MJD="+mjd)

if Nline > 0:
    plt.plot(hwx,hwy,'b-')

asmax = max(amplAverage)
plt.ylim([am0,am9]) # Y range
cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
print "MAX amplitude <R,L> Jy ",asmax

if os.path.isdir(figFolder):  
    pylab.savefig(figFolder+filename[0:-4]+'av.eps')
    pylab.savefig(figFolder+filename[0:-4]+'av.jpg')
else:
    print "Please create folder to save plot images: ",figFolder
    print "Warning: Plot images are NOT saved!"
plt.show()
###################################################
# Maxim found script imports
import numpy
import peakutils
from peakutils.plot import plot as pplot
##############################################
# Plot Astropy fit vith local maxim POINTS
zz=np.array(z)
vv=np.array(velocity)
#smart thres fiding!!!
smartthres=0.9*sigma_1u/asmax
print "smartthres ", smartthres
print("% smaler as max vill be serched",smartthres)
indexes = peakutils.indexes(zz, thres=smartthres, min_dist=3)
print(indexes)
print(vv[indexes], zz[indexes])
fig = plt.figure()
ax = fig.add_subplot(111)
for xy in zip(vv[indexes], zz[indexes]):                        # <--
    ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data') # <--
fig.add_subplot(111) #,xticks=[],yticks=[])
plt.plot(velocity,z,'k-')
plt.plot(vv[indexes], zz[indexes], 'ro')
plt.xlabel ('Velocity (km sec$^{-1}$)', fontsize=20)
plt.ylabel ('Flux density (Jy)', fontsize=20)
plt.title (filename[:-4]+"  "+Source+"Fit... MJD="+mjd)

if Nline > 0:
    plt.plot(hwx,hwy,'b-')

asmax = max(amplAverage)

plt.ylim([am0,am9]) # Y range
cid = fig.canvas.mpl_connect('key_press_event',onkeyclick)
print "MAX amplitude <R,L> Jy ",asmax

if os.path.isdir(figFolder):  
    pylab.savefig(figFolder+filename[0:-4]+'fit.eps')
    pylab.savefig(figFolder+filename[0:-4]+'fit.jpg')
else:
    print "Please create folder to save plot images: ",figFolder
    print "Warning: Plot images are NOT saved!"
plt.show()
