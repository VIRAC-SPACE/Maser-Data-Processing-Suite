import numpy as np
from astropy.time import Time
import os

#cepa: 225617.90, 620149.7, 2000.0

ra = 225617.90
dec = 620149.7
epoch = 2000.0

year = 2019
month = 1
day = 22
timeUT = '08:30:30'

vsun = 11.192600995434107
ra_vsun = "18:00:00"
dec_vsun = "30:00:00"
epoch_vsun = 1900
vobs = 0.11203614556028638

lat = '57:33:12.3' 
lat = 57.553417
altitude = 87.30
long = '-021:51:17'
long = -21.854722

dlat = -(11. * 60. + 32.743000) * np.sin (2*lat) + 1.163300 * np.sin (4*lat) - 0.002600 * np.sin (6*lat)
r = 6378160.0 * (0.998327073 + 0.00167643800 * np.cos(2*lat) - 0.00000351 * np.cos(4*lat) + 0.000000008 * np.cos(6*dlat)) + altitude
TWOPI = np.pi * 2
times = ['1999-01-01T00:00:00.123456789', '2010-01-01T00:00:00']
t = Time(times, format='isot', scale='utc')
lmst =  7.9627442150#t.sidereal_time('apparent', 'greenwich')
v = TWOPI * (r / 1000.) / (23.934469591229 * 3600.)
vdiurnal = v * np.cos (lat) * np.cos (dec) * np.sin (ra-lmst)

JD = 58330.212685185186
t = (JD - 2415020) / 36525
oblq = 23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503))
omega = 259.183275-t*(1934.142008+t*(0.002078+t*0.000002))
mlong = 270.434164+t*(481267.88315+t*(-0.001133+t*0.0000019))-omega
lperi = 334.329556+t*(4069.034029-t*(0.010325+t*0.000012))-omega
em = 0.054900489
inclin = 5.1453964
manom = mlong - lperi
tanom = manom + (2 * em - 0.25 * em**3) * np.sin (manom) + 1.25 * em**2 * np.sin (2 * manom) + 13/12 * em**3 * np.sin (3 * manom)
tlong = tanom + lperi

vmoon = (TWOPI * 384403.12040) / (27.321661 * 86400) /np.sqrt (1. - em**2)
v = vmoon * np.cos (long) * (np.sin (tlong-lat) - em*np.sin (lperi-lat))
vlunar = v / 81.53

t = (JD - 2415020) / 36525
manom = 358.47583+t*(35999.04975-t*(0.000150+t*0.000003))
oblq = 23.452294-t*(0.0130125+t*(0.00000164-t*0.000000503))
lperi = 101.22083+t*(1.7191733+t*(0.000453+t*0.000003))
eccen = 0.01675104-t*(0.00004180+t*0.000000126)
tanom = manom + (2 * eccen - 0.25 * eccen**3) * np.sin (manom) +1.25 * eccen**2 * np.sin (2 * manom) + 13./12. * eccen**3 * np.sin (3 * manom)
v = ((TWOPI * 149598500.) / (365.2564 * 86400.)) /np.sqrt (1. - eccen**2)
slong = lperi + tanom + 180
vannual = v * np.cos(em) * (np.sin(slong-lperi) - eccen*np.sin(lperi-lperi))
print(vannual)

#cepa: 225617.90, 620149.7, 2000.0
RaStr = "22 56 17.90"
DecStr = "62 01 49.7"
dopsetPar = "2019 01 22 08 30 30 22 56 17.90 62 01 49.7"
print(dopsetPar)
os.system("code/dopsetpy_v1.5 " + dopsetPar)

# dopsetpy parametru nolasisana
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
            print ("MJD: \t", mjd)
        elif vards == "Vobs":
            Vobs = Header[1]
            print ("Vobs: \t", Vobs)
        elif vards == "AtFreq":
            AtFreq = Header[1]
            print ("At Freq: \t", AtFreq)
        elif vards == "FreqShift":
            FreqShift = Header[1]
            print ("FreqShift: \t", FreqShift)
        elif vards == "VelTotal":
            VelTotal = float(Header[1])
            print ("VelTotal: \t", VelTotal)
            #Header +=1









