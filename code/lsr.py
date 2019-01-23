from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u
from astropy.coordinates import Angle, HeliocentricTrueEcliptic
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, BarycentricTrueEcliptic, Angle
import astropy.coordinates as coord
import numpy as np
import datetime

import os

#cepa: 225617.90, 620149.7, 2000.0

ra = 225617.90
dec = 620149.7
sc = SkyCoord(ra='22h56m17.90s', dec='+62d01m49.7s', unit=u.deg, frame=FK5, equinox='J2000.0')
sc.transform_to(sc) 
ra=sc.ra
dec=sc.dec
ra = Angle(ra, u.radian)
dec = Angle(dec, u.radian)
dec = float(dec.to_string(unit=u.radian, decimal=True))
ra = float(ra.to_string(unit=u.radian, decimal=True))

epoch = 2000.0

year = 2019
month = 1
day = 22
timeUT = '08:30:30'

lat = '57:33:12.3' 
lat = 57.553417
altitude = 87.30
long = '-021:51:17'
long = -21.854722

v_sun = coord.Galactocentric.galcen_v_sun.to_cartesian()
print("v_sun", v_sun)
t = datetime.datetime(2019,1,22,8,30,30)
irbene = EarthLocation.from_geodetic(lat=lat*u.deg, lon=long*u.deg, height=altitude*u.m)
sc = SkyCoord(ra=ra*u.radian, dec=dec*u.radian)
barycorr = sc.radial_velocity_correction(obstime=Time(t), location=irbene)
barycorr = barycorr.to(u.km/u.s)
print ("barycorr", barycorr)
heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(t), location=irbene)  
heliocorr = heliocorr.to(u.km/u.s) 
print("heliocorr", heliocorr)

l = float(Angle(BarycentricTrueEcliptic(lon=ra*u.radian, lat=dec*u.radian).lon).to_string(unit=u.radian, decimal=True))
b = float(Angle(BarycentricTrueEcliptic(lon=ra*u.radian, lat=dec*u.radian).lat).to_string(unit=u.radian, decimal=True))
print("l", l, "b", b)


vLSR = barycorr.value + 9 * np.cos(l) * np.cos(b) + 12 * np.sin(l) * np.cos(b) + 7 * np.sin(b)
print("vLSR", vLSR)

vGSR = vLSR + 220 * np.sin(l)*np.cos(b)
print("vGSR", vGSR)

vLGSR = vGSR - 62 * np. cos(l)* np.cos(b) + 40 * np.sin(l) * np.cos(b) - 35 * np.sin(b)
print("vLGSR", vLGSR)

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









