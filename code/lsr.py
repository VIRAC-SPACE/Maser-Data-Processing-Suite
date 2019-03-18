import os
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, FK5, EarthLocation, AltAz, Angle
from astropy.time import Time
import math
import datetime
from jplephem.spk import SPK

import time
from novas.tests.test_example import height

os.environ['TZ'] = 'UTC'
time.tzset()

kernel = SPK.open('/home/janis/Downloads/de435.bsp')

def dopsety(dopsetPar):
    
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
            elif vards == "Vobs":
                Vobs = Header[1]
            elif vards == "AtFreq":
                AtFreq = Header[1]
            elif vards == "FreqShift":
                FreqShift = Header[1]
            elif vards == 'VSun':
                Vsun = Header[1]
            elif vards == 'vEarth':
                VEarth = Header[1]
            elif vards == "VelTotal":
                VelTotal = float(Header[1])
                #Header +=1
    print("Vobs f", Vobs)
    #print ("Vearth f", float (VEarth) )
    print ("Vobs f", float(Vobs))
    #print ("Vsun f", float(Vsun))
    #return float(VEarth) 
    #return(float(Vobs))
    #return float(Vsun)          
    return VelTotal 

@u.quantity_input(canon_velocity=(u.meter / u.second))
def v_sun(source, apex = "18h03m50.29s +30d00m16.8s", canon_velocity = 20.0 * 1000 * u.meter / u.second):
    if isinstance(source.obstime.value, str):
        equinox = source.obstime
    else:
        equinox = source.obstime[0]
        
    #SkyCoord.from_name("Cep A")      
    c = SkyCoord(apex, equinox="J2000", frame=FK5)
    c = c.transform_to(FK5(equinox=equinox, representation='cartesian'))
    c = c.cartesian
    source = source.cartesian
    #print("vsun", canon_velocity*np.dot(c.xyz,source.xyz)/1000)
    return canon_velocity*np.dot(c.xyz,source.xyz)

def v_earth(source):
    time = source.obstime.jd
    position, velocity = kernel[0,3].compute_and_differentiate(time)
    position2, velocity2 = kernel[3,399].compute_and_differentiate(time)
    velocity = velocity - velocity2
    source = source.cartesian
    #print("vearth ", np.mean(np.dot(((velocity * 1000 * u.meter / u.day).to(u.meter / u.second)).T, source.xyz)).value /1000)
    return np.dot(((velocity * 1000 * u.meter / u.day).to(u.meter / u.second)).T, source.xyz) 

def convertDatetimeObjectToJD(time):
    time=time.isoformat()
    t=Time(time, format='isot')
    return t.jd

def vbos(dec, ra, time):
    def dmod(a, b):
        return np.double(a) % np.double(b)
    
    dj = convertDatetimeObjectToJD(time)
    print ("djp", dj)
    ro=3430.18601252541
    vhor=2*np.pi*ro/(24*3600)*1.002737909350795
    print("vhorp", vhor)
    angleDEC = SkyCoord(ra='22h56m17.90s', dec='+62d01m49.7s').dec.radian
    cdec=np.cos(angleDEC)
    #cdec=np.cos(dec)
    print("cdecp", cdec)
    albypi=0.121415123729
    t=dj-2451545 
    print("t", t)
    print("tp", t)
    #t=2
    #albypi=3 
    print("t", t)
    aLMST=np.pi*math.fmod(np.double(5.55811454652)+math.fmod(t+t,np.double(2))+ t*(np.double(0.547581870159)-2+t*(np.double(1.61549)-15-t*np.double(1.473)-24))+albypi,np.double(2))
    #aLMST = np.rad2deg(aLMST)
    print("aLMSTP",aLMST)
    #location = EarthLocation(x=3183.661, y=1276.902, z=5359.291).value
    location = EarthLocation(x=3183661 * u.m, y=1276902 * u.m, z=5359291 * u.m)
    print("location", location)
    #('-21d51m17s', '57d33m12.3s','87.30m')
    sdTIme = Time('2019-02-18 08:30:30', scale='utc',location=location)
    aLMST2 = np.deg2rad(sdTIme.sidereal_time('apparent', 'greenwich')).value
    print("aLMST2p", aLMST2 , "sdTIME", sdTIme)
    print("math.fmod(t+t,2.0)", math.fmod(t+t,2.0))
    
    aLMSTf = 5.1923676562989565
    ravhor=aLMST+np.pi/2
    print("ravhorp", ravhor)
    print("Angle", SkyCoord(ra='22h56m17.90s', dec='+62d01m49.7s').ra.radian)
    print("ra python", np.deg2rad(ra))
    angleRA = SkyCoord(ra='22h56m17.90s', dec='+62d01m49.7s').ra.radian
    print("np.cos(np.deg2rad(ra)-ravhor)p", np.cos(angleRA-ravhor))
    vobs=vhor*cdec*np.cos(angleRA-ravhor)
    #vobs=vhor*cdec*np.cos(ra-ravhor)
    print("vobsp", vobs) 
    return (vobs)

year = 2019
month = 2
day = 18

def v_lsr(source, time):
    ra = 22 + 56/60 + 17.90/3600
    dec = 62 + 1/60 + 49.7/3600
    rar0  = (22 + 56/60.0 + 17.90/3600.0)/12*np.pi
    decr0 = (np.abs(62)+np.abs(1)/60.0+np.abs(49.7)/3600.0)/180*np.pi
    v =  np.mean(v_sun(source)).value / 1000  + vbos(decr0, rar0, time)  + np.mean(v_earth(source)).value /1000
    #v = np.mean(v_earth(source)).value /1000
    #v =  np.mean(v_sun(source)).value / 1000
    #v =  np.mean(v_earth(source)).value /1000
    #v = vbos(dec, ra, time)
    return v

lat = '57:33:12.3' 
lat = 57.553417
altitude = 87.30
long = '-021:51:17'
long = -21.854722

timeUT = '08:30:30'

#cepa: 225617.90, 620149.7, 2000.0
ra = 225617.90
dec = 620149.7
epoch = 2000.0

def lsr(ra, dec, date, t):
    start_time = Time(date, format='isot', scale='utc')
    time_range = np.linspace(0, 10000, 10000)*u.second
    times = start_time + time_range
    source = SkyCoord(ra=ra, dec=dec, frame=FK5, equinox='J2000.0', obstime=times)
    source.transform_to(source)
    V_lsr= v_lsr(source, t)
    return V_lsr

t = datetime.datetime(year, month, day, 8, 30, 30)
time=t.isoformat()
date = Time(time, format='isot', scale='utc')


vTotal = lsr('22h56m17.90s', '+62d01m49.7s', date, t)
print("vTotal p", vTotal)

dopsetPar = "2019 02 18 08 30 30 22 56 17.90 62 01 49.7"
print("vTotal f", dopsety(dopsetPar))


def convertDatetimeObjectToMJD(time):
    time=time.isoformat()
    t=Time(time, format='isot')
    return t.mjd

'''
for i in range (0, 512):
    time = datetime.datetime.today() - datetime.timedelta(weeks = i)
    timeStr = datetime.datetime.strftime(time , '%Y %m %d %H %M %S')
    dopsetPar = timeStr + " 22 56 17.90 62 01 49.7"
    lsrDatePar = datetime.datetime.strftime(time , '%Y-%m-%dT%H:%M:%S')
    t=time.isoformat()
    date = Time(t, format='isot', scale='utc')
    #print(type(lsrDatePar))
    #lsr('22h56m17.90s', '+62d01m49.7s', lsrDatePar, time)
    #print(convertDatetimeObjectToMJD(time), dopsety(dopsetPar), lsr('22h56m17.90s', '+62d01m49.7s', lsrDatePar), np.abs(dopsety(dopsetPar)) - np.abs(lsr('22h56m17.90s', '+62d01m49.7s', lsrDatePar)))
    print(convertDatetimeObjectToMJD(time), dopsety(dopsetPar), lsr('22h56m17.90s', '+62d01m49.7s', t, time), np.abs(dopsety(dopsetPar)) - np.abs(lsr('22h56m17.90s', '+62d01m49.7s', t, time)))
'''

