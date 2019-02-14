from __future__ import print_function
from __future__ import division



import jplephem
from jplephem.spk import SPK
import astropy
from astropy.coordinates import SkyCoord, ICRS, FK5, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import get_body_barycentric, get_body
import math
import numpy as np
import os
import datetime

from barycorrpy import get_BC_vel

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
                     
def dopsety2(dopsetPar):
    
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
            elif vards == "VelTotal":
                VelTotal = float(Header[1])
                #Header +=1
                
    return VelTotal 

@u.quantity_input(canon_velocity=(u.meter / u.second))
def v_sun(source, apex = "18h03m50.29s +30d00m16.8s", canon_velocity = 20.0 * 1000 * u.meter / u.second):
    """
    Calculate the velocity of the sun relative to a point in the sky.
    
    Parameters
    ----------
    source : `astropy.coordinate.SkyCoord` object
        The point on the sky which the Earth's projected velocity should be 
        calculated for.
        
    apex : str, optional
        The point in the sky which the Sun is travelling towards. This should 
        be given as a string in the format "ra dec". The default value is
        18h03m50.29s +30d00m16.8s.
        
    canon_velocity : float
        The velocity at which the sun is travelling towards `apex`. This 
        defaults to 20 000 m/s. The value must be specified using astropy 
        units.
        
    Returns
    -------
    v_sun : `astropy.quantity`
        The velocity of the sun, relative to the source, with correct units.
    
    Notes
    -----
    The Sun is moving towards a point in the sky close to Vega at a velocity
    of around 20 km/s. In order to calculate the velocity of the Earth towards 
    any other point in the sky we need to project the velocity vector onto
    the position vector of that point in the sky, using a dot product.
    """
    if isinstance(source.obstime.value, str):
        equinox = source.obstime
    else:
        equinox = source.obstime[0]
        
    #SkyCoord.from_name("Cep A")      
    c = SkyCoord(apex, equinox="J2000", frame=FK5)
    c = c.transform_to(FK5(equinox=equinox, representation='cartesian'))
    c = c.cartesian
    source = source.cartesian
    v = canon_velocity * c.xyz
    return canon_velocity*np.dot(c.xyz,source.xyz)

def v_obs(source, observer):
    """
    Calculate the velocity of the observatory due to its rotation about the 
    Earth's barycentre in the direction of a point on the sky.
    
    Parameters
    ----------
    source : `astropy.coordinate.SkyCoord` object
        The point on the sky which the Earth's projected velocity should be 
        calculated for.
    observer : `astropy.coordinates.EarthLocation` object
        The location on the Earth of the observatory.
    
    Returns
    -------
    v_obs : `astropy.quantity` object
        The velocity of the observer relative to the source.
    """
    pos = [observer.x.value, observer.y.value, observer.z.value] * u.meter
    vel = [0, 0,  (2*np.pi)/(24.*3600.)]  / u.second # 459
    source = source.transform_to(AltAz(location=observer))
    return np.dot(np.cross(pos, vel), source.cartesian.xyz)*u.meter/u.second

def v_earth(source):
    """
    Calculate the velocity of the Earth relative to some point on the sky.
    
    Parameters
    ----------
    source : `astropy.coordinate.SkyCoord`
        The point on the sky which the Earth's projected velocity should be 
        calculated for.
        
    Returns
    -------
    v_earth : `astropy.quantity` object
        The velocity of the Earth relative to the source.
    """
    time = source.obstime.jd
    position, velocity = kernel[0,3].compute_and_differentiate(time)
    position2, velocity2 = kernel[3,399].compute_and_differentiate(time)
    velocity = velocity - velocity2
    source = source.cartesian
    return np.dot(((velocity * 1000 * u.meter / u.day).to(u.meter / u.second)).T, source.xyz) 

def v_lsr(source, observer):
    """
    Calculate the velocity of the local standard of rest for an observatory in
    the direction of a source on the sky.
    
    Parameters
    ----------
    source : `astropy.coordinate.SkyCoord` object
        The point on the sky which the Earth's projected velocity should be 
        calculated for.
    observer : `astropy.coordinates.EarthLocation` object
        The location on the Earth of the observatory.
    
    Returns
    -------
    v_lsr : `astropy.quantity` object
        The velocity of the local standard of rest relative to the source.
    """
    
    print("v_sun", np.mean(v_sun(source))/1000)
    print("v_earth", np.mean(v_earth(source))/1000)
    print("v_obs",  np.mean(v_obs(source, observer))/1000)
    v =  v_sun(source) + (v_earth(source) + v_obs(source, observer))
    return v

def lsr(ra, dec, date):
    start_time = Time(date, format='isot', scale='utc')
    time_range = np.linspace(0, 10000, 10000)*u.second
    times = start_time + time_range
    acre_road = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=altitude*u.m)
    source = SkyCoord(ra=ra, dec=dec, frame=FK5, equinox='J2000.0', obstime=times)
    source.transform_to(source) 
    V_lsr= v_lsr(source, acre_road)
    Vavg = np.mean(V_lsr)/1000
    #print("Vavg", Vavg)
    return Vavg.value
    
lat = '57:33:12.3' 
lat = 57.553417
altitude = 87.30
long = '-021:51:17'
long = -21.854722

year = 2019
month = 1
day = 22
timeUT = '08:30:30'

#cepa: 225617.90, 620149.7, 2000.0
ra = 225617.90
dec = 620149.7
epoch = 2000.0

dopsetPar = "2019 01 22 08 30 30 22 56 17.90 62 01 49.7"
print(dopsety(dopsetPar))
#lsr('22h56m17.90s', '+62d01m49.7s')

def convertDatetimeObjectToMJD(time):
    time=time.isoformat()
    t=Time(time, format='isot')
    return t.mjd

time = datetime.datetime(year, month, day, 8, 30, 30)
time=time.isoformat()
t = Time(time, format='isot', scale='utc')
a = get_BC_vel (t, lat=lat, longi=long, alt=altitude, ephemeris='de430', zmeas=0.0, ra=ra, dec=dec, leap_update=True, epoch=2451545.0)
b = a[0] /1000
print("barry", b)

lsr('22h56m17.90s', '+62d01m49.7s', t)

'''  
for i in range (0, 520):
    time = datetime.datetime.today() - datetime.timedelta(weeks = i)
    timeStr = datetime.datetime.strftime(time , '%Y %m %d %H %M %S')
    dopsetPar = timeStr + " 22 56 17.90 62 01 49.7"
    lsrDatePar = datetime.datetime.strftime(time , '%Y-%m-%dT%H:%M:%S')
    lsr('22h56m17.90s', '+62d01m49.7s', lsrDatePar)
    print(convertDatetimeObjectToMJD(time), dopsety2(dopsetPar), lsr('22h56m17.90s', '+62d01m49.7s', lsrDatePar), np.abs(dopsety2(dopsetPar)) - np.abs(lsr('22h56m17.90s', '+62d01m49.7s', lsrDatePar)))
'''

'''
#g107p3: 222126.81, 635137.14, 2000
print ("\n" + "g107p3", "\n")
dopsetPar = "2019 01 22 08 30 30 22 21 26.81 63 51 37.14"
#print(dopsety(dopsetPar))
lsr('22h21m26.81s', '+63d51m37.14s')
'''






