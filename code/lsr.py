from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5
import astropy.units as u
from astropy.coordinates import Angle, HeliocentricTrueEcliptic
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, BarycentricTrueEcliptic, Angle, Galactic, ICRS, LSR
import astropy.coordinates as coord
import numpy as np
import datetime

import jplephem
from jplephem.spk import SPK
kernel = SPK.open('/home/janis/Downloads/de430.bsp')

import astropy
from astropy.coordinates import SkyCoord, ICRS, FK5, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

import math
import numpy as np

import os

@u.quantity_input(canon_velocity=(u.meter / u.second))
def v_sun(source, apex = "18h03m50.29s +30d00m16.8s", 
                  canon_velocity = 20.0 * 1000 * u.meter / u.second):
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
    v =  v_sun(source) + (v_earth(source) + v_obs(source, observer))
    return v


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


start_time = Time('2019-01-22T08:30:30', format='isot', scale='utc')
time_range = np.linspace(0, 100, 100)*u.second
times = start_time + time_range
acre_road = EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=altitude*u.m)
source = SkyCoord(ra=ra*u.radian, dec=dec*u.radian, unit=u.radian, frame=FK5, equinox='J2000.0', obstime=times)
#source = SkyCoord(306.51141024*u.degree, 40.84802122*u.degree, frame=FK5, obstime=times)
'''
source.transform_to(source) 
'''

V= v_lsr(source, acre_road)
print("V", V)

Vavg = np.mean(V)/1000

print("Vavg", Vavg)

'''
v_sun = coord.Galactocentric.galcen_v_sun.to_cartesian()
print("v_sun", v_sun)
t = datetime.datetime(2019,1,22,8,30,30)
irbene = EarthLocation.from_geodetic(lat=lat*u.deg, lon=long*u.deg, height=altitude*u.m)
sc = SkyCoord(ra=ra*u.radian, dec=dec*u.radian)
barycorr = sc.radial_velocity_correction('barycentric', obstime=Time(t), location=irbene).to('km/s')
barycorr = barycorr.to(u.km/u.s)
print ("barycorr", barycorr)
heliocorr = sc.radial_velocity_correction('heliocentric', obstime=Time(t), location=irbene)  
heliocorr = heliocorr.to(u.km/u.s) 
print("heliocorr", heliocorr)

icrs = ICRS(ra=ra*u.radian, dec=dec*u.radian)  
icrs = icrs.transform_to(Galactic)

icrsL = icrs.l
icrsB = icrs.b
icrsL = float(Angle(icrsL).to_string(unit=u.radian, decimal=True))
icrsB = float(Angle(icrsB).to_string(unit=u.radian, decimal=True))
print("icrsL", icrsL, "icrsB ", icrsB)
vLSRicrs = barycorr.value + 9 * np.cos(icrsL) * np.cos(icrsB) + 12 * np.sin(icrsL) * np.cos(icrsB) + 7 * np.sin(icrsB)
print("vLSRicrs", vLSRicrs)
#LSR = icrs.transform_to(LSR)  

l = float(Angle(BarycentricTrueEcliptic(lon=ra*u.radian, lat=dec*u.radian).lon).to_string(unit=u.radian, decimal=True))
b = float(Angle(BarycentricTrueEcliptic(lon=ra*u.radian, lat=dec*u.radian).lat).to_string(unit=u.radian, decimal=True))
print("l", l, "b", b)

vLSR = barycorr.value + 9 * np.cos(l) * np.cos(b) + 12 * np.sin(l) * np.cos(b) + 7 * np.sin(b)
print("vLSR", vLSR)

vGSR = vLSR + 220 * np.sin(l)*np.cos(b)
print("vGSR", vGSR)

vLGSR = vGSR - 62 * np. cos(l)* np.cos(b) + 40 * np.sin(l) * np.cos(b) - 35 * np.sin(b)
print("vLGSR", vLGSR)
'''
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
