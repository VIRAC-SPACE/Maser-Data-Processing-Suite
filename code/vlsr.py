import os
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, ICRS, FK5, EarthLocation, AltAz, Angle
from astropy.time import Time
import datetime
from jplephem.spk import SPK
import time

@u.quantity_input(canon_velocity=(u.meter / u.second))
def v_sun(source, apex="18h03m50.29s +30d00m16.8s", canon_velocity=20.0 * 1000 * u.meter / u.second):
    if isinstance(source.obstime.value, str):
        equinox = source.obstime
    else:
        equinox = source.obstime[0]

    c = SkyCoord(apex, equinox="J2000", frame=FK5)
    c = c.transform_to(FK5(equinox=equinox, representation='cartesian'))
    c = c.cartesian
    source = source.cartesian
    return canon_velocity * np.dot(c.xyz, source.xyz)

def v_earth(source):
    kernel = SPK.open('de435.bsp')
    time = source.obstime.jd
    position, velocity = kernel[0, 3].compute_and_differentiate(time)
    position2, velocity2 = kernel[3, 399].compute_and_differentiate(time)
    velocity = velocity - velocity2
    source = source.cartesian
    return np.dot(((velocity * 1000 * u.meter / u.day).to(u.meter / u.second)).T, source.xyz)

def convertDatetimeObjectToJD(time):
    time = time.isoformat()
    t = Time(time, format='isot')
    return t.jd

def vobs(dec, ra, time, stringTime):
    dj = convertDatetimeObjectToJD(time)
    x = 3183661
    y = 1276902
    ro = np.sqrt(x ** 2 + y ** 2) / 1000.0
    vhor = 2 * np.pi * ro / (24 * 3600) * 1.002737909350795
    angleDEC = SkyCoord(ra='22h56m17.90s', dec='+62d01m49.7s').dec.radian
    cdec = np.cos(angleDEC)
    x = x / 1000.0
    y = y / 1000.0
    along = np.arctan2(y, x)
    albypi = along / np.pi
    t = dj - 2451545
    location = EarthLocation(x=3183661 * u.m, y=1276902 * u.m, z=5359291 * u.)
    sdTIme = Time(stringTime, scale='utc', location=location)
    aLMST = np.deg2rad(sdTIme.sidereal_time('apparent')).value
    ravhor = aLMST + np.pi / 2
    angleRA = SkyCoord(ra='22h56m17.90s', dec='+62d01m49.7s').ra.radian
    vobs = vhor * cdec * np.cos(angleRA - ravhor)
    return vobs

def v_lsr(source, time, stringTime):
    rar0 = (22 + 56 / 60.0 + 17.90 / 3600.0) / 12 * np.pi
    decr0 = (np.abs(62) + np.abs(1) / 60.0 + np.abs(49.7) / 3600.0) / 180 * np.pi
    v = np.mean(v_sun(source)).value / 1000 + vobs(decr0, rar0, time, stringTime) + np.mean(v_earth(source)).value / 1000
    return v

def lsr(ra, dec, date, t, stringTime):
    start_time = Time(date, format='isot', scale='utc')
    time_range = np.linspace(0, 10000, 10000) * u.second
    times = start_time + time_range
    source = SkyCoord(ra=ra, dec=dec, frame=FK5, equinox='J2000.0', obstime=times)
    source.transform_to(source)
    V_lsr = v_lsr(source, t, stringTime)
    return V_lsr