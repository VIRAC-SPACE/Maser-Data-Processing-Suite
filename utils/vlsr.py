"""
 Compute Local Standard of Rest
"""
import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, FK5, EarthLocation
from astropy.time import Time
from jplephem.spk import SPK


@u.quantity_input(canon_velocity=(u.meter / u.second))
def v_sun(source, apex="18h03m50.29s +30d00m16.8s",
          canon_velocity=19.954 * 1000 * u.meter / u.second):
    """

    :param source: observed source
    :param apex: apex
    :param canon_velocity: canon velocity
    :return: sun velocity
    """
    if isinstance(source.obstime.value, str):
        equinox = source.obstime
    else:
        equinox = source.obstime[0]

    sky_coord = SkyCoord(apex, equinox="J2000", frame=FK5)
    sky_coord = sky_coord.transform_to(FK5(equinox=equinox, representation='cartesian'))
    sky_coord = sky_coord.cartesian
    source = source.cartesian
    return canon_velocity * np.dot(sky_coord.xyz, source.xyz)


def v_earth(source):
    """

    :param source: observed source
    :return: Earth velocity
    """
    kernel = SPK.open('de435.bsp')
    time = source.obstime.jd
    _, velocity = kernel[0, 3].compute_and_differentiate(time)
    _, velocity2 = kernel[3, 399].compute_and_differentiate(time)
    velocity = velocity - velocity2
    source = source.cartesian
    return np.dot(((velocity * 1000 * u.meter / u.day).to(u.meter / u.second)).T, source.xyz)


def vobs(ra, dec, string_time, x, y, z):
    """

    :param ra: Right ascension
    :param dec: Declination
    :param string_time: observation start time
    :param x: x coordinate
    :param y: y coordinate
    :param z: z coordinate
    :return: observatory velocity
    """
    ro = np.sqrt(x ** 2 + y ** 2) / 1000.0
    vhor = 2 * np.pi * ro / (24 * 3600) * 1.002737909350795
    angleDEC = SkyCoord(ra=ra, dec=dec).dec.radian
    cdec = np.cos(angleDEC)
    location = EarthLocation(x=x * u.m, y=y * u.m, z=z * u.m)
    sd_time = Time(string_time, scale='utc', location=location)
    aLMST = np.deg2rad(sd_time.sidereal_time('apparent')).value
    ravhor = aLMST + np.pi / 2
    angleRA = SkyCoord(ra=ra, dec=dec).ra.radian
    vobs = vhor * cdec * np.cos(angleRA - ravhor)
    return vobs


def v_lsr(ra, dec, source, string_time, x, y, z):
    """

    :param ra: Right ascension
    :param dec: Declination
    :param source: observed source
    :param string_time: observation start time
    :param x: x coordinate
    :param y: y coordinate
    :param z: z coordinate
    :return: Local Standard of Rest
    """
    v = np.mean(v_sun(source)).value / 1000 + \
        vobs(ra, dec, string_time, x, y, z) + np.mean(v_earth(source)).value / 1000
    return v


def lsr(ra, dec, date, string_time, x, y, z):
    """

    :param ra: Right ascension
    :param dec: Declination
    :param date: date
    :param string_time: observation start time
    :param x: x coordinate
    :param y: y coordinate
    :param z: z coordinate
    :return: Local Standard of Rest
    """
    start_time = Time(date, format='isot', scale='utc')
    time_range = np.linspace(0, 1, 1) * u.second
    times = start_time + time_range
    source = SkyCoord(ra=ra, dec=dec, frame=FK5, equinox='J2000.0', obstime=times)
    source.transform_to(source)
    V_lsr = v_lsr(ra, dec, source, string_time, x, y, z)
    return V_lsr
