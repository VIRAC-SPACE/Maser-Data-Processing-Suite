"""
common used functions
"""

from functools import reduce
from numpy import trapz
import numpy as np
from astropy.time import Time
from astropy.modeling import models
from astropy.modeling.fitting import LevMarLSQFitter


class Experiment:
    """
     Experiment class
    """

    def __init__(self, **entries):
        self.flag = None
        self.modifiedJulianDays = None
        self.Iteration_number = None
        self.polarizationAVG = None
        self.polarizationU9 = None
        self.polarizationU1 = None
        self.location = None
        self.Date = None
        self.specie = None
        self.type = None
        self.__dict__.update(entries)


def file_len(file_name):
    """

    :param file_name: file name
    :return: rows of file
    """
    i = 0
    with open(file_name) as file:
        for i, _ in enumerate(file):
            pass
    return i + 1


def get_iteration_from_output_file(file):
    return int(file.split(".")[1].split("_")[-1])


def indexies(array, value):
    """

    :param array: nupmy array
    :param value: value to search
    :return: indexes in array for value
    """
    indexes = list()
    for i in range(0, len(array) - 1):
        if array[i] == value:
            indexes.append(i)
    return indexes


def find_nearest_index(array, value):
    """

    :param array: array
    :param value: value to search
    :return: nearest index of value for array
    """
    array = np.asarray(array)
    index = (np.abs(array - value)).argmin()
    return index


def correct_numpy_read_data(data):
    """

    :param data: numpy array
    :return: corrected numpy array
    """
    corrected_data = list()
    for dd_ in data:
        corrected_data.append(dd_[0])
    return np.array(corrected_data)


def convert_datetime_object_to_mjd(time):
    """

    :param time: datetime object
    :return: MJD
    """
    time = time.isoformat()
    tt_ = Time(time, format='isot')
    return tt_.mjd


def compute_gauss(xarray, avg_y_not_smooht_data, gauss_lines):
    """

    :param xarray: xarray
    :param avg_y_not_smooht_data:
    :param gauss_lines: gauss_lines
    :return: area under signal from spectre
    """
    ampvid = avg_y_not_smooht_data
    velocity = xarray
    indexs = [(np.abs(velocity - float(line))).argmin() for line in gauss_lines]
    mons = [max(ampvid[index - 5:index + 5]) for index in indexs]
    gaussian = [models.Gaussian1D(mons[index], gauss_lines[index],
                                  0.05, bounds={'stddev': (None, 0.15)})
                for index in range(0, len(mons))]

    gg_init = reduce(lambda a, b : a + b, gaussian)
    fit = LevMarLSQFitter()
    gg_fit = fit(gg_init, velocity, ampvid)
    if len(gauss_lines) > 1:
        sts = [models.Gaussian1D(gg_fit[index].amplitude,
                                 gg_fit[index].mean,
                                 gg_fit[index].stddev) for index in range(0, len(gauss_lines))]
        gaussian_areas = []
        gaussiana_amplitudes = [str(gg_fit[index].amplitude).split("=")[-1].replace(")", "")
                                for index in range(0, len(gauss_lines))]
        gaussiana_mean = [str(gg_fit[index].mean).split("=")[-1].replace(")", "")
                          for index in range(0, len(gauss_lines))]
        gaussiana_std = [str(gg_fit[index].stddev).split(",")[1].split("=")[-1]
                         for index in range(0, len(gauss_lines))]

        for st_ in sts:
            gaussian_areas.append(trapz(st_(velocity), velocity))

        return (gaussian_areas, sts, gg_fit, velocity, ampvid, gauss_lines,
                gaussiana_amplitudes, gaussiana_mean, gaussiana_std)

    else:
        sts = models.Gaussian1D(gg_fit.amplitude, gg_fit.mean, gg_fit.stddev)
        gaussian_areas = []
        gaussiana_amplitudes = str(gg_fit.amplitude).split("=")[-1].replace(")", "")
        gaussiana_mean = str(gg_fit.mean).split("=")[-1].replace(")", "")
        gaussiana_std = str(gg_fit.stddev).split(",")[1].split("=")[-1]

        gaussian_areas.append(trapz(sts(velocity), velocity))

        return (gaussian_areas, sts, gg_fit, velocity, ampvid,
                gauss_lines, gaussiana_amplitudes, gaussiana_mean, gaussiana_std)
