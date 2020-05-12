"""
common used fuctions
"""

import numpy as np
from astropy.time import Time


def file_len(file_name):
    """

    :param fname: file name
    :return: rows of file
    """
    i = 0
    with open(file_name) as file:
        for i, _ in enumerate(file):
            pass
    return i + 1


def indexies(array, value):
    """

    :param array: nupmy array
    :param value: value to search
    :return: indexes in array for value
    """
    indexes = list()
    for i in range(0, len(array) -1):
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
