#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
For given output files compute compute spectral density
"""
import h5py
import numpy as np
import argparse


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''For given output files compute compute spectral density.''')
    parser.add_argument("file_names", help="full path to files", type=str, nargs='+')
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args


def get_args(key):
    """

    :param key: argument key
    :return: to script passed argument value
    """
    return str(parse_arguments().__dict__[key])


def spectral_density(file):
    file = h5py.File(file, 'r')
    if "amplitude_corrected_not_smooht" in file:
        spectr_data = file['amplitude_corrected_not_smooht'][()]
        xdata = spectr_data[:, 0]
        ydata = spectr_data[:, 3]
        area = np.trapz(ydata, xdata)
        return area


def main():
    files = get_args("file_names").replace("[", "").replace("]", "").replace("'", "").split(",")
    files = [file.strip() for file in files]
    for file in files:
        print(file, spectral_density(file))


if __name__ == "__main__":
    main()
