#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
For definite (low, middle, high) output files for give sources compute spectral density and display histogram
"""
import sys
import os
import argparse
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import h5py

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''For definite (low, middle, high) output files for give sources 
    compute spectral density and display histogram. ''')
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="../config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args


def get_args(key):
    """

    :param key: argument key
    :return: to script passed argument value
    """
    return str(parse_arguments().__dict__[key])


def get_configs(section, key):
    """

    :param section: configuration file section
    :param key: configuration file sections key
    :return: configuration file section key value
    """
    config_file_path = get_args("config")
    config = ConfigParser(config_file_path)
    return config.get_config(section, key)


def get_configs_items():
    """

    :return: None
    """
    config_file_path = "../config/plot.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items("main")


def main():
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=20)
    font_properties = FontProperties()
    font_properties.set_size('small')

    config_items = get_configs_items()
    for key, value in config_items.items():
        rcParams[key] = value

    spectr_files_cepa = ["cepa_58891.388703703706_IRBENE_1290.h5"]
    spectr_files_for_all_sources = [spectr_files_cepa]

    density = []

    for source in spectr_files_for_all_sources:
        for file_name in source:
            file = get_configs("paths", "outputFilePath") + "/6668/" + file_name
            spectr_data = h5py.File(file, 'r')['amplitude_corrected_not_smooht'][()]

            xdata_ = spectr_data[:, 0]
            ydata_ = spectr_data[:, 3]

            area = np.trapz(ydata_, xdata_)
            print(area)

            density.append(area)

    plt.hist(density)
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
