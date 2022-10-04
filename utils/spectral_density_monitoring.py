#! /usr/bin/python3
# -*- coding: utf-8 -*-
import sys

import os

import argparse
from matplotlib import pyplot as plt

PACKAGE_PARENT = '..'

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from compute_spectral_density import spectral_density


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring areas of a spectrum. ''')
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="line", type=int)
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


def main():
    source = get_args("source")
    line = str(get_args("line"))

    output_path = get_configs("paths", "outputFilePath") + str(line) + "/" + source + "/"
    output_files = os.listdir(output_path)
    areas = []
    mjds = []
    for output_file in output_files:
        areas.append(spectral_density(output_path + output_file))
        output_file_data = output_file.replace(".h5", "").split("_")
        mjds.append(float(output_file_data[1]))

    plt.xlabel('MJD')
    plt.ylabel('Integral flux density (Jy km sec$^{-1}$)')
    plt.scatter(mjds, areas)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
