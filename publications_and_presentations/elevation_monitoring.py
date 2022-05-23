#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
create a plot for publications where elevation and monitoring plot are plotted in one graph
"""

import sys
import os
import argparse

import datetime
from datetime import datetime

from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time

from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes


PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import find_nearest_index


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''')
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


def get_configs_items():
    """

    :return: None
    """
    config_file_path = "../config/plot.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items("main")


def main():
    """

    :return: None
    """
    font_properties = FontProperties()
    font_properties.set_size('small')

    config_items = get_configs_items()
    for key, value in config_items.items():
        rcParams[key] = value

    source = get_args("source")
    line = str(get_args("line"))
    component_count = len(get_configs("velocities", get_args("source") + "_" +
                          get_args("line")).replace(" ", "").split(","))

    components = [i for i in range(1, component_count + 1)]

    log_path = get_configs("paths", "logPath") + "SDR/"
    log_files = [file for file in os.listdir(log_path) if file.startswith(source + "_") and "f" + line in file]
    #log_files = [log_files[0]]

    log_file_data = []
    for log_file in log_files:
        iteration_number = log_file.split(".")[0].split("_")[-1]
        try:
            logs = LogReaderFactory.getLogReader(LogTypes.SDR, get_configs("paths", "logPath") + "SDR/" + log_file,
                                                 get_configs("paths", "prettyLogsPath") + source + "_" +
                                                 str(iteration_number)).getLogs()
            log_file_data.append(logs)
        except IndexError:
            print("Wrong log file")

    elevations = []
    mjd = []

    for log in log_file_data:
        try:
            mjd.append(Time(datetime.strptime(log["1s0"]["date"], '%Y-%m-%dT%H:%M:%S').isoformat(),
                            format='isot', scale='utc').mjd)
            elevations.append(np.float(log["1s0"]["AzEl"][1]))
        except KeyError:
            print("Wrong log file")

    new_monitoring_file = get_configs("paths", "monitoringFilePath") + \
                          get_args("source") + "_" + get_args("line") + ".npy"
    new_data = np.load(new_monitoring_file, allow_pickle=True)
    new_x = new_data[0][0]

    both_index = []
    for m in mjd:
        both_index.append(find_nearest_index(new_x, m))

    mjd_both = []
    for index in both_index:
        mjd_both.append(mjd[index])

    output = np.zeros((len(mjd_both), len(components) + 2))
    output[:,0] = mjd_both
    output[:,len(components) + 1] = elevations

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)

    symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#FF33E6']
    velocity = get_configs("velocities", get_args("source") + "_" + get_args("line")).replace(" ", "").split(",")
    for component in components:
        index = components.index(component)
        y = new_data[index + 1].astype('float64')

        for index_ in range(0, len(both_index)):
            output[index_, index + 1] = y[both_index[index_]]

        ax1.plot(new_x, y, symbols[index], color=colors[index], linewidth=0.5,
                 markersize=5, label=str(velocity[index]))

    ax2 = ax1.twinx()
    ax2.scatter(mjd, elevations, label="elevation")

    ax1.set_xlabel("MJD")
    ax1.set_ylabel("Velocity (km sec$^{-1}$)")
    ax1.legend()

    ax2.set_ylabel("Elevation")
    ax2.legend()

    plt.show()

    np.savetxt("output.dat", output)


if __name__ == "__main__":
    main()
    sys.exit(0)
