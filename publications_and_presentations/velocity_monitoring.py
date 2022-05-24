#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
create a plot for publications where velocity and monitoring plot are plotted in one graph
"""

import sys
import os
import argparse

from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import h5py


PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


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
    velocities = get_configs("velocities", get_args("source") + "_" + get_args("line")).replace(" ", "").split(",")
    output_path = get_configs("paths", "outputFilePath") + str(line) + "/" + source + "/"
    output_files = os.listdir(output_path)
    index_range_for_local_maxima = int(get_configs('parameters', "index_range_for_local_maxima"))

    vel_monitoring = {component: [] for component in velocities}
    mjd = []
    for output_file in output_files:
        output_data = h5py.File(output_path + output_file, "r")
        if "amplitude_corrected" in output_data:
            mjd.append(float(output_file.split("_")[1]))
            total_results = output_data["amplitude_corrected"]
            velocity = total_results[: , 0]
            amplitude = total_results[: , 1]
            for component in velocities:
                index_for_local_maxima = (np.abs(velocity - float(component))).argmin()
                max_amplitudes = []
                for i in range(index_for_local_maxima - index_range_for_local_maxima,
                               index_for_local_maxima + index_range_for_local_maxima):
                    max_amplitudes.append(amplitude[i])
                max_amplitude = np.max(max_amplitudes)
                vel_monitoring[component].append(velocity[np.where(amplitude == max_amplitude)[0][0]])

        else:
            print(output_file + " wrong output file")

    symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#FF33E6']
    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(16, 16), dpi=150)
    for vel in velocities:
        index = velocities.index(vel)
        ax1.plot(mjd, vel_monitoring[vel], symbols[index], color=colors[index], linewidth=0.5, markersize=5,
                 label=str(vel))

    ax1.set_xlabel("MJD")
    ax1.set_ylabel("Velocity (km sec$^{-1}$)")
    ax1.legend()
    plt.show()


if __name__ == "__main__":
    main()
    sys.exit(0)
