#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
create a plot for publications where spectra and monitoring plot are viewed side by side
"""
import datetime
import json

import sys
import os
import argparse

from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import numpy as np
import matplotlib.pyplot as plt
import h5py
from matplotlib.ticker import StrMethodFormatter

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import convert_datetime_object_to_mjd


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''')
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("-b", "--base", help="Base component", type=int, default=-1)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="../config/config.cfg")
    parser.add_argument("-n", "--not_show", type=str, default="")
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

    old_dates = []
    old_log_file_dir = get_configs("paths", "oldprettylogpath")
    for old_log_file in os.listdir(old_log_file_dir):
        with open(old_log_file_dir + old_log_file, "r") as old_log_file2:
            lines = old_log_file2.readlines()
            station = lines[1].split(";")[1]
            if station == "IRBENE":
                line_index = -1
                for line in lines:
                    line_index += 1
                    if line.startswith("Source;"):
                        source = line.split(";")[1].split(",")[0]
                        if source == get_args("source"):
                            data = lines[line_index + 1].split(";")[1]
                            time = lines[line_index + 2].split(";")[1]
                            time_stamp = data + "_" + time
                            old_dates.append(time_stamp)

    old_dates = [convert_datetime_object_to_mjd(datetime.datetime.strptime(d, '%Y-%m-%d_%H:%M:%S')) for d in old_dates]
    component_count = len(get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(","))
    data_file = get_configs("paths", "monitoringFilePath") + get_args("source") + "_" + get_args("line") + ".npy"
    data = np.load(data_file, allow_pickle=True)
    xdata = data[0][0]
    print("Number of observations", len(xdata))

    components = [i for i in range( 1, component_count + 1 )]
    symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#FF33E6']
    velocity = get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(",")

    f, (ax1, ax2) = plt.subplots( 1, 2, gridspec_kw={'width_ratios': [2, 4]}, figsize=(19.8, 9.3) )
    f.tight_layout(pad=2, h_pad=2, w_pad=2, rect=None)
    f.subplots_adjust( wspace=0.1, right=0.876, bottom=0.076, top=0.945)
    f.suptitle(get_configs("Full_source_name", get_args("source")),
               horizontalalignment="center", verticalalignment="center", x=0.4)

    spectr_files = ["cepa_58891.388703703706_IRBENE_1290.h5",
                     "cepa_58836.87358796296_IRBENE_1256.h5",
                     "cepa_58836.889375_IRBENE_1260.h5"]

    symbols_2 = [".", "-.", "-"]
    i = 0

    for spectr_file in spectr_files:
        file = get_configs("paths", "outputFilePath") + "/" + get_args("line") + "/" + spectr_file
        spectr_data = h5py.File(file, 'r')['amplitude_corrected'][()]

        xdata_ = spectr_data[:, 0]
        ydata_ = spectr_data[:, 3]

        ax1.plot(xdata_, ydata_, symbols_2[i])
        i += 1

    ax1.set_xlim(-42.0, -33.0)
    ax1.set_xlabel("Velocity (km sec$^{-1}$)")
    ax1.set_ylabel("Flux density (Jy)")
    ax1.set_ylim(-1.5, 6 )

    ax2.plot([], [], ' ', label="km sec$^{-1}$")

    with open( get_configs("paths", "resultFilePath") + "/" + get_args("source") + "_6668.json", "r") as result_data:
        result_data = json.load(result_data)

    rt32_observation_dates = old_dates

    for observation in result_data:
        if observation.split("_")[-2] == "IRBENE":
            mjd = result_data[observation]["modifiedJulianDays"]
            rt32_observation_dates.append(mjd)

    rt32_x = [rt32 for rt32 in xdata if rt32 in rt32_observation_dates]
    rt32_x_indexies = [xdata.tolist().index(rt32) for rt32 in xdata if rt32 in rt32_observation_dates]

    for component in components:
        index = components.index(component)
        if len(get_args("not_show")) != 0:
            if index + 1 != 3:
                y = data[:, index + 1]
                ax2.plot(xdata, y, symbols[index], color=colors[index],
                          linewidth=0.5,
                          markersize=5,
                          label=str(velocity[index]))
                ax2.errorbar(xdata[0], y[0],
                             yerr=1.5 + 0.05 * y[0],
                             xerr=None, ls='none',
                             ecolor='k')  # 1st point error bar
        else:
            y = data[index + 1]
            rt32_y = np.array(y)[rt32_x_indexies]
            ax2.plot(xdata, y, symbols[index], c=colors[index], linewidth=0.5, markersize=4, label=str( velocity[index]))
            ax2.plot(rt32_x, rt32_y, symbols[index], c=colors[index], linewidth=0.5, markersize=10)
            ax2.errorbar(xdata[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')

    ax1.yaxis.set_ticks_position( 'both' )
    ax1.xaxis.set_ticks_position( 'both' )
    ax1.tick_params(axis="x", direction="in",
                    which="both", length=16,
                    width=2, labelsize=12,
                    rotation=0)
    ax1.tick_params(axis="y", direction="in",
                    which="major", length=16,
                    width=2, labelsize=12)
    ax1.tick_params(axis="y", direction="in",
                    which="minor", length=10,
                    width=1.5)
    ax2.set_ylim(1)
    ax2.set_yscale("log")
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter('{x:.0f}'))
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.tick_params(axis="x", direction="in",
                    which="both", length=16,
                    width=2, labelsize=13,
                    rotation=0)
    ax2.tick_params(axis="y", direction="in",
                    which="major", length=16,
                    width=2, labelsize=12)
    ax2.tick_params(axis="y", direction="in",
                    which="minor", length=10,
                    width=1.5)
    ax2.set_xlabel("MJD")
    ax2.legend( bbox_to_anchor=(1, 0.5, 0.3, 0.3), loc='upper left', borderaxespad=0.5 )
    plt.show()
    f.savefig("monitoring_" + get_args("source") + ".eps", format="eps", dpi=5000, papertype="a4")


if __name__ == "__main__":
    main()
