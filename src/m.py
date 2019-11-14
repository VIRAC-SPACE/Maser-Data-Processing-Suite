#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os

import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.ticker import StrMethodFormatter, NullFormatter
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import rcParams

import numpy as np
import argparse
from datetime import datetime
from astropy.time import Time
from functools import reduce

from help import *
from parsers._configparser import ConfigParser


def get_configs_items():
    config_file_path = "config/plot.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(config_file_path)
    return config.getItems("main")

def parse_arguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("-b", "--base", help="Base component", type=int, default=-1)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args


def get_args(key):
    args = parse_arguments()
    return str(args.__dict__[key])


def get_configs(key, value):
    config_file_path = get_args("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(config_file_path)
    return config.getConfig(key, value)


def main():
    config_items = get_configs_items()
    for key, value in config_items.items():
        rcParams[key] = value

    component_count = len(get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(","))
    data_file = get_configs("paths", "monitoringFilePath") + get_args("source") + ".out"
    data = np.fromfile(data_file, dtype="float64", count=-1, sep=" ").reshape((file_len(data_file), component_count + 1))
    x = correctNumpyReadData(data[:, [0]])

    components = [i for i in range(1, component_count + 1)]
    symbols = ["*-", "o-", "v-", "^-", "<-", ">-", "1-", "2-", "3-", "4-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']


    f, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 4]}, dpi=75, figsize=(20, 12))
    f.tight_layout(pad=2, h_pad=2, w_pad=2, rect=None)
    f.suptitle(get_configs("Full_source_name", get_args("source")), horizontalalignment="center", verticalalignment="center", x=0.4)

    #ax1 = plt.subplot(121)
    specter_files = ["cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat", "cepa_22_24_38_23_Jul_2018_IRBENE16_46.dat", "cepa_22_08_44_29_Aug_2018_IRBENE16_121.dat"]

    symbols_2 = [".", "-.", "-"]
    i = 0

    for specter in specter_files:
        file = get_configs("paths", "outputFilePath") + "/" + get_args("source") + "/6668/" + specter
        spectre_data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file), 4))

        x_ = correctNumpyReadData(spectre_data[:, [0]])
        y_ = correctNumpyReadData(spectre_data[:, [3]])

        ax1.plot(x_, y_, symbols_2[i])
        i += 1

    ax1.set_xlabel("Velocity (km sec$^{-1}$)")
    ax1.set_ylabel("Flux density (Jy)")


    #ax2 = plt.subplot(122)
    for component in components:
        index = components.index(component)
        y = correctNumpyReadData(data[:, [index + 1]])
        ax2.plot(x, y, symbols[index] + colors[index], linewidth=0.5, markersize=5)
        ax2.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar

    ax1.yaxis.set_minor_formatter(NullFormatter())
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=12,rotation=0)  # MJD atzimju formâts
    ax1.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax1.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats

    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))  # clasic log scale!
    ax2.yaxis.set_minor_formatter(NullFormatter())
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=12,rotation=0)  # MJD atzimju formâts
    ax2.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax2.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats

    ax2.set_yscale("log")
    ax2.set_xlabel("MJD")

    plt.show()
    f.savefig("/home/janis/Desktop/monitoring.eps", format="eps", dpi=5000)

if __name__ == "__main__":
    main()