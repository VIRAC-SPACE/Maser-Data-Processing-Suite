#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Qt5Agg")
import numpy as np
import argparse
from datetime import datetime

from help import *
from parsers._configparser import ConfigParser


def parse_arguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("source", help="Source Name", type=str)
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
    old_monitoring_file = "old_monitoring/" + get_args("source") + ".dat"
    new_monitoring_file = "monitoring/" + get_args("source") + "_6668" + ".txt"
    compunet_count = len(get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(","))
    velocity = get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(",")
    old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape((file_len(old_monitoring_file), compunet_count + 1))
    new_data = np.fromfile(new_monitoring_file, dtype="float64", count=-1, sep=" ").reshape((file_len(new_monitoring_file), compunet_count + 1))
    new_x = correctNumpyReadData(new_data[:, [0]])
    old_x = correctNumpyReadData(old_data[:, [0]])
    old_x = [convertDatetimeObjectToMJD(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
    x = old_x + list(new_x)
    components = [i for i in range(1, compunet_count + 1)]
    old_data[:, [0]] = old_x[0]

    data = np.array(np.concatenate((old_data, new_data), axis=0), dtype="float64")
    print("total time in years", (np.max(x) - np.min(x)) / 365)

    fig = plt.figure("Monitoring")
    ax1 = fig.add_subplot(111)
    for component in components:
        index = components.index(component)
        ax1.plot(x, correctNumpyReadData(data[:, [index + 1]]), label=velocity[index])
        ax1.legend()
        ax1.grid(True)

    plt.title(get_configs("Full_source_name", get_args("source")))
    ax1.set_xlabel("MJD")
    ax1.set_ylabel("Flux density (Jy)")
    ax2 = ax1.twiny()
    #ax2.set_xticks(new_tick_locations)
    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()
