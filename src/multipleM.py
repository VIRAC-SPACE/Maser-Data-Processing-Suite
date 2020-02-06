#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

from parsers._configparser import ConfigParser


def parse_arguments():
    parser = argparse.ArgumentParser(description='''Monitoring multiple sources. ''', epilog="""MMonitor.""")
    parser.add_argument("sources", help="Sources Names", type=str, nargs='+')
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


def read_monitoring_files(monitoring_files, sources):
    lines = {}

    for source in sources:
        lines[source] = {"y_data": []}

    for file in monitoring_files:
        date = np.loadtxt(file, usecols=(0,), unpack=True)
        source = file.split("/")[-1].split(".")[0]
        lines[source]["date"] = date
        column_nr = len(get_configs("velocities", source).split(","))

        for c in range(1, column_nr+1):
            lines[source]["y_data"].append(np.loadtxt(file, usecols=(c,), unpack=True))

    return lines


def main():
    sources = get_args("sources").replace("[", "").replace("]", "").replace("'", "").split(",")
    sources = [s.strip() for s in sources]
    monitoring_files = []
    monitoring_dir = get_configs("paths", "monitoringFilePath")

    for source in sources:
        monitoring_files.append(monitoring_dir + source + ".out")

    lines = read_monitoring_files(monitoring_files, sources)

    for source in sources:
        date = lines[source]["date"]
        velocities = get_configs("velocities", source).split(",")
        velocities = [v.strip() for v in velocities]

        i = 0
        for y_data in lines[source]["y_data"]:
            z = np.polyfit(date, y_data, 1 )
            p = np.poly1d(z)
            plt.plot(date, y_data, label=source+" velocity " + velocities[i] + " " + "y=%.6fx+(%.6f)"%(z[0],z[1]) + " " + "R^2 = " + str(r2_score(y_data, p(date))))
            plt.plot(date, p(date), "r--" )
            i += 1

    plt.legend()
    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()
