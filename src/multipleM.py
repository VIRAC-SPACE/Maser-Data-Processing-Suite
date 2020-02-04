#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

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
        lines[source] = {"y_data":[]}

    for file in monitoring_files:
        date = np.loadtxt(file, usecols=(0,), unpack=True)
        source = file.split("/")[-1].split(".")[0]
        lines[source]["date"] = date
        colomn_nr = len(get_configs("velocities", source).split(","))

        for c in range(1, colomn_nr+1):
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

        for y_data in lines[source]["y_data"]:
            plt.plot(date, y_data)

    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()