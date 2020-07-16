#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Display monitoring for multiple sources
"""
import argparse
import re
from tabulate import tabulate
import numpy as np

from parsers.configparser_ import ConfigParser
from utils.help import divide_list_with_num


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring multiple sources. ''')
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("--sources", help="Sources Names", type=str, nargs='+',
                        default="[g32p745, w51, g59p783, on1, s252, ngc7538, w3oh]")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str,
                        default="config/config.cfg")
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


def print_stats(stats, labels):
    headers = ["source", "nobs", "min", "max", "mean", "variance", "skewness", "kurtosis"]
    data = [headers]
    i = 0
    for s in stats:
        stats_results = [labels[stats.index(s)]]
        stats_tmp = str(s).replace("DescribeResult", "").replace("(", "").replace(")", "").replace(",", "").split(" ")
        stats_tmp = [re.sub("[^0-9.]", "", s.split("=")[-1])[0:5] for s in stats_tmp]
        stats_results.extend(stats_tmp)
        data.append(stats_results)
        i += 1

    print(tabulate(data))


def read_monitoring_files(monitoring_files, sources):
    lines = {}

    velocities_to_plot_for_source = {}
    for source in sources:
        lines[source] = {"y_data": []}

        if source == "g32p745":
            velocities_tmp = ["30.49",  "39.18"]
        elif source == "w51":
            velocities_tmp = ["59.29"]
        elif source == "g59p783":
            velocities_tmp = ["19.2"]
        elif source == "on1":
            velocities_tmp = ["14.64"]
        elif source == "s252":
            velocities_tmp = ["10.84"]
        elif source == "ngc7538":
            velocities_tmp = ["-58.04"]
        else:
            velocities_tmp = ["-44.6"]

        velocities_to_plot_for_source[source] = velocities_tmp

    for file in monitoring_files:
        date = np.load(file, allow_pickle=True)
        source = file.split("/")[-1].split(".")[0].split("_")[0]
        lines_to_plot = velocities_to_plot_for_source[source]
        tmp = get_configs("velocities", source + "_6668").split(",")
        tmp = [v.strip().replace(" ", "") for v in tmp]
        idexie_for_lines_to_plot = []
        for tmpl in lines_to_plot:
            idexie_for_lines_to_plot.append(tmp.index(tmpl))

        lines[source]["date"] = date
        column_nr = len(get_configs("velocities", source + "_6668").split(","))

        for c in range(1, column_nr+1):
            if c - 1 in idexie_for_lines_to_plot:
                tmp = divide_list_with_num(np.load(file, allow_pickle=True), np.mean(np.load(file, allow_pickle=True)))
                lines[source]["y_data"].append(tmp)

    return lines


def main():
    sources = get_args("sources").replace("[", "").replace("]", "").replace("'", "").split(",")
    sources = [s.strip() for s in sources]
    monitoring_files = []
    monitoring_dir = get_configs("paths", "monitoringFilePath")

    for source in sources:
        monitoring_files.append(monitoring_dir + source + "_" + get_args("line") + ".npy")

    lines = read_monitoring_files(monitoring_files, sources)


if __name__ == "__main__":
    main()
