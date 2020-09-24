#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Create LATEX tables for publications
"""
import sys
import os
import argparse
from datetime import datetime
from functools import reduce
import matplotlib.pyplot as plt
from astropy.time import Time
from matplotlib import rcParams
from matplotlib.ticker import StrMethodFormatter
import numpy as np

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import file_len, correct_numpy_read_data, convert_datetime_object_to_mjd


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Create LATEX tables for publications. ''')
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("-b", "--base", help="Base component", type=int, default=-1)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="../config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args


def get_configs_items():
    """

    :return: None
    """
    config_file_path = "../config/plot.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items("main")


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
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    old_monitoring_file = get_configs("paths", "oldMonitoringFilePath") + get_args("source") + ".dat"
    new_monitoring_file = get_configs("paths", "monitoringFilePath") + \
                          get_args("source") + "_" + \
                          get_args("line") + ".npy"

    print(new_monitoring_file)
    component_count = len(get_configs("velocities",
                                      get_args("source") + "_" +
                                      get_args("line")).replace(" ", "").split(","))
    velocity = get_configs("velocities",
                           get_args("source") + "_" +
                           get_args("line")).replace(" ", "").split(",")
    components = [i for i in range(1, component_count + 1)]

    if os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
        new_data = np.load(new_monitoring_file, allow_pickle=True)
        new_x = new_data[0][0]
        old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
            (file_len(old_monitoring_file), component_count + 1))
        old_x = correct_numpy_read_data(old_data[:, [0]])
        old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
        old_data[:, [0]] = old_x[0]
        x = old_x + list(new_x)
        old_data_tmp = []
        for tmp in range(0, len(new_data)):
            old_data_tmp.append([])
        tmp2 = 0
        for j in range(0, old_data.shape[1]):
            for i in range(0, old_data.shape[0]):
                old_data_tmp[tmp2].append(old_data[i][j])
            old_data_tmp[tmp2] = np.array(old_data_tmp[tmp2]).reshape(old_data.shape[0],)
            tmp2 += 1
        old_data = np.array(old_data_tmp)
        new_data[0] = new_data[0][0]
        data = []

        for tmp3 in range(0, old_data.shape[0]):
            data_tmp = np.concatenate((old_data[tmp3], new_data[tmp3]), axis=0)
            data.append(data_tmp)

        data = np.array(data)
        old = True

    elif os.path.isfile(old_monitoring_file):
        old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
            (file_len(old_monitoring_file), component_count + 1))
        old_x = correct_numpy_read_data(old_data[:, [0]])
        old_x = [convert_datetime_object_to_mjd( datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
        old_data[:, [0]] = old_x[0]
        data = old_data
        x = list(old_x)
        old = True

    else:
        new_data = np.load(new_monitoring_file, allow_pickle=True)
        new_x = new_data[0][0]
        data = new_data
        x = list(new_x)
        old = False

    print("total time in years", (np.max(x) - np.min(x)) / 365)
    print("Nmbers of observations", len(x))
    print("Observations per month", (len(x) / ((np.max(x) - np.min(x)) / 365)) / 12)

    fig = plt.figure("Monitoring", figsize=(4, 3), dpi=75)
    ax1 = fig.add_subplot(111)

    symbols = ["*-", "o-", "v-", "^-", "<-", ">-", "1-", "2-", "3-", "4-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    y_ticks = []

    variances = dict()
    variability_indexies = dict()
    variances_normal = dict()
    fuction_index = dict()
    result_org = [x]

    print("\n")
    print("\hline")
    print("\multicolumn{4}{l}{" + get_configs("Full_source_name", get_args("source")) +
          "MJD\\textsubscript{{s}}= {:.3f} $T\\textsubscript{{s}}$= "
          "{:.3f},}}\\\\".format(np.min(x), (np.max(x) - np.min(x)) / 365))
    print("\multicolumn{{4}}{{l}}{{$N$={:d}, $C(month^{{-1}})$={:.3f})}}"
          " \\\\".format(len(x), (len(x) / ((np.max(x) - np.min(x)) / 365)) / 12))
    print("\hline")

    for component in components:
        index = components.index(component)
        y = data[index + 1]
        y = [np.float128(yi) for yi in y]
        N = len(y)
        ax1.plot(x, y, symbols[index] + colors[index], linewidth=0.5, markersize=5)
        ax1.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar
        result_org.append(y)
        variances[component] = reduce(lambda x, y: x + y, [((i - np.mean(y)) / np.std(y)) ** 2 for i in y] )
        variability_indexies[component] = ((np.max(y) - np.std(y)) - (np.min(y) + np.std(y)))\
                                          / ((np.max(y) - np.std(y)) + (np.min(y) + np.std(y)))
        variances_normal[component] = variances[component] * (1 / N - 1)
        fuction_index[component] = np.sqrt((N / reduce(lambda x, y: x + y, [(1.5 + 0.05 * i) ** 2 for i in y] )) *
                                           ((reduce(lambda x, y: x + y, [i ** 2 * (1.5 + 0.05 * i) ** 2 for i in y]) - np.mean(y) *
                                            reduce(lambda x, y: x + y, [i * (1.5 + 0.05 * i) ** 2 for i in y]))
                                            / (N - 1)) - 1) / np.mean(y)
        v = velocity[index]
        print("{:3} &  {:.3f} & {:.3f} & {:.3f}\\\\".
              format(v, np.mean(y), variability_indexies[component], fuction_index[component]))

        y_min = np.min(y)
        if y_min < 0:
            ymin = 10

        ytick = 1
        n = 1
        while ytick < np.max(y) + 100:
            if ytick > y_min:
                y_ticks.append(ytick)
            n += 1
            ytick = 10 ** n

        y_ticks.append(np.max(y))

    print("\hline")
    print("\n")
    print("variances", variances)
    print("variances_normal", variances_normal)
    print("variability indexies", variability_indexies)
    print("Start time", np.min(x))

    plt.title(get_configs("Full_source_name", get_args("source")))
    ax1.set_xlabel("MJD" )
    ax1.set_ylabel("Flux density (Jy)")

    t = Time(x, format='mjd', scale='utc', out_subfmt='date')
    t.format = 'isot'

    fig.autofmt_xdate()
    plt.yscale("log")
    ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(axis="x", direction="in", which="both",
                    length=16, width=2, labelsize=12, rotation=0)
    ax1.tick_params(axis="y", direction="in", which="major",
                    length=16, width=2, labelsize=12)
    ax1.tick_params(axis="y", direction="in",
                    which="minor", length=10, width=1.5)

    ax1.grid(False)
    plt.xticks(rotation=0, ha='right')
    plt.show()

    result_calib = [x]
    if int(get_args("base")) >= 0:

        for component in components:
            index = components.index(component)
            y = np.float128(data[index + 1])[0]
            base = np.float128(data[[int(get_args("base"))]])[0]
            plt.plot(x, y/base, symbols[index] + colors[index], linewidth=0.5, markersize=5)
            result_calib.append(y / base)

        plt.grid(False)
        plt.show()

    result_norm_file_name = get_configs("paths", "monitoringFilePath") + \
                            get_args("source") + "_" + get_args("line") + "_norm"
    np.save(result_norm_file_name, np.transpose(result_calib))

    sys.exit(0)


if __name__ == "__main__":
    main()
