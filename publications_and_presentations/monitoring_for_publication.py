#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Create LATEX tables for publications
"""
import json
import sys
import os
import argparse
from datetime import datetime
from functools import reduce
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.ticker import StrMethodFormatter
import numpy as np


PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import file_len, correct_numpy_read_data, convert_datetime_object_to_mjd, find_nearest_index


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Create LATEX tables for publications. ''')
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("-t0", "--start", help="start date in mjd", type=int, default=-1)
    parser.add_argument("-tn", "--stop", help="stop date in mjd", type=int, default=-1)
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
    both = False
    new = False
    old = False
    configuration_items = get_configs_items()
    for key, value in configuration_items.items():
        rcParams[key] = value

    result_file_name = get_args("source") + "_" + get_args("line") + ".json"
    result_file_path = get_configs("paths", "resultFilePath")

    correction_file = get_configs("parameters", "amplitude_correction_file") + ".npy"
    correction_data = np.load(correction_file)
    correction_mjd = correction_data[:, 0]
    correction_factor = correction_data[:, 1]

    rms = []
    with open(result_file_path + result_file_name) as result_data:
        result = json.load(result_data)
        for experiment in result:
            rms.append(result[experiment]['rms_avg'])

    old_monitoring_file = get_configs("paths", "oldMonitoringFilePath") + get_args("source") + ".dat"
    new_monitoring_file = get_configs("paths", "monitoringFilePath") + \
                          get_args("source") + "_" + \
                          get_args("line") + ".npy"
    component_count = len(get_configs("velocities",
                                      get_args("source") + "_" +
                                      get_args("line")).replace(" ", "").split(","))
    velocity = get_configs("velocities",
                           get_args("source") + "_" +
                           get_args("line")).replace(" ", "").split(",")
    components = [i for i in range(1, component_count + 1)]

    old_x = []
    if os.path.isfile(old_monitoring_file) and os.path.isfile(new_monitoring_file):
        new_data = np.load(new_monitoring_file, allow_pickle=True)
        new_x = new_data[0][0]
        old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
            (file_len(old_monitoring_file), component_count + 1))
        old_x = correct_numpy_read_data(old_data[:, [0]])
        old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
        old_data[:, 0] = old_x
        x = old_x + list(new_x)
        old_data_tmp = []
        for tmp in range(0, len(new_data)):
            old_data_tmp.append([])
        tmp2 = 0
        for j in range(0, old_data.shape[1]):
            for i in range(0, old_data.shape[0]):
                old_data_tmp[tmp2].append(old_data[i][j])
            old_data_tmp[tmp2] = np.array(old_data_tmp[tmp2]).reshape(old_data.shape[0], )
            tmp2 += 1
        old_data = np.array(old_data_tmp)
        new_data[0] = new_data[0][0]
        data = []

        for tmp3 in range(0, old_data.shape[0]):
            data_tmp = np.concatenate((old_data[tmp3], new_data[tmp3]), axis=0)
            data.append(data_tmp)
        data = np.array(data)
        both = True

    elif os.path.isfile(old_monitoring_file):
        old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape(
            (file_len(old_monitoring_file), component_count + 1))
        old_x = correct_numpy_read_data(old_data[:, [0]])
        old_x = [convert_datetime_object_to_mjd(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
        old_data[:, 0] = old_x
        data = old_data
        x = list(old_x)
        old = True

    else:
        new_data = np.load(new_monitoring_file, allow_pickle=True)
        new_x = new_data[0][0]
        data = new_data
        x = list(new_x)
        new = True

    if int(get_args("start")) != -1 and int(get_args("stop")) != -1:
        a = (np.abs(np.array(x) - int(get_args("start")))).argmin()
        b = (np.abs(np.array(x) - int(get_args("stop")))).argmin()
        x = x[a:b]

    print("total time in years", (np.max(x) - np.min(x)) / 365)
    print("Nmbers of observations", len(x))
    print("Observations per month", (len(x) / ((np.max(x) - np.min(x)) / 365)) / 12)

    fig = plt.figure("Monitoring", figsize=(4, 3), dpi=75)
    ax1 = fig.add_subplot(111)

    symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    y_ticks = []

    variances = dict()
    variability_index = dict()
    variances_normal = dict()
    fluctuation_index = dict()
    xhi_sqer_red = dict()
    xhi_sqer = dict()
    result_org = [x]

    print("\n")
    print("\hline")
    print("\multicolumn{4}{l}{" + get_configs("Full_source_name", get_args("source")) +
          " (MJD\\textsubscript{{s}}={:d} $T\\textsubscript{{s}}$= "
          "{:.3f},}}\\\\".format(int(np.min(x)), (np.max(x) - np.min(x)) / 365))

    print("\multicolumn{{4}}{{l}}{{$N$={:d}, $C(month^{{-1}})$={:.3f})}}"
          " \\\\".format(len(x), (len(x) / ((np.max(x) - np.min(x)) / 365)) / 12))
    print("\hline")

    ax1.plot([], [], ' ', label="km sec$^{-1}$")
    if len(old_x) > 0:
        rms_old = [np.mean(rms)] * len(old_x)
        rms_old.extend(rms)
        rms = rms_old

    factor = []
    for m in range(0, len(x)):
        correction_index = find_nearest_index(correction_mjd, x[m])
        factor.append(correction_factor[correction_index])

    for component in components:

        index = components.index(component)
        if old:
            y = data[:, index + 1]
        elif both:
            y = data[index + 1, :]
        else:
            y = data[index + 1]

        y = np.array([np.float128(yi) for yi in y]).clip(min=0)
        if int(get_args("start")) != -1 and int(get_args("stop")) != -1:
            y = y[a:b]
        N = len(y)

        # ax1.scatter(x, y/y[1], color=colors[index], marker=symbols[index])
        # ax1.scatter(x, y/y[1], color=colors[index], marker=symbols[index], label=str(velocity[index]))

        ax1.scatter(x, y, color=colors[index], marker=symbols[index])
        ax1.scatter(x, y, color=colors[index], marker=symbols[index], label=str(velocity[index]))

        np.savetxt(get_args("source") + "_" + str(index) + ".txt", np.vstack((x, y)).T, fmt='%8.1f  %8.1f')
        ax1.errorbar(x[0], y[0], yerr=1.5 + 0.1 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar
        result_org.append(y)
        variances[component] = reduce(lambda x_, y_: x_ + y_, [((i - np.mean(y)) / np.std(y)) ** 2 for i in y])

        error = []
        for i in range(0, len(y)):
            if 1 - factor[i] < 0.05:
                error.append(2 * rms[i] + y[i] * 0.05)
            else:
                error.append(2 * rms[i] + y[i] * (1 - factor[i]))

        variability_index[component] = ((np.max(y) - error[list(y).index(np.max(y))]) - (
                    np.min(y) + error[list(y).index(np.min(y))])) \
                                       / ((np.max(y) - error[list(y).index(np.max(y))]) + (
                    np.min(y) + error[list(y).index(np.min(y))]))

        variances_normal[component] = variances[component] * (1 / N - 1)

        fluctuation_index[component] = np.sqrt(
            np.abs((N / reduce(lambda x_, y_: x_ + y_, [error[list(y).index(i)] ** 2 for i in y])) *
                   ((reduce(lambda x_, y_: x_ + y_,
                            [i ** 2 * (error[list(y).index(i)]) ** 2 for i in y]) -
                     np.mean(y) * reduce(lambda x_, y_: x_ + y_,
                                         [i * error[list(y).index(i)] ** 2 for i in y]))
                    / (N - 1)) - 1)) / np.mean(y)

        xhi_sqer_red[component] = reduce(lambda x_, y_: x_ + y_,
                                         [((i - np.mean(y)) / error[list(y).index(i)]) ** 2 for i in y]) / (N - 1)

        xhi_sqer[component] = reduce(lambda x_, y_: x_ + y_, [((i - np.mean(y)) / (1.9 + 0.2 * i)) ** 2 for i in y])

        xhi_sqer[component] = sum(((i - np.mean(y)) / (1.9 + 0.2 * i)) ** 2 for i in y)

        v = velocity[index]
        # print(get_configs("Full_source_name", get_args("source")) + " & " + "{} &  {} & {} & {:.1f} & {:3} &  {:.3f} & {:.3f} & {:.3f}\\\\".
        # format(int(x[0]), int(x[-1]), N,  len(x) / ((np.max(x) - np.min(x)) / 365) / 12,  v, np.mean(y), variability_index[component], fluctuation_index[component]))

        print(get_configs("Full_source_name",
                          get_args("source")) + "  {}  {}  {}  {:.1f}  {:3}  {:.3f}  {:.3f}  {:.3f}  {:.3f} {:.3f}".
              format(int(x[0]), int(x[-1]), N, len(x) / ((np.max(x) - np.min(x)) / 365) / 12, v, np.mean(y),
                     variability_index[component], fluctuation_index[component], xhi_sqer_red[component],
                     xhi_sqer[component]))

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
    print("variability indexies", variability_index)
    print("Start time", np.min(x))

    plt.title(get_configs("Full_source_name", get_args("source")))
    ax1.set_xlabel("MJD")
    ax1.set_ylabel("Flux density (Jy)")

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
    plt.legend(bbox_to_anchor=(1, 0.5, 0.3, 0.3), loc='upper left', borderaxespad=0.5)
    plt.show()

    result_calib = [x]
    if int(get_args("base")) >= 0:

        for component in components:
            index = components.index(component)
            y = np.float128(data[index + 1])[0]
            base = np.float128(data[[int(get_args("base"))]])[0]
            plt.plot(x, y / base, symbols[index] + colors[index], linewidth=0.5, markersize=5)
            result_calib.append(y / base)

        plt.grid(False)
        plt.show()

    result_norm_file_name = get_configs("paths", "monitoringFilePath") + \
                            get_args("source") + "_" + get_args("line") + "_norm"
    np.save(result_norm_file_name, np.transpose(result_calib))

    sys.exit(0)


if __name__ == "__main__":
    main()