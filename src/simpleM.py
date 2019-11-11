#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os

import matplotlib
matplotlib.use("Qt5Agg")

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


def getConfigsItems():
    configFilePath = "config/plot.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getItems("main")

configItems = getConfigsItems()
for key, value in configItems.items():
    rcParams[key] = value

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
    years = mdates.YearLocator()
    months = mdates.MonthLocator()
    weeks = mdates.WeekdayLocator()

    years_fmt = mdates.DateFormatter('%Y')

    old_monitoring_file = get_configs("paths", "oldMonitoringFilePath") + get_args("source") + ".dat"
    new_monitoring_file = get_configs("paths", "monitoringFilePath") + get_args("source") + "_6668" + ".txt"
    compunet_count = len(get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(","))
    velocity = get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(",")

    new_data = np.fromfile(new_monitoring_file, dtype="float64", count=-1, sep=" ").reshape((file_len(new_monitoring_file), compunet_count + 1))
    new_x = correctNumpyReadData(new_data[:, [0]])
    components = [i for i in range(1, compunet_count + 1)]

    if os.path.isfile(old_monitoring_file):
        old_data = np.loadtxt(old_monitoring_file, dtype=str).reshape((file_len(old_monitoring_file), compunet_count + 1))
        old_x = correctNumpyReadData(old_data[:, [0]])
        old_x = [convertDatetimeObjectToMJD(datetime.strptime(x, "%Y-%m-%d%H:%M:%S")) for x in old_x]
        old_data[:, [0]] = old_x[0]
        x = old_x + list(new_x)

        data = np.array(np.concatenate((old_data, new_data), axis=0), dtype="float64")
        old = True

    else:
        data = new_data
        x = list(new_x)
        old = False

    #x = [int(str(int(round(qwerty)))[0:3].ljust(5, '0')) for qwerty in x]
    #x = list(set(x))
    #print(x)
    print("total time in years", (np.max(x) - np.min(x)) / 365)

    fig = plt.figure("Monitoring", figsize=(16, 9))
    ax1 = fig.add_subplot(111)

    symbols = ["*-", "o-", "v-", "^-", "<-", ">-", "1-", "2-", "3-", "4-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    yTicks = []

    variances = dict()
    variability_indexies = dict()

    result_org = [x]
    ax1.plot([], [], ' ', label="km sec$^{-1}$")  # viltîba lai dabûtu km/s pirms series labels
    for component in components:
        index = components.index(component)
        y = correctNumpyReadData(data[:, [index + 1]])
        ax1.plot(x, y, symbols[index] + colors[index], label=velocity[index], linewidth=0.5, markersize=5)
        ax1.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar
        result_org.append(y)
        sum = lambda x, y: x + y
        variances[component] = reduce(sum, [((i - np.mean(y)) / np.std(y)) ** 2 for i in y])
        variability_indexies[component] = ((np.max(y) - np.std(y)) - (np.min(y) + np.std(y))) / ((np.max(y) - np.std(y)) + (np.min(y) + np.std(y)))

        y_min = np.min(y)
        if y_min < 0:
            ymin = 10

        ytick = 1
        n = 1
        while ytick < np.max(y) + 100:
            if ytick > y_min:
                yTicks.append(ytick)
            n += 1
            ytick = 10 ** n

        yTicks.append(np.max(y))

    print("variances", variances)
    print("variability indexies", variability_indexies)
    plt.title(get_configs("Full_source_name", get_args("source")), fontsize=18)
    ax1.set_xlabel("MJD", fontsize=15)
    ax1.set_ylabel("Flux density (Jy)", fontsize=15)

    ax2 = ax1.twiny()
    #ax2.set_xlabel("Date", fontsize=15)
    t = Time(x,  format='mjd', scale='utc', out_subfmt='date')
    t.format = 'isot'
    newvalues = [datetime.strptime(i.value, "%Y-%m-%d") for i in t]

    newpos = [p for p in range(0, len(t))]
    ax2.set_xticks(newpos)
    ax2.set_xticklabels(newvalues)
    ax2.xaxis.set_major_locator(years)
    ax2.xaxis.set_major_formatter(years_fmt)
    # ax2.xaxis.set_minor_locator(weeks)
    ax2.autoscale()

    # round to nearest years.
    datemin = np.datetime64(newvalues[0], 'D')
    datemax = np.datetime64(newvalues[-1], 'D')
    ax2.set_xlim(datemin, datemax)

    # format the coords message box
    ax2.format_xdata = mdates.DateFormatter('%Y-%m')
    ax2.grid(False)
    ax1.autoscale()
    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    fig.autofmt_xdate()

    plt.yscale("log")
    ax1.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))  # clasic log scale!
    #ax1.yaxis.set_minor_formatter(NullFormatter())
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=12, rotation=0)  # MJD atzimju formâts
    ax1.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax1.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats
    yTicks = list(set(yTicks))
    # ax1.set_yticks(yTicks)  # vecais scale variants
    # if old:
    # ax1.axvline(x=new_data[0][0], linewidth=4, color='r', linestyle='--', label="New Monitoring") # rakstam liniju neradam

    ax1.legend(fontsize=15)
    ax1.grid(False)
    ax2.grid(False)
    plt.xticks(rotation=0, ha='right')
    plt.savefig("/home/janis/Desktop/monitoring.eps", format="eps", dpi=5000)
    plt.show()

    result_calib = [x]
    if int(get_args("base")) >= 0:
        fig = plt.figure("Monitoring")
        ax1 = fig.add_subplot(111)

        for component in components:
            index = components.index(component)
            y = correctNumpyReadData(data[:, [index + 1]])
            base = correctNumpyReadData(data[:, [int(get_args("base"))]])
            ax1.plot(x, y / base, symbols[index] + colors[index], label=velocity[index], linewidth=0.5, markersize=5)
            result_calib.append(y / base)

        ax1.legend()
        ax1.grid(False)
        plt.show()

    result_org_file_name = get_configs("paths", "monitoringFilePath") + get_args("source") + ".out"
    result_calib_file_name = get_configs("paths", "monitoringFilePath") + get_args("source") + "_calib" + ".out"

    np.savetxt(result_org_file_name, np.transpose(result_org))
    np.savetxt(result_calib_file_name, np.transpose(result_calib))

    sys.exit(0)


if __name__ == "__main__":
    main()