#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from matplotlib.widgets import Slider

from parsers._configparser import ConfigParser
from help import *


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

    slider_min = []
    slider_max = []
    trends = []
    lines2 = []
    fig1 = plt.figure()
    for source in sources:
        date = lines[source]["date"]
        slider_min.append(min(date))
        slider_max.append(max(date))
        velocities = get_configs("velocities", source).split(",")
        velocities = [v.strip() for v in velocities]

        i = 0
        for y_data in lines[source]["y_data"]:
            z = np.polyfit(date, y_data, 1)
            p = np.poly1d(z)
            line = plt.plot(date, y_data, label=source+" velocity " + velocities[i] + " " + "y=%.6fx+(%.6f)"%(z[0],z[1]) + " " + "R^2 = " + str(r2_score(y_data, p(date))))
            lines2.append(line)
            trend = plt.plot(date, p(date), "r--")
            trends.append(trend)
            i += 1
    plt.legend()

    plt.figure()

    axcolor = 'lightgoldenrodyellow'
    a_axis = plt.axes([0.25, 0.6, 0.65, 0.03], facecolor=axcolor )
    b_axis = plt.axes([0.25, 0.5, 0.65, 0.03], facecolor=axcolor )
    a_slider = Slider(a_axis, 'a', max(slider_min), min(slider_max) - 10, max(slider_min) + 10)
    b_slider = Slider(b_axis, 'b', max(slider_min) + 10, min(slider_max) - 10, min(slider_max) - 20)

    def update(val):
        a = a_slider.val
        b = b_slider.val

        for t in range(0, len(trends)):
            x = lines2[t][0].get_xdata()
            y = lines2[t][0].get_ydata()

            i = findNearestIndex(x, a)
            j = findNearestIndex(x, b)

            z = np.polyfit(x[i:j], y[i:j], 1)
            p = np.poly1d(z)
            trends[t][0].set_xdata(x[i:j])
            trends[t][0].set_ydata(p(x[i:j]))

        fig1.canvas.draw_idle()

    a_slider.on_changed(update)
    b_slider.on_changed(update)

    plt.show()

    sys.exit(0)


if __name__ == "__main__":
    main()
