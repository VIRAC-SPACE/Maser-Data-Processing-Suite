#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
from matplotlib.widgets import Slider, TextBox, Button
from pandas import DataFrame
from scipy import stats
from sympy import *
import re
from tabulate import tabulate

from parsers._configparser import ConfigParser
from help import *


def parse_arguments():
    parser = argparse.ArgumentParser(description='''Monitoring multiple sources. ''', epilog="""MMonitor.""")
    parser.add_argument("--sources", help="Sources Names", type=str, nargs='+', default="[g32p745, w51, g59p783, on1, s252, ngc7538, w3oh]")
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


def print_stats(stats):
    #DescribeResult(nobs=1, minmax=(1.2230937354864293, 1.2230937354864293), mean=1.2230937354864293, variance=nan, skewness=0.0, kurtosis=-3.0)
    headers = ["nobs", "min", "max", "mean", "variance", "skewness", "kurtosis"]
    data = [headers]
    for s in stats:
        stats = str(stats).replace("DescribeResult", "").replace("(", "").replace(")", "").replace(",", "").split(" ")
        stats = [re.sub("[^0-9.]", "", s.split("=")[-1])[0:5] for s in stats]
        data.append(stats)

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
        date = np.loadtxt(file, usecols=(0,), unpack=True)
        source = file.split("/")[-1].split(".")[0]
        lines_to_plot = velocities_to_plot_for_source[source]
        TMP = get_configs("velocities", source + "_6668" ).split(",")
        TMP = [v.strip().replace(" ", "") for v in TMP]
        idexie_for_lines_to_plot = []
        for tmpl in lines_to_plot:
            idexie_for_lines_to_plot.append(TMP.index(tmpl))

        lines[source]["date"] = date
        column_nr = len(get_configs("velocities", source + "_6668").split(","))

        for c in range(1, column_nr+1):
            if c - 1 in idexie_for_lines_to_plot:
                lines[source]["y_data"].append(np.loadtxt(file, usecols=(c,), unpack=True)/np.mean(np.loadtxt(file, usecols=(c,), unpack=True)))

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
    ax = gca()
    for source in sources:
        date = lines[source]["date"]
        slider_min.append(min(date))
        slider_max.append(max(date))
        velocities = get_configs("velocities", source + "_6668").split(",")
        velocities = [v.strip() for v in velocities]

        if source == "g32p745":
            velocities_tmp = ["30.49",  "39.18"]
        elif source == "w51":
            velocities_tmp = ["59.29"]
        elif source == "g59783":
            velocities_tmp = ["19.2"]
        elif source == "on1":
            velocities_tmp = ["14.64"]
        elif source == "s252":
            velocities_tmp = ["10.84"]
        elif source == "ngc7538":
            velocities_tmp = ["-58.04"]
        else:
            velocities_tmp = ["-44.6"]

        i = 0
        for y_data in lines[source]["y_data"]:
            z = np.polyfit(date, y_data, 2)
            p = np.poly1d(z)
            fit = str( z[0] ) + " * x1**2 "

            if "-" in str( z[1] ):
                fit += " + " + str( z[1] ) + " * x1 "
            else:
                fit += " + " + str( z[1] ) + " * x1"

            if "-" in str( z[2] ):
                fit += " + " + str( z[2] )
            else:
                fit += " + " + str( z[2] )
            line = plt.plot(date, y_data, "*", label=source+" velocity " + velocities_tmp[i] + " " + "y= " + fit)
            lines2.append(line)
            trend = plt.plot(date, p(date), "r--")
            trends.append(trend)
            i += 1

    ax.legend()
    fig2 = plt.figure()

    axcolor = 'lightgoldenrodyellow'
    a_axis = plt.axes([0.25, 0.6, 0.65, 0.03], facecolor=axcolor)
    b_axis = plt.axes([0.25, 0.5, 0.65, 0.03], facecolor=axcolor)
    a_slider = Slider(a_axis, 'a', max(slider_min), min(slider_max) -10, max(slider_min) + 10)
    b_slider = Slider(b_axis, 'b', max(slider_min) -10, min(slider_max), min(slider_max) - 10)
    a_slider.drawon = False
    b_slider.drawon = False

    axbox_l = plt.axes([0.25, 0.4, 0.65, 0.03])
    axbox_r = plt.axes([0.25, 0.3, 0.65, 0.03])
    text_box_l = TextBox(axbox_l, 'left', initial=str(max(slider_min)))
    text_box_r = TextBox(axbox_r, 'right', initial=str(min(slider_max)))
    axbox_b = plt.axes( [0.25, 0.2, 0.65, 0.03] )
    b = Button(axbox_b, "Compute statistical values", image=None, color='0.85', hovercolor='0.95')

    def compute_statistical_values(val):
        a = float(text_box_l.text)
        b = float(text_box_r.text)
        stats_list = []
        for t in range(0, len(trends)):
            x = lines2[t][0].get_xdata()
            y = lines2[t][0].get_ydata()
            i = findNearestIndex(x, a)
            j = findNearestIndex(x, b)
            y_tmp = y[i:j]
            label = lines2[t][0].get_label()
            stats_list.append(stats.describe(y_tmp))
        print_stats(stats_list)

    def update(val):
        a = a_slider.val
        b = b_slider.val

        xtmp = np.arange(a,b)
        data = {}
        columns = []

        for t in range(0, len(trends)):
            x = lines2[t][0].get_xdata()
            y = lines2[t][0].get_ydata()
            old_label = lines2[t][0].get_label()

            i = findNearestIndex(x, a)
            j = findNearestIndex(x, b)

            z = np.polyfit(x[i:j], y[i:j], 2)
            p = np.poly1d(z)

            fit = str(z[0]) + " * x1**2 "

            if "-" in str(z[1]):
                fit += " + " + str(z[1]) + " * x1 "
            else:
                fit += " + " + str(z[1]) + " * x1"

            if "-" in str(z[2]):
                fit += " + " + str(z[2])
            else:
                fit += " + " + str(z[2])

            equation = lambdify('x1', fit, 'numpy')
            ytmp = equation(xtmp)
            data["_".join(old_label.split(" ")[0:3])] = ytmp
            columns.append("_".join(old_label.split(" ")[0:3]))
            new_label = " ".join(old_label.split(" ")[0:3]) + " " + "y=" + fit
            lines2[t][0].set_label(new_label)
            trends[t][0].set_xdata(x[i:j])
            trends[t][0].set_ydata(p(x[i:j]))

        df = DataFrame(data, columns=columns)
        corr = df.corr()
        print(corr)
        ax.get_legend().remove()
        ax.legend()
        plt.draw()
        fig1.canvas.draw_idle()
        fig1.canvas.draw()
        fig1.canvas.flush_events()

    b.on_clicked(compute_statistical_values)
    a_slider.on_changed(update)
    b_slider.on_changed(update)

    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
