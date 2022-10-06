#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Display monitoring for multiple sources
"""
import sys
import argparse
import re
import json

from scipy.ndimage import uniform_filter1d, gaussian_filter1d
from sympy import lambdify
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import gca
from matplotlib.widgets import Slider, TextBox, Button
from pandas import DataFrame
from scipy import stats

from parsers.configparser_ import ConfigParser
from utils.help import find_nearest_index


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring multiple sources. ''')
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("--sources", help="Sources Names", type=str, nargs='+',
                        default="[g32p745, w51, g59p783, on1, s252, ngc7538, w3oh, w49n, g85p41, g78p12, g75p78]")
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


def get_velocities_tmp(source):
    """

    :param source: observed source
    :return: stable velocity for source
    """
    if source == "g32p745":
        velocities_tmp = ["30.49", "39.18"]
    elif source == "w51":
        velocities_tmp = ["59.29"]
    elif source == "g59p783":
        velocities_tmp = ["27.15"]
    elif source == "on1":
        velocities_tmp = ["14.64"]
    elif source == "s252":
        velocities_tmp = ["10.84"]
    elif source == "ngc7538":
        velocities_tmp = ["-61.31"]
    elif source == "w49n":
        velocities_tmp = ["9.27"]
    elif source == "g85p41":
        velocities_tmp = ["-31.65"]
    elif source == "g78p12":
        velocities_tmp = ["-6.13"]
    elif source == "g75p78":
        velocities_tmp = ["-2.57"]
    else:
        velocities_tmp = ["-44.6"]
    return velocities_tmp


def get_time_cut(source):
    """

    :param source: observed source
    :return: stable time cut  for source
    """
    if source == "g32p745":
        time_cut = [58400, 59800]
    elif source == "w51":
        time_cut = [58400, 59800]
    elif source == "g59p783":
        time_cut = [58400, 59800]
    elif source == "on1":
        time_cut = [58400, 59800]
    elif source == "s252":
        time_cut = [58400, 59800]
    elif source == "ngc7538":
        time_cut = [58400, 59800]
    elif source == "w49n":
        time_cut = [58400, 59800]
    elif source == "g85p41":
        time_cut = [58400, 59800]
    elif source == "g78p12":
        time_cut = [58400, 59800]
    elif source == "g75p78":
        time_cut = [58400, 59800]
    else:
        time_cut = [58400, 59800]
    return time_cut


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
    """

    :param monitoring_files: all monitoring files
    :param sources: all stable sources
    :return:
    """
    lines = {}

    velocities_to_plot_for_source = {}
    for source in sources:
        lines[source] = {"y_data": []}
        velocities_to_plot_for_source[source] = get_velocities_tmp(source)

    for file in monitoring_files:
        data = np.load(file, allow_pickle=True)
        date = data[0][0]
        source = file.split("/")[-1].split(".")[0].split("_")[0]
        lines_to_plot = velocities_to_plot_for_source[source]
        tmp = get_configs("velocities", source + "_" + get_args("line")).split(",")
        tmp = [v.strip().replace(" ", "") for v in tmp]
        idexie_for_lines_to_plot = []
        for tmpl in lines_to_plot:
            idexie_for_lines_to_plot.append(tmp.index(tmpl))

        lines[source]["date"] = date
        column_nr = len(get_configs("velocities", source + "_" + get_args("line")).split(","))

        for c in range(1, column_nr+1):
            if c - 1 in idexie_for_lines_to_plot:
                tmp = data[c] / np.mean(data[c])
                lines[source]["y_data"].append(tmp)

    return lines


def get_iterations_from_mjd(star_time, stop_time):
    iterations = dict()
    source_list = get_args("sources").replace("[", "").replace("]", "").replace("'", "").split(",")
    result_file_path = get_configs("paths", "resultFilePath")

    for source in source_list:
        result_file = result_file_path + source + "_" + get_args("line") + ".json"

        with open(result_file) as result_data:
            result = json.load(result_data)
            modified_julian_days = []
            iteration_numbers = []
            for observation in result:
                modified_julian_days.append(result[observation]["modifiedJulianDays"])
                iteration_numbers.append(result[observation]["Iteration_number"])

            if len(modified_julian_days) > 0:
                left_index = find_nearest_index(modified_julian_days, star_time)
                right_index = find_nearest_index(modified_julian_days, stop_time)
                iterations[source] = iteration_numbers[min(right_index,left_index):max(right_index,left_index)]

    return iterations


def main():
    sources = get_args("sources").replace("[", "").replace("]", "").replace("'", "").split(",")
    sources = [s.strip() for s in sources]
    monitoring_files = []
    monitoring_dir = get_configs("paths", "monitoringFilePath")

    for source in sources:
        monitoring_files.append(monitoring_dir + source + "_" + get_args("line") + ".npy")

    lines = read_monitoring_files(monitoring_files, sources)

    slider_min = []
    slider_max = []
    trends = []
    lines2 = []
    fig1 = plt.figure()
    ax = gca()

    out = np.array([[],[]])

    for source in sources:
        cut_index = get_time_cut(source)
        date = lines[source]["date"]
        slider_min.append(min(date))
        slider_max.append(max(date))
        velocities_tmp = get_velocities_tmp(source)

        i = 0
        for y_data in lines[source]["y_data"]:
            z = np.polyfit(date, y_data, 2)
            p = np.poly1d(z)
            fit = str(z[0]) + " * x**2 "

            if "-" in str(z[1]):
                fit += " + " + str(z[1]) + " * x "
            else:
                fit += " + " + str(z[1]) + " * x "

            if "-" in str(z[2]):
                fit += " + " + str(z[2])
            else:
                fit += " + " + str(z[2])
            line = plt.plot(date, y_data, "*", label=source + " velocity " + velocities_tmp[i] + " " + "y= " + fit)
            lines2.append(line)
            trend = plt.plot(date, p(date), "r--", visible=False)
            trends.append(trend)

            start_cut = find_nearest_index(date, cut_index[0])
            end_cut = find_nearest_index(date, cut_index[1])
            out = np.concatenate((out, [date[start_cut:end_cut],
                                        y_data[start_cut:end_cut]/np.mean(y_data[start_cut:end_cut])]), axis=1)
            i += 1

    dtype = [('mjd', float), ('amp', float)]
    values = []

    for i in range(0, out.shape[1]):
        values.append((out[0][i], out[1][i]))

    out = np.array(values, dtype=dtype)
    out = np.sort(out, order='mjd')
    out_filtered = gaussian_filter1d(out['amp'], 10, mode='nearest')

    fig1, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), dpi=90)
    ax1.scatter(out['mjd'], out['amp'])
    ax1.plot(out['mjd'], out_filtered, color='red')

    fig2, ax2 = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), dpi=90)
    out_filtered_2 = (out['amp'] + np.max(out_filtered))/out_filtered
    ax2.scatter(out['mjd'], out_filtered_2, color='yellow')

    ax.legend()
    fig2 = plt.figure()

    axcolor = 'lightgoldenrodyellow'
    a_axis = plt.axes([0.25, 0.6, 0.65, 0.03], facecolor=axcolor)
    b_axis = plt.axes([0.25, 0.5, 0.65, 0.03], facecolor=axcolor)
    a_slider = Slider(a_axis, 'a', max(slider_min), min(slider_max) - 10, max(slider_min) + 10)
    b_slider = Slider(b_axis, 'b', max(slider_min) - 10, min(slider_max), min(slider_max) - 10)
    a_slider.drawon = False
    b_slider.drawon = False

    axbox_l = plt.axes([0.25, 0.4, 0.65, 0.03])
    axbox_r = plt.axes([0.25, 0.3, 0.65, 0.03])
    text_box_l = TextBox(axbox_l, 'left', initial=str( max( slider_min)))
    text_box_r = TextBox(axbox_r, 'right', initial=str( min( slider_max)))
    axbox_b = plt.axes([0.25, 0.2, 0.65, 0.03])
    b = Button(axbox_b, "Compute statistical values", image=None, color='0.85', hovercolor='0.95')

    def compute_statistical_values(val):
        a = float(text_box_l.text)
        b = float(text_box_r.text)
        stats_list = []
        labels = []
        for t in range(0, len(trends)):
            x = lines2[t][0].get_xdata()
            y = lines2[t][0].get_ydata()
            i = find_nearest_index(x, a)
            j = find_nearest_index(x, b)
            y_tmp = y[i:j]
            label = lines2[t][0].get_label()
            labels.append(label)
            stats_list.append(stats.describe(y_tmp))
        print_stats(stats_list, labels)

    def update(val):
        a = a_slider.val
        b = b_slider.val

        print("Iteration list for time range ", a, b, get_iterations_from_mjd(a, b))

        xtmp = np.arange(a, b)
        data = {}
        columns = []

        for t in range(0, len(trends)):
            x = lines2[t][0].get_xdata()
            y = lines2[t][0].get_ydata()
            old_label = lines2[t][0].get_label()

            i = find_nearest_index(x, a)
            j = find_nearest_index(x, b)

            z = np.polyfit(x[i:j], y[i:j], 2)
            p = np.poly1d(z)

            fit = str(z[0]) + " * x**2 "

            if "-" in str(z[1]):
                fit += " + " + str(z[1]) + " * x "
            else:
                fit += " + " + str(z[1]) + " * x"

            if "-" in str(z[2]):
                fit += " + " + str(z[2])
            else:
                fit += " + " + str(z[2])

            equation = lambdify('x', fit, 'numpy')
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
