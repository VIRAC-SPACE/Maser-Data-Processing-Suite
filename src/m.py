#! /usr/bin/python3
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Qt5Agg')

from matplotlib.ticker import StrMethodFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import numpy as np
import argparse

from help import *
from parsers._configparser import ConfigParser


def get_configs_items():
    config_file_path = "config/plot.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(config_file_path)
    return config.getItems("main")


def parse_arguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("-b", "--base", help="Base component", type=int, default=-1)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-n", "--not_show", type=str,  default="")
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
    fontP = FontProperties()
    fontP.set_size('small')

    config_items = get_configs_items()
    for key, value in config_items.items():
        rcParams[key] = value

    component_count = len(get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(","))
    data_file = get_configs("paths", "monitoringFilePath") + get_args("source") + ".out"
    data = np.fromfile(data_file, dtype="float64", count=-1, sep=" ").reshape((file_len(data_file), component_count + 1))
    x = correctNumpyReadData(data[:, [0]])
    print("Number of observations", len(x))

    components = [i for i in range(1, component_count + 1)]
    symbols = ["*-", "o-", "v-", "^-", "<-", ">-", "1-", "2-", "3-", "4-"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    velocity = get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(",")

    f, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 4]}, dpi=75, figsize=(20.7, 9.3))
    f.tight_layout(pad=2, h_pad=2, w_pad=2, rect=None)
    f.subplots_adjust(wspace=0.1, right=0.876, bottom=0.076, top=0.945)
    f.suptitle(get_configs("Full_source_name", get_args("source")), horizontalalignment="center", verticalalignment="center", x=0.4)

    specter_files = ["w51_18_04_58_24_Jul_2018_IRBENE16_2.dat", "w51_18_43_13_11_Aug_2019_IRBENE_18.dat","w51_10_56_46_24_Nov_2019_IRBENE16_30.dat" ]
    #specter_files = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    #specter_files = ["s255_2017-04-21_m21_n25.dat.out", "s255_2018-03-05_m119_n41.dat.out", "s255_09_41_03_20_Oct_2019_IRBENE16_21.dat"]# S255
    #specter_files = ["w3oh_16_38_18_23_Mar_2019_IRBENE16_53.dat","w3oh_07_39_00_28_Nov_2018_IRBENE16_28.dat", "w3oh_20_22_47_08_Jun_2018_IRBENE16_9.dat"]# w3oh
    #specter_files = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]# g90p92

    symbols_2 = [".", "-.", "-"]
    i = 0

    for specter in specter_files:
        file = get_configs("paths", "outputFilePath") + "/" + get_args("source") + "/6668/" + specter
        spectre_data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file), 4))

        x_ = correctNumpyReadData(spectre_data[:, [0]])
        y_ = correctNumpyReadData(spectre_data[:, [3]])

        ax1.plot(x_, y_, symbols_2[i])
        i += 1

    ax1.set_xlim(46.0, 65.0)
    ax1.set_xlabel("Velocity (km sec$^{-1}$)")
    ax1.set_ylabel("Flux density (Jy)")

    ax2.plot([], [], ' ', label="km sec$^{-1}$")
    for component in components:
        index = components.index(component)

        if len(get_args("not_show")) !=0:

            if index + 1 != 3:
                y = correctNumpyReadData(data[:, [index + 1]])
                ax2.plot(x, y, symbols[index] + colors[index], linewidth=0.5, markersize=5, label=str(velocity[index]))
                ax2.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar

        else:
            y = correctNumpyReadData(data[:, [index + 1]])
            ax2.plot(x, y, symbols[index] + colors[index], linewidth=0.5, markersize=5, label=str(velocity[index]))
            ax2.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar

    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=12,rotation=0)  # MJD atzimju formâts
    ax1.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax1.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats

    ax2.set_yscale("log")
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))  # clasic log scale!
    #ax2.yaxis.set_minor_formatter(StrMethodFormatter('{x:.0f}'))
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=13,rotation=0)  # MJD atzimju formâts
    ax2.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax2.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats
    ax2.set_xlabel("MJD")

    ax2.legend(bbox_to_anchor=(1, 0.5, 0.3, 0.3), loc='upper left', borderaxespad=0.5)

    plt.show()
    f.savefig("/home/artis/monitoring_w51.eps", format="eps", dpi=5000, papertype="a4")


if __name__ == "__main__":
    main()