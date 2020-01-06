#! /usr/bin/python3
# -*- coding: utf-8 -*-
import matplotlib

matplotlib.use('Qt5Agg')

from matplotlib.ticker import StrMethodFormatter
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd
import argparse
import json
import datetime
import os

from monitor.months import Months
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
    parser.add_argument("-n", "--not_show", type=str, default="")
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

    old_dates = []
    old_log_file_dir = "/home/janis/Documents/maser/old_prety_log/"
    for old_log_file in os.listdir(old_log_file_dir):
        with open(old_log_file_dir + old_log_file, "r") as old_log_file2:
            lines = old_log_file2.readlines()
            station = lines[1].split(";")[1]
            if station == "IRBENE":
                line_index = -1
                for line in lines:
                    line_index += 1
                    if line.startswith("Source;"):
                        source = line.split(";")[1].split(",")[0]
                        if source == get_args("source"):
                            data = lines[line_index + 1].split(";")[1]
                            time = lines[line_index + 2].split(";")[1]
                            time_stamp = data + "_" + time
                            old_dates.append(time_stamp)

    old_dates = [datetime.datetime.strptime(d, '%Y-%m-%d_%H:%M:%S') for d in old_dates]

    component_count = len(get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(","))
    data_file = get_configs("paths", "monitoringFilePath") + get_args("source") + ".out"
    data = np.fromfile(data_file, dtype="float64", count=-1, sep=" ").reshape((file_len(data_file), component_count + 1))
    x = correctNumpyReadData(data[:, [0]])
    print("Number of observations", len(x))

    components = [i for i in range(1, component_count + 1)]
    symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', '#FF33E6']
    velocity = get_configs("velocities", get_args("source") + "_6668").replace(" ", "").split(",")

    # f, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 4]}, dpi=75, figsize=(19.8, 9.3))
    f, (ax1, ax2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [2, 4]},  figsize=(19.8, 9.3))
    f.tight_layout(pad=2, h_pad=2, w_pad=2, rect=None)
    f.subplots_adjust(wspace=0.1, right=0.876, bottom=0.076, top=0.945)
    f.suptitle(get_configs("Full_source_name", get_args("source")), horizontalalignment="center", verticalalignment="center", x=0.4)

    #ax1.set_aspect(aspect=1)
    #ax2.set_aspect(aspect=100000)
    
    ax1.set_aspect(aspect=0.5)
    ax2.set_aspect(aspect=0.3)

    #specter_files = ["w51_18_04_58_24_Jul_2018_IRBENE16_2.dat", "w51_18_43_13_11_Aug_2019_IRBENE_18.dat","w51_10_56_46_24_Nov_2019_IRBENE16_30.dat" ]
    #specter_files = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    #specter_files = ["s255_2017-04-21_m21_n25.dat.out", "s255_2018-03-05_m119_n41.dat.out", "s255_09_41_03_20_Oct_2019_IRBENE16_21.dat"]# S255
    #specter_files = ["w3oh_16_38_18_23_Mar_2019_IRBENE16_53.dat","w3oh_07_39_00_28_Nov_2018_IRBENE16_28.dat", "w3oh_20_22_47_08_Jun_2018_IRBENE16_9.dat"]# w3oh
    # specter_files = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]# g90p92

    # specter_files = ["w51_18_04_58_24_Jul_2018_IRBENE16_2.dat", "w51_18_43_13_11_Aug_2019_IRBENE_18.dat","w51_10_56_46_24_Nov_2019_IRBENE16_30.dat" ]
    # specter_files = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    # specter_files = ["s255_2017-04-21_m21_n25.dat.out", "s255_2018-03-05_m119_n41.dat.out", "s255_09_41_03_20_Oct_2019_IRBENE16_21.dat"]# S255
    # specter_files = ["w3oh_16_38_18_23_Mar_2019_IRBENE16_53.dat","w3oh_07_39_00_28_Nov_2018_IRBENE16_28.dat", "w3oh_20_22_47_08_Jun_2018_IRBENE16_9.dat"]# w3oh
    # specter_files = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]# g90p92
    specter_files = ["cepa_17_55_39_04_Oct_2019_IRBENE16_214.dat", "cepa_2018-06-04_m242_n29.dat.out", "cepa_07_05_58_10_Sep_2018_IRBENE16_144.dat"]
    # specter_files = ["g33p64_20_39_39_24_Jul_2018_IRBENE16_2.dat", "g33p64_00_26_34_23_Jun_2019_IRBENE_8.dat", "g33p64_13_20_17_14_Nov_2019_IRBENE16_33.dat"]
    # specter_files = ["g32p745_20_18_24_10_Jul_2018_IRBENE16_5.dat", "g32p745_21_18_07_15_Aug_2018_IRBENE16_18.dat", "g32p745_12_24_15_10_Dec_2018_IRBENE16_36.dat"]
    # specter_files = ["g25p65_20_38_40_07_Aug_2018_IRBENE16_8.dat","g25p65_00_23_27_25_May_2019_IRBENE_4.dat","g25p65_17_31_29_09_Oct_2019_IRBENE16_22.dat"]
    # specter_files = ["afgl5142_2017-03-20_m10_n15.dat.out","afgl5142_15_23_59_10_Jan_2019_IRBENE16_24.dat","afgl5142_05_50_02_27_Jul_2018_IRBENE16_4.dat"]
    # specter_files = ["g37p55_07_41_40_15_Mar_2019_IRBENE16_47.dat","g37p55_04_54_14_07_Jun_2019_IRBENE_6.dat","g37p55_11_36_00_22_Nov_2019_IRBENE16_31.dat"]
    # specter_files = ["g75p78_20_35_17_25_Aug_2018_IRBENE16_10.dat","g75p78_07_29_17_05_Jan_2019_IRBENE16_36.dat","g75p78_10_19_19_11_Nov_2019_IRBENE16_32.dat"]
    # specter_files = ["g78p12_2017-03-22_m11_n05.dat.out","g78p12_21_55_26_31_Oct_2019_IRBENE16_44.dat","g78p12_2018-07-12_m311_n25.dat.out"]
    # specter_files = ["ngc7538_05_53_01_27_Nov_2018_IRBENE16_15.dat","ngc7538_06_22_08_21_Dec_2018_IRBENE16_21.dat","ngc7538_15_21_34_24_Oct_2019_IRBENE16_35.dat"]
    # specter_files = ["ngc281_2017-03-27_m12_n12.dat.out","ngc281_10_59_24_20_Nov_2019_IRBENE16_40.dat","ngc281_07_18_55_12_Apr_2019_IRBENE16_51.dat"]
    # specter_files = ["s231_2017-04-24_m22_n26.dat.out","s231_16_04_42_10_Jan_2019_IRBENE16_26.dat","s231_10_13_42_20_Oct_2019_IRBENE16_22.dat"]
    # specter_files = ["w75n_18_26_02_08_Sep_2019_IRBENE16_25.dat","w75n_05_00_03_02_Oct_2018_IRBENE16_14.dat","w75n_12_17_07_18_Nov_2019_IRBENE16_35.dat"]
    # specter_files = ["w49n_20_19_51_03_Aug_2018_IRBENE16_5.dat","w49n_21_25_48_07_Oct_2019_IRBENE16_24.dat","w49n_05_52_43_09_Apr_2019_IRBENE16_50.dat"]
    # specter_files = ["v645_2018-07-05_m297_n13.dat.out","v645_06_15_15_15_Dec_2018_IRBENE16_125.dat","v645_09_40_25_07_Dec_2018_IRBENE16_123.dat"]
    # specter_files = ["s252_08_13_26_05_Oct_2019_IRBENE16_21.dat","s252_09_13_07_07_Jul_2018_IRBENE16_2.dat","s252_07_26_10_17_Sep_2018_IRBENE16_14.dat"]
    # specter_files = ["l1206_2017-04-13_m18_n11.dat.out","l1206_05_35_44_21_Aug_2018_IRBENE16_7.dat","l1206_17_59_05_15_Sep_2019_IRBENE16_26.dat"]
    # specter_files = ["l1287_11_15_36_20_Nov_2019_IRBENE16_37.dat","l1287_2017-05-05_m25_n15.dat.out","l1287_05_04_50_30_Sep_2018_IRBENE16_10.dat"]
    # specter_files = ["g37p479_19_44_36_31_Jul_2018_IRBENE16_4.dat","g37p479_12_14_00_17_Nov_2019_IRBENE16_29.dat","g37p479_07_00_54_15_Mar_2019_IRBENE16_47.dat"]
    # specter_files = ["g107p3_12_59_26_04_Apr_2019_IRBENE16_938.dat","g107p3_10_47_23_27_Mar_2019_IRBENE16_919.dat","g107p3_06_52_31_01_Apr_2019_IRBENE16_927.dat"]
    # specter_files = ["afgl6366_05_47_10_22_Aug_2018_IRBENE16_7.dat","afgl6366_18_43_29_18_Dec_2018_IRBENE16_21.dat","afgl6366_13_45_13_02_Apr_2019_IRBENE16_33.dat"]
    # specter_files = ["g35p20_2017-08-01_m56_n36.dat.out","g35p20_05_58_24_27_Dec_2018_IRBENE16_32.dat","g35p20_18_45_16_19_Aug_2018_IRBENE16_9.dat"]
    # specter_files = ["w3oh_05_41_38_29_Aug_2018_IRBENE16_21.dat","w3oh_17_16_16_13_Oct_2019_IRBENE16_31.dat","w3oh_21_08_41_07_Oct_2019_IRBENE16_30.dat"]
    # specter_files = ["g43p79_18_05_19_12_Aug_2018_IRBENE16_6.dat","g43p79_06_59_44_08_Dec_2018_IRBENE16_23.dat","g43p79_16_28_40_06_Nov_2019_IRBENE16_30.dat"]
    # specter_files = ["on1_06_10_52_01_Apr_2019_IRBENE16_55.dat","on1_18_38_21_14_Aug_2018_IRBENE16_7.dat","on1_11_35_30_19_Nov_2019_IRBENE16_33.dat"]
    # specter_files = ["g34p403_12_39_51_19_Nov_2019_IRBENE16_34.dat","g34p403_19_44_01_29_Jul_2018_IRBENE16_2.dat","g34p403_06_39_07_27_Dec_2018_IRBENE16_27.dat"]
    #specter_files = ["g50p03_2017-06-24_m43_n29.dat.out", "g50p03_2018-06-11_m257_n30.dat.out", "g50p03_2017-09-28_m73_n29.dat.out"]
    #specter_files = ["g59p783_2017-03-20_m10_n04.dat.out", "g59p783_2017-11-13_m89_n29.dat.out", "g59p783_2017-05-18_m29_n09.dat.out"]
    #specter_files = ["g60p57_2017-05-08_m26_n09.dat.out", "g60p57_2017-04-13_m18_n06.dat.out", "g60p57_2017-03-22_m11_n03.dat.out"]
    #specter_files = ["dr20_2017-04-21_m21_n17.dat.out", "dr20_2018-03-18_m122_n27.dat.out", "dr20_2017-11-10_m88_n27.dat.out"]
    #specter_files = ["dr21_2017-07-25_m54_n16.dat.out", "dr21_2017-05-05_m25_n20.dat.out", "dr21_2017-04-03_m15_n11.dat.out"]
    #specter_files = ["g136p84_2018-06-01_m239_n30.dat.out", "g136p84_2017-06-27_m44_n20.dat.out", "g136p84_2017-04-21_m21_n11.dat.out"]
    #specter_files = ["g183p35_2018-01-25_m107_n37.dat.out", "g183p35_2017-09-13_m68_n08.dat.out", "g183p35_2017-05-05_m25_n12.dat.out"]
    #specter_files = ["g188p79_2018-05-30_m235_n34.dat.out", "g188p79_2017-07-08_m47_n09.dat.out", "g188p79_2017-05-24_m32_n19.dat.out"]
    #specter_files = ["g30p99_19_14_48_04_Sep_2019_IRBENE16_19.dat", "g30p99_09_39_23_08_Jan_2019_IRBENE16_32.dat", "g30p99_12_07_53_19_Nov_2019_IRBENE16_34.dat"]
    specter_files = ["g30p99_19_14_48_04_Sep_2019_IRBENE16_19.dat", "g30p99_09_39_23_08_Jan_2019_IRBENE16_32.dat", "g30p99_12_07_53_19_Nov_2019_IRBENE16_34.dat"]

    symbols_2 = [".", "-.", "-"]
    i = 0

    for specter in specter_files:
        file = get_configs("paths", "outputFilePath") + "/" + get_args("source") + "/6668/" + specter
        spectre_data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file), 4))

        x_ = correctNumpyReadData(spectre_data[:, [0]])
        y_ = correctNumpyReadData(spectre_data[:, [3]])

        ax1.plot(x_, y_, symbols_2[i])
        i += 1

    ax1.set_xlim(52.0, 65.0)
    ax1.set_xlabel("Velocity (km sec$^{-1}$)")
    ax1.set_ylabel("Flux density (Jy)")
    # ax1.set_ylim(-2,9)

    ax2.plot([], [], ' ', label="km sec$^{-1}$")

    print(get_configs("paths", "resultFilePath") + get_args("source") + "_6668.json")
    with open(get_configs("paths", "resultFilePath") + "/" + get_args("source") + "_6668.json", "r") as r_data:
        result_data = json.load(r_data)

    months = Months()
    RT32_observation_dates = old_dates

    for observation in result_data:
        if observation.split("_")[-2] == "IRBENE":
            date = result_data[observation]["Date"]
            dates = date.split("_")
            monthsNumber = dates[1]
            dates[1] = months.getMonthNumber([monthsNumber][0])
            date = result_data[observation]["time"].replace(":", " ") + " " + " ".join(dates)
            date = datetime.datetime.strptime(date, '%H %M %S %d %m %Y')
            date = convertDatetimeObjectToMJD(date)
            RT32_observation_dates.append(date)

    rt32_x = [rt32 for rt32 in x if rt32 in RT32_observation_dates]
    rt32_x_indexies = [x.tolist().index(rt32) for rt32 in x if rt32 in RT32_observation_dates]

    markersize = [10 for i in range(0, len(x))]
    for component in components:
        index = components.index(component)

        if len(get_args("not_show")) != 0:

            if index + 1 != 3:
                y = correctNumpyReadData(data[:, [index + 1]])
                ax2.plot(x, y, symbols[index], color=colors[index], linewidth=0.5, markersize=5,label=str(velocity[index]))
                ax2.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar

        else:
            y = correctNumpyReadData(data[:, [index + 1]])
            rt32_y = y[rt32_x_indexies]
            ax2.plot(x, y, symbols[index], c=colors[index], linewidth=0.5, markersize=4, label=str(velocity[index]))
            ax2.plot(rt32_x, rt32_y, symbols[index], c=colors[index], linewidth=0.5, markersize=10, label=str(velocity[index]))
            ax2.errorbar(x[0], y[0], yerr=1.5 + 0.05 * y[0], xerr=None, ls='none', ecolor='k')  # 1st poiont error bar

    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=12, rotation=0)  # MJD atzimju formâts
    ax1.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax1.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats
    # ax2.set_ylim(1)
    ax2.set_yscale("log")
    ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))  # clasic log scale!

    #ax2.yaxis.set_minor_formatter(StrMethodFormatter('{x:.0f}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter('{x:.0f}'))
    

    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.tick_params(axis="x", direction="in", which="both", length=16, width=2, labelsize=13, rotation=0)  # MJD atzimju formâts
    ax2.tick_params(axis="y", direction="in", which="major", length=16, width=2, labelsize=12)  # Flux atzimju formats
    ax2.tick_params(axis="y", direction="in", which="minor", length=10, width=1.5)  # Flux atzimju formats
    ax2.set_xlabel("MJD")

    ax2.legend(bbox_to_anchor=(1, 0.5, 0.3, 0.3), loc='upper left', borderaxespad=0.5)

    plt.show()


    #f.savefig("/home/artis/monitoring_g34p403.eps", format="eps", dpi=5000)  # , papertype="a4")
    f.savefig("monitoring_" + get_args("source") + ".eps", format="eps", dpi=5000, papertype="a4")



if __name__ == "__main__":
    main()