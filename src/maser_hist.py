#! /usr/bin/python3
# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Qt5Agg')

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import numpy as np
import argparse
import sys

from help import *
from parsers._configparser import ConfigParser


def get_configs_items():
    config_file_path = "config/plot.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(config_file_path)
    return config.getItems("main")


def parse_arguments():
    parser = argparse.ArgumentParser(description='''Compute flux density. ''', epilog="""density.""")
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
    fontP = FontProperties()
    fontP.set_size('small')

    config_items = get_configs_items()
    for key, value in config_items.items():
        rcParams[key] = value

    # Pirmais fails ir ar mazako amplitudu, otrais ir ar videjo amplitudu un tresai ir ar augstako amplitudu
    specter_files_w51 = ["w51_18_04_58_24_Jul_2018_IRBENE16_2.dat", "w51_18_43_13_11_Aug_2019_IRBENE_18.dat","w51_10_56_46_24_Nov_2019_IRBENE16_30.dat" ]
    specter_files_s255 = ["s255_2017-04-21_m21_n25.dat.out", "s255_2018-03-05_m119_n41.dat.out", "s255_09_41_03_20_Oct_2019_IRBENE16_21.dat"]
    specter_files_w3oh = ["w3oh_16_38_18_23_Mar_2019_IRBENE16_53.dat","w3oh_07_39_00_28_Nov_2018_IRBENE16_28.dat", "w3oh_20_22_47_08_Jun_2018_IRBENE16_9.dat"]
    specter_files_g90p92 = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat", "g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat","g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]

    specter_files_for_all_sourses = [specter_files_w51, specter_files_s255, specter_files_w3oh, specter_files_g90p92]

    density_low = []
    density_middle = []
    density_high = []
    file_name_index = 0

    for source in specter_files_for_all_sourses:
        source_name = source[0].split("_")[0]
        for file_name in source:
            file = get_configs("paths", "outputFilePath") + "/" + source_name + "/6668/" + file_name
            spectre_data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file), 4))

            x_ = correctNumpyReadData(spectre_data[:, [0]])
            y_ = correctNumpyReadData(spectre_data[:, [3]])

            area = np.trapz(y_, x_)
            print(area)

            if file_name_index == 0:
                density_low.append(area)
                file_name_index += 1

            elif file_name_index == 1:
                density_middle.append(area)
                file_name_index += 1

            elif file_name_index == 2:
                density_high.append(area)
                file_name_index = 0

    plt.figure("low")
    plt.hist(density_low)

    plt.figure("midle")
    plt.hist(density_middle)

    plt.figure("high")
    plt.hist(density_high)

    plt.show()
    sys.exit(0)


if __name__ == main():
    main()