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
    specter_files_S255 = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    specter_files_w3oh = ["s255_2017-04-21_m21_n25.dat.out", "s255_2018-03-05_m119_n41.dat.out", "s255_09_41_03_20_Oct_2019_IRBENE16_21.dat"]
    specter_files_g90p92 = ["w3oh_16_38_18_23_Mar_2019_IRBENE16_53.dat","w3oh_07_39_00_28_Nov_2018_IRBENE16_28.dat", "w3oh_20_22_47_08_Jun_2018_IRBENE16_9.dat"]
    specter_files_g90p92 = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    specter_files_w51 = ["w51_18_04_58_24_Jul_2018_IRBENE16_2.dat", "w51_18_43_13_11_Aug_2019_IRBENE_18.dat","w51_10_56_46_24_Nov_2019_IRBENE16_30.dat" ]
    specter_files_g90p92 = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    specter_files_S255 = ["s255_2017-04-21_m21_n25.dat.out", "s255_2018-03-05_m119_n41.dat.out", "s255_09_41_03_20_Oct_2019_IRBENE16_21.dat"]
    specter_files_w3oh = ["w3oh_16_38_18_23_Mar_2019_IRBENE16_53.dat","w3oh_07_39_00_28_Nov_2018_IRBENE16_28.dat", "w3oh_20_22_47_08_Jun_2018_IRBENE16_9.dat"]
    specter_files_g90p92 = ["g90p92_15_58_46_23_Mar_2019_IRBENE16_42.dat","g90p92_09_09_22_18_Jan_2019_IRBENE16_27.dat", "g90p92_11_10_44_23_Nov_2019_IRBENE16_45.dat"]
    specter_files_cepa = ["cepa_17_55_39_04_Oct_2019_IRBENE16_214.dat", "cepa_2018-06-04_m242_n29.dat.out", "cepa_07_05_58_10_Sep_2018_IRBENE16_144.dat"]
    specter_files_g33p64 = ["g33p64_20_39_39_24_Jul_2018_IRBENE16_2.dat", "g33p64_00_26_34_23_Jun_2019_IRBENE_8.dat", "g33p64_13_20_17_14_Nov_2019_IRBENE16_33.dat"]
    specter_files_g32p745 = ["g32p745_20_18_24_10_Jul_2018_IRBENE16_5.dat", "g32p745_21_18_07_15_Aug_2018_IRBENE16_18.dat", "g32p745_12_24_15_10_Dec_2018_IRBENE16_36.dat"]
    specter_files_g25p65 = ["g25p65_20_38_40_07_Aug_2018_IRBENE16_8.dat","g25p65_00_23_27_25_May_2019_IRBENE_4.dat","g25p65_17_31_29_09_Oct_2019_IRBENE16_22.dat"]
    specter_files_afgl5142 = ["afgl5142_2017-03-20_m10_n15.dat.out","afgl5142_15_23_59_10_Jan_2019_IRBENE16_24.dat","afgl5142_05_50_02_27_Jul_2018_IRBENE16_4.dat"]
    specter_files_g37p55 = ["g37p55_07_41_40_15_Mar_2019_IRBENE16_47.dat","g37p55_04_54_14_07_Jun_2019_IRBENE_6.dat","g37p55_11_36_00_22_Nov_2019_IRBENE16_31.dat"]
    specter_files_g75p78 = ["g75p78_20_35_17_25_Aug_2018_IRBENE16_10.dat","g75p78_07_29_17_05_Jan_2019_IRBENE16_36.dat","g75p78_10_19_19_11_Nov_2019_IRBENE16_32.dat"]
    specter_files_g78p122 = ["g78p122_2017-03-22_m11_n05.dat.out","g78p12_21_55_26_31_Oct_2019_IRBENE16_44.dat","g78p122_2018-07-12_m311_n25.dat.out"]
    specter_files_ngc7538 = ["ngc7538_05_53_01_27_Nov_2018_IRBENE16_15.dat","ngc7538_06_22_08_21_Dec_2018_IRBENE16_21.dat","ngc7538_15_21_34_24_Oct_2019_IRBENE16_35.dat"]
    specter_files_ngc281 = ["ngc281_2017-03-27_m12_n12.dat.out","ngc281_10_59_24_20_Nov_2019_IRBENE16_40.dat","ngc281_07_18_55_12_Apr_2019_IRBENE16_51.dat"]
    specter_files_s231 = ["s231_2017-04-24_m22_n26.dat.out","s231_16_04_42_10_Jan_2019_IRBENE16_26.dat","s231_10_13_42_20_Oct_2019_IRBENE16_22.dat"]
    specter_files_w75n = ["w75n_18_26_02_08_Sep_2019_IRBENE16_25.dat","w75n_05_00_03_02_Oct_2018_IRBENE16_14.dat","w75n_12_17_07_18_Nov_2019_IRBENE16_35.dat"]
    specter_files_w49n = ["w49n_20_19_51_03_Aug_2018_IRBENE16_5.dat","w49n_21_25_48_07_Oct_2019_IRBENE16_24.dat","w49n_05_52_43_09_Apr_2019_IRBENE16_50.dat"]
    specter_files_v645 = ["v645_2018-07-05_m297_n13.dat.out","v645_06_15_15_15_Dec_2018_IRBENE16_125.dat","v645_09_40_25_07_Dec_2018_IRBENE16_123.dat"]
    specter_files_s252 = ["s252_08_13_26_05_Oct_2019_IRBENE16_21.dat","s252_09_13_07_07_Jul_2018_IRBENE16_2.dat","s252_07_26_10_17_Sep_2018_IRBENE16_14.dat"]
    specter_files_l1206 = ["l1206_2017-04-13_m18_n11.dat.out","l1206_05_35_44_21_Aug_2018_IRBENE16_7.dat","l1206_17_59_05_15_Sep_2019_IRBENE16_26.dat"]
    specter_files_l1287 = ["l1287_11_15_36_20_Nov_2019_IRBENE16_37.dat","l1287_2017-05-05_m25_n15.dat.out","l1287_05_04_50_30_Sep_2018_IRBENE16_10.dat"]
    specter_files_g37p479 = ["g37p479_19_44_36_31_Jul_2018_IRBENE16_4.dat","g37p479_12_14_00_17_Nov_2019_IRBENE16_29.dat","g37p479_07_00_54_15_Mar_2019_IRBENE16_47.dat"]
    specter_files_g107p3 = ["g107p3_12_59_26_04_Apr_2019_IRBENE16_938.dat","g107p3_10_47_23_27_Mar_2019_IRBENE16_919.dat","g107p3_06_52_31_01_Apr_2019_IRBENE16_927.dat"]
    specter_files_afgl6366 = ["afgl6366_05_47_10_22_Aug_2018_IRBENE16_7.dat","afgl6366_18_43_29_18_Dec_2018_IRBENE16_21.dat","afgl6366_13_45_13_02_Apr_2019_IRBENE16_33.dat"]
    specter_files_g35p20 = ["g35p20_2017-08-01_m56_n36.dat.out","g35p20_05_58_24_27_Dec_2018_IRBENE16_32.dat","g35p20_18_45_16_19_Aug_2018_IRBENE16_9.dat"]
    specter_files_w3oh = ["w3oh_05_41_38_29_Aug_2018_IRBENE16_21.dat","w3oh_17_16_16_13_Oct_2019_IRBENE16_31.dat","w3oh_21_08_41_07_Oct_2019_IRBENE16_30.dat"]
    specter_files_g43p79 = ["g43p79_18_05_19_12_Aug_2018_IRBENE16_6.dat","g43p79_06_59_44_08_Dec_2018_IRBENE16_23.dat","g43p79_16_28_40_06_Nov_2019_IRBENE16_30.dat"]
    specter_files_on1 = ["on1_06_10_52_01_Apr_2019_IRBENE16_55.dat","on1_18_38_21_14_Aug_2018_IRBENE16_7.dat","on1_11_35_30_19_Nov_2019_IRBENE16_33.dat"]
    specter_files_g34p403 = ["g34p403_12_39_51_19_Nov_2019_IRBENE16_34.dat","g34p403_19_44_01_29_Jul_2018_IRBENE16_2.dat","g34p403_06_39_07_27_Dec_2018_IRBENE16_27.dat"]
    specter_files_g34p403 = ["g34p403_12_39_51_19_Nov_2019_IRBENE16_34.dat", "g34p403_19_44_01_29_Jul_2018_IRBENE16_2.dat", "g34p403_06_39_07_27_Dec_2018_IRBENE16_27.dat"]

    specter_files_for_all_sourses = [specter_files_g37p55, specter_files_S255, specter_files_ngc281, specter_files_g33p64, specter_files_l1287, specter_files_g78p122, specter_files_v645, specter_files_on1, specter_files_s252, specter_files_w51, specter_files_g25p65, specter_files_w75n, specter_files_g90p92,specter_files_g90p92, specter_files_g107p3, specter_files_afgl6366, specter_files_afgl5142, specter_files_w49n, specter_files_g90p92, specter_files_ngc7538, specter_files_g37p479, specter_files_s255, specter_files_g34p403, specter_files_s231, specter_files_l1206, specter_files_g35p20, specter_files_g75p78, specter_files_w3oh, specter_files_cepa, specter_files_g32p745, specter_files_g43p79]

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