#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
For definite (low, middle, high) output files for give sources compute spectral density and display histogram
"""
import sys
import os
import argparse
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
from matplotlib import rc
import numpy as np
import h5py

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''For definite (low, middle, high) output files for give sources 
    compute spectral density and display histogram. ''')
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="../config/config.cfg")
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


def get_configs_items():
    """

    :return: None
    """
    config_file_path = "../config/plot.cfg"
    config = ConfigParser(config_file_path)
    return config.get_items("main")


def main():
    rc('font', family='serif', style='normal', variant='normal', weight='normal', stretch='normal', size=20)
    font_properties = FontProperties()
    font_properties.set_size('small')

    config_items = get_configs_items()
    for key, value in config_items.items():
        rcParams[key] = value

    spectr_files_for_all_sources =["cepa_58441.31113425926_IRBENE16_241.h5", "afgl5142_58574.31668981481_IRBENE16_34.h5",
                                   "g107p3_58989.78732638889_IRBENE16_150.h5", "g111p26_58775.535532407404_IRBENE16_33.h5",
                                   "g196p454_59069.27920138889_IRBENE16_29.h5", "g22p357_58726.78696759259_IRBENE16_15.h5",
                                   "g24p33_58861.44650462963_IRBENE16_97.h5", "g25p65_58875.316655092596_IRBENE16_37.h5",
                                   "g25p710_58945.19710648148_IRBENE16_14.h5", "g30p99_58980.01193287037_IRBENE16_61.h5",
                                   "g32p045_58991.96428240741_IRBENE16_23.h5", "g32p745_58801.53355324074_IRBENE16_115.h5",
                                   "g33p64_59024.93394675926_IRBENE16_232.h5", "g34p403_58502.41447916667_IRBENE16_34.h5",
                                   "g35p20_58813.50585648148_IRBENE16_32.h5", "g36p7_58876.327627314815_IRBENE16_41.h5",
                                   "g37p43_59058.85618055556_IRBENE16_30.h5", "g37p479_59077.82991898148_IRBENE16_71.h5",
                                   "g37p55_59070.817199074074_IRBENE16_74.h5", "g43p79_58804.52097222222_IRBENE16_32.h5",
                                   "g45p071_59097.74915509259_IRBENE16_35.h5", "g49p04_58808.47662037037_IRBENE16_28.h5",
                                   "g59p783_59057.817141203705_IRBENE16_31.h5", "g73p06_59013.925625_IRBENE16_29.h5",
                                    "g75p78_58500.43746527778_IRBENE16_39.h5", "g78p12_59085.67820601852_IRBENE16_126.h5",
                                   "g85p41_59124.64989583333_IRBENE16_219.h5", "g90p92_58489.323159722226_IRBENE16_24.h5",
                                   "l1206_59125.47690972222_IRBENE16_82.h5", "l1287_59118.506840277776_IRBENE16_83.h5",
                                   "ngc281_59062.78135416667_IRBENE16_81.h5", "ngc7538_59125.46569444444_IRBENE16_90.h5",
                                   "on1_59119.714375_IRBENE16_82.h5", "s231_59125.44200231481_IRBENE16_59.h5",
                                   "s252_59125.419583333336_IRBENE16_60.h5", "s255_59125.40837962963_IRBENE16_53.h5",
                                   "v645_59120.74517361111_IRBENE16_92.h5", "w3oh_59119.739224537036_IRBENE16_86.h5",
                                   "w49n_59119.73453703704_IRBENE16_80.h5", "w51_59119.72194444444_IRBENE16_77.h5",
                                   "w75n_59060.821608796294_IRBENE16_73.h5"]

    density = []

    for file_name in spectr_files_for_all_sources:
        source = file_name.split("_")[0]
        file = get_configs("paths", "outputFilePath") + "/6668/" + source + "/"+ file_name
        spectr_data = h5py.File(file, 'r')['amplitude_corrected_not_smooht'][()]

        xdata_ = spectr_data[:, 0]
        ydata_ = spectr_data[:, 3]

        area = np.trapz(ydata_, xdata_)
        print(area)

        density.append(area)

    plt.hist(density)
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
