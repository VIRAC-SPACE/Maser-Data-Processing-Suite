#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
'Check config file and add amplitude for new velocity.
"""
import json
import sys
import os
import argparse

import h5py
import numpy as np
from astropy.time import Time

PACKAGE_PARENT = '..'

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


def get_args(key):
    """

    :param key: argument key
    :return: to script passed argument value
    """
    return str(parse_arguments().__dict__[key])


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Check config file and add amplitude for new velocity. ''')
    parser.add_argument("source", help="source name", type=str, default="")
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')
    args = parser.parse_args()
    return args


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
    source = get_args("source")
    line = str(get_args("line"))
    source_velocities = get_configs('velocities', source + "_" + str(line)).replace(" ", "").split(",")
    output_path = get_configs("paths", "outputFilePath") + str(line) + "/" + source + "/"
    output_files = os.listdir(output_path)
    index_range_for_local_maxima = int(get_configs('parameters', "index_range_for_local_maxima"))

    result_file_name = source + "_" + str(line) + ".json"
    result_file_path = get_configs("paths", "resultFilePath")

    if os.path.isfile(result_file_path + result_file_name):
        pass
    else:
        os.system("touch " + result_file_path + result_file_name)

        result_file = open(result_file_path + result_file_name, "w")
        result_file.write("{ \n" + "\n}")
        result_file.close()

    with open(result_file_path + result_file_name) as result_data:
        result = json.load(result_data)

    for output_file in output_files[0:2]:
        output_file_data = output_file.replace(".h5", "").split("_")
        mjd = output_file_data[1]

        experiment_name = None
        experiment_names = result.keys()
        for experiment_name_ in experiment_names:
            if mjd == result[experiment_name_]['modifiedJulianDays']:
                experiment_name = experiment_name_
                break

        if experiment_name is not None:
            output_data = h5py.File(output_path + output_file, "r")
            if "amplitude_corrected_not_smooht" in output_data:

                total_results = output_data["amplitude_corrected_not_smooht"]
                velocity = total_results[:, 0]
                z1_not_smooht_data = total_results[:, 1]
                z2_not_smooht_data = total_results[:, 2]
                avg_y_not_smooht_data= total_results[:, 2]

                indexies_for_source_velocities = [0] * len(source_velocities)
                for index in range(0, len(source_velocities)):
                    indexies_for_source_velocities[index] = (
                        np.abs(velocity- float(source_velocities[index]))).argmin()

                max_amplitude_list_u1 = list()
                max_amplitude_list_u9 = list()
                max_amplitude_list_uavg = list()
                for index in indexies_for_source_velocities:
                    max_amplitude_list_tmp_u1 = list()
                    max_amplitude_list_tmp_u9 = list()
                    max_amplitude_list_tmp_uavg = list()
                    for i in range(index - index_range_for_local_maxima,
                                   index + index_range_for_local_maxima):
                        max_amplitude_list_tmp_u1.append(z1_not_smooht_data[i])
                        max_amplitude_list_tmp_u9.append(z2_not_smooht_data[i])
                        max_amplitude_list_tmp_uavg.append(avg_y_not_smooht_data[i])
                    max_amplitude_list_u1.append(max_amplitude_list_tmp_u1)
                    max_amplitude_list_u9.append(max_amplitude_list_tmp_u9)
                    max_amplitude_list_uavg.append(max_amplitude_list_tmp_uavg)

                max_apmlitudes_u1 = [np.max(value) for value in max_amplitude_list_u1]
                max_apmlitudes_u9 = [np.max(value) for value in max_amplitude_list_u9]
                max_apmlitudes_uavg = [np.max(value) for value in max_amplitude_list_uavg]

                for maximum in range(0, len(max_apmlitudes_u1)):
                    max_apmlitudes_u1[maximum] = [source_velocities[maximum], max_apmlitudes_u1[maximum]]
                    max_apmlitudes_u9[maximum] = [source_velocities[maximum], max_apmlitudes_u9[maximum]]
                    max_apmlitudes_uavg[maximum] = \
                        [source_velocities[maximum], max_apmlitudes_uavg[maximum]]

                result[experiment_name]["polarizationU1"] = max_apmlitudes_u1
                result[experiment_name]["polarizationU9"] = max_apmlitudes_u9
                result[experiment_name]["polarizationAVG"] = max_apmlitudes_uavg

            else:
                print(output_file + " wrong output file")

    sys.exit()


if __name__ == "__main__":
    main()