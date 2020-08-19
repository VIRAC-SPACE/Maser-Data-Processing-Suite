#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
If output file is changed manually, run this script to change result file
"""
import sys
import os
import argparse
import json
import h5py
import numpy as np

PACKAGE_PARENT = '..'

SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser
from utils.help import get_iteration_from_output_file


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Fix amplitudes in result file. ''')
    parser.add_argument("source", help="source name", type=str, default="")
    parser.add_argument("line", help="line", type=int)
    parser.add_argument("type", help="back-end type sdr or dbbc", type=str, default="SDR")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')
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


def get_local_max(data_tmp, source_velocities, index_range_for_local_maxima):
    data = data_tmp['amplitude_corrected_not_smooht'][()]
    x = data[:, 0]
    y_u1 = data[:, 1]
    y_u9 = data[:, 2]
    y_avg = data[:, 3]
    indexies_for_source_velocities = [0] * len(source_velocities)
    for index in range(0, len( source_velocities )):
        indexies_for_source_velocities[index] = (np.abs( x - float(source_velocities[index]))).argmin()

    max_amplitude_list_u1 = list()
    max_amplitude_list_u9 = list()
    max_amplitude_list_uavg = list()

    for index in indexies_for_source_velocities:
        max_amplitude_list_tmp_u1 = list()
        max_amplitude_list_tmp_u9 = list()
        max_amplitude_list_tmp_uavg = list()

        for i in range(index - index_range_for_local_maxima, index + index_range_for_local_maxima):
            max_amplitude_list_tmp_u1.append(y_u1[i])
            max_amplitude_list_tmp_u9.append(y_u9[i])
            max_amplitude_list_tmp_uavg.append(y_avg[i])

        max_amplitude_list_u1.append(max_amplitude_list_tmp_u1)
        max_amplitude_list_u9.append(max_amplitude_list_tmp_u9)
        max_amplitude_list_uavg.append(max_amplitude_list_tmp_uavg)

    max_amplitudes_u1 = [np.max(value) for value in max_amplitude_list_u1]
    max_amplitudes_u9 = [np.max(value) for value in max_amplitude_list_u9]
    max_amplitudes_uavg = [np.max(value) for value in max_amplitude_list_uavg]

    for max in range( 0, len( max_amplitudes_u1 ) ):
        max_amplitudes_u1[max] = [source_velocities[max], max_amplitudes_u1[max]]
        max_amplitudes_u9[max] = [source_velocities[max], max_amplitudes_u9[max]]
        max_amplitudes_uavg[max] = [source_velocities[max], max_amplitudes_uavg[max]]

    return (max_amplitudes_u1, max_amplitudes_u9, max_amplitudes_uavg)


def change_result_amplitudes(output_files, result_file, output_file_path, source_velocities, index_range_for_local_maxima, source, type_of_observation, line):
    with open(result_file) as result:
        result_json = json.load(result)

    for output_file in output_files:
        iteration = get_iteration_from_output_file(output_file)
        data_tmp = h5py.File(output_file_path + line + "/" + source + "/" +output_file, "r")
        if "amplitude_corrected_not_smooht" in data_tmp:
            max_apmlitudes_u1, max_apmlitudes_u9, max_apmlitudes_uavg = get_local_max(data_tmp,
                                                                                source_velocities, index_range_for_local_maxima)
            for key in result_json.keys():
                if key.endswith( "_" + str( iteration ) ):
                    if result_json[key]["type"] == type_of_observation:
                        result_json[key]["polarizationU1"] = max_apmlitudes_u1
                        result_json[key]["polarizationU9"] = max_apmlitudes_u9
                        result_json[key]["polarizationAVG"] = max_apmlitudes_uavg

        else:
            print( "Output " + output_file + " file has no amplitude_corrected_not_smooht colomm" )

    result_file = open(result_file, "w")
    result_file.write(json.dumps(result_json, indent=2))
    result_file.close()


def main():
    output_file_path = get_configs('paths', "outputFilePath")
    result_file_path = get_configs('paths', "resultFilePath")
    source_velocities = get_configs('velocities', get_args("source") + "_" + get_args("line")).replace(" ", "").split(",")
    index_range_for_local_maxima = int( get_configs('parameters', "index_range_for_local_maxima"))
    result_file = result_file_path + get_args("source") + "_" + get_args("line") + ".json"
    output_files = os.listdir(output_file_path + get_args("line") + "/" + get_args("source"))
    change_result_amplitudes(output_files, result_file, output_file_path, source_velocities, index_range_for_local_maxima, get_args("source"), get_args("type"), get_args("line"))

    sys.exit()


if __name__ == "__main__":
    main()
