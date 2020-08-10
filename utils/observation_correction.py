#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Correct observation by given factor
"""
import sys
import os
import argparse
import json
import h5py
from help import Experiment, get_iteration_from_output_file

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Correct observation by given factor. ''')
    parser.add_argument("source", help="source name", type=str)
    parser.add_argument("line", help="observed frequency", type=str)
    parser.add_argument("factor", help="factor", type=float)
    parser.add_argument("station", help="station", type=str)
    parser.add_argument("type", help="back-end type sdr or dbbc", type=str)
    parser.add_argument("iteration_list", help="numbers of iterations to fix", type=str, nargs='+')
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
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


def get_mjd_from_output_file(file):
    return file.split("_")[1]


def correct_output_file(output_dir, file, factor):
    output_data = h5py.File(output_dir + file, "a")
    amplitude_corrected = output_data["amplitude_corrected"]
    amplitude_corrected_not_smooht = output_data["amplitude_corrected_not_smooht"]
    amplitude_corrected[...] = output_data["amplitude_corrected"][()] * factor
    amplitude_corrected_not_smooht[...] = output_data["amplitude_corrected_not_smooht"][()] * factor
    output_data.close()


def correct_result_file(result_file_name, iteration_to_fix, factor, station, type):
    with open(result_file_name, "r") as result_file:
        result_data = json.load(result_file)

    for experiment in result_data:
        if experiment.split("_")[-1] in iteration_to_fix and \
                station in experiment and \
                result_data[experiment]["type"] == type:
            for v in range(0, len(result_data[experiment]["polarizationU1"])):
                result_data[experiment]["polarizationU1"][v][1] = result_data[experiment]["polarizationU1"][v][1] * factor
                result_data[experiment]["polarizationU9"][v][1] = result_data[experiment]["polarizationU9"][v][1] * factor
                result_data[experiment]["polarizationAVG"][v][1] = result_data[experiment]["polarizationAVG"][v][1] * factor

    result_file = open(result_file_name, "w")
    result_file.write(json.dumps(result_data, indent=2))
    result_file.close()


def main():
    iteration_to_fix = get_args("iteration_list"). \
        replace("[", "").replace("]", ""). \
        replace("'", "").split(",")

    type = get_args("type").upper()
    iteration_to_fix = [int(i.strip()) for i in iteration_to_fix]
    result_file_path = get_configs("paths", "resultFilePath")
    source = get_args("source")
    line = get_args("line")
    station = get_args("station").upper()
    factor = float(get_args("factor"))
    result_file_name = result_file_path + source + "_" + line + ".json"
    with open(result_file_name, "r") as result_file:
        result_data = json.load(result_file)

    output_dir = get_configs("paths", "outputFilePath") + "/" + get_args("line") + "/" + get_args("source") + "/"
    experiments = [Experiment(**result_data[experiment]) for experiment in result_data]

    mdj_for_experiments_to_fix = [experiment.modifiedJulianDays for experiment in experiments
                                  if experiment.Iteration_number in iteration_to_fix and
                                  experiment.location == station and
                                  experiment.type == type
                                  ]

    for file in os.listdir(output_dir):
        if get_iteration_from_output_file(file) in iteration_to_fix and \
                get_mjd_from_output_file(file) in mdj_for_experiments_to_fix:
            correct_output_file(output_dir, file,factor)

    correct_result_file( result_file_name, iteration_to_fix, factor, station, type )
    sys.exit()


if __name__ == "__main__":
    main()
