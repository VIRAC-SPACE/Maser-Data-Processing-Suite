#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
Delete flagged observations
"""
import sys
import os
import argparse
import shutil
import json

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join( SCRIPT_DIR, PACKAGE_PARENT)))

from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Delete flagged observations. ''')
    parser.add_argument("source", help="source name", type=str)
    parser.add_argument("line", help="Observed frequency", type=int)
    parser.add_argument("-c", "--config", help="Configuration cfg file",
                        type=str, default="config/config.cfg")
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


def get_iteration_from_output_file(file):
    return int(file.split(".")[1].split("_")[-1])


def get_station_from_output_file(file):
    return file.split(".")[1].split("_")[-2]


def create_info_dict_from_output_file(file):
    return {"iteration_number": get_iteration_from_output_file(file),
            "station": get_station_from_output_file(file)}


def get_flagged_result_name(result_data, iteration, station):
    exper = ""
    for experiment in result_data:
        if str(iteration) in experiment and station in experiment:
            exper = experiment
            break
    return exper


class Experiment:
    """
     Experiment class
    """

    def __init__(self, **entries):
        self.flag = None
        self.modifiedJulianDays = None
        self.Iteration_number = None
        self.polarizationAVG = None
        self.polarizationU9 = None
        self.polarizationU1 = None
        self.location = None
        self.Date = None
        self.specie = None
        self.type = None
        self.__dict__.update(entries)


def main():
    """
    :return: None
    """

    source = get_args("source")
    line = get_args("line")
    result_file_path = get_configs("paths", "resultFilePath")
    result_file_name = result_file_path + source + "_" + line + ".json"
    with open(result_file_name, "r") as result_file:
        result_data = json.load(result_file)

    experiments = [Experiment(**result_data[experiment]) for experiment in result_data]
    flagged_experiment_info = [{"iteration_number": experiment.Iteration_number, "station": experiment.location} for
                               experiment in experiments if experiment.flag]
    output_file_dir = get_configs("paths", "outputFilePath") + line + "/" + source + "/"
    flag_output_files = [file for file in os.listdir(output_file_dir)
                         if create_info_dict_from_output_file(file) in flagged_experiment_info]

    for file in flag_output_files:
        del_outpu_file = output_file_dir + file
        choice = input("Should file " + del_outpu_file + " be deleted Y/n ")
        if choice == "Y" or choice == "y":
            os.remove(del_outpu_file)
            print(del_outpu_file + " are deleted")

    for flag_info in flagged_experiment_info:
        iteration = flag_info["iteration_number"]
        station =  flag_info["station"]
        exper = get_flagged_result_name(result_data, iteration, station)
        choice2 = input("Should this experiment  " + exper + " be deleted  from result file Y/n " )
        if choice2 == "Y" or choice2 == "y":
            del result_data[exper]
            print("experiment  " + exper + " are deleted from result file " + result_file_name)

        if station == "IRBENE":
            st = "ir"
        elif station == "IRBENE16":
            st = "ib"
        data_file_dir = get_configs("paths", "dataFilePath") + source + "_f" + line + "_" + st + "_" + str(iteration)
        choice3 = input("Should this data file directory "  + data_file_dir + " be deleted Y/n ")
        if choice3 == "Y" or choice3 == "y":
            shutil.rmtree(data_file_dir)
            print("Data file directory " + data_file_dir + " are deleted")

    sys.exit(0)


if __name__ == "__main__":
    main()
