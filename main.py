#! /usr/bin/python3
# -*- coding: utf-8 -*-

"""
run scripts sdr_fs.py and total_spectrum_analyzer_qt5.py for all unprocessed experiments
"""
import os
import sys
import argparse
import json
import logging
import coloredlogs
import h5py
from parsers.configparser_ import ConfigParser

coloredlogs.install(level='PRODUCTION')
LOGGER = logging.getLogger('Main')


def parse_arguments():
    """
    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''automatically call sdr_fs.py 
    and total_spectrum_analyzer_qt5.py.''', epilog="""Main program.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="frequency", type=int)
    parser.add_argument("-c", "--config", help="Configuration "
                                                "cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 3.0')
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


def find_log_file(log_list, iteration, station):
    """

    :param log_list: list of log files
    :param iteration: iteration of observations
    :param station: station of observation
    :return:
    """
    tmpl = ""
    for log in log_list:
        if "_" + str(iteration) + ".log" in log and station in log:
            tmpl = log
            break

    else:
        tmpl = log_list[-1]

    if tmpl == "":
        LOGGER.warning("Warning " + "log for iteration " +
                       iteration + " do not exist log file " +
                       log_list[-1] + " will be used instead!")
    return tmpl


def get_iteration(dir_name):
    """

    :param dir_name:
    :return: iteration number
    """
    return dir_name.split("_")[-1]


def get_station(dir_name):
    """

    :param dir_name:
    :return: station
    """
    return dir_name.split("_")[-2]


def create_iteration_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: iterations list
    """
    stations = list(set(create_station_list(path, source, line)))
    iterations_for_source_and_line = [file for file in os.listdir(path)
                                      if source in file and line in file and os.path.isdir(path + file)]
    iterations_for_station = {station: [] for station in stations} #{get_station(iteration):get_iteration(iteration) for iteration in os.listdir(path) if  source in iteration and line in iteration and os.path.isdir(path + iteration)}
    for iteration in iterations_for_source_and_line:
        iterations_for_station[get_station(iteration)].append(get_iteration(iteration))

    for station in stations:
        iterations_for_station[station].sort(key=int, reverse=False)

    return iterations_for_station


def create_station_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: stations list
    """
    iterations = [file for file in os.listdir(path) if source in file and line in file and os.path.isdir(path + file)]
    iterations.sort(key=get_iteration, reverse=False)
    stations = [get_station(iteration) for iteration in iterations]
    return stations


def create_log_file_list(path, source, line):
    """

    :param line: frequency
    :param path: log file path
    :param source: source
    :return: all log files for source
    """
    return [log for log in os.listdir(path) if log.startswith(source + "_") and line in log]


def main():
    """
    :return: None
    """
    source_name = get_args("source")
    line = get_args("line")
    data_files_path = get_configs('paths', "dataFilePath")
    result_path = get_configs('paths', "resultFilePath")
    log_path = get_configs('paths', "logPath")
    output_path = get_configs('paths', "outputFilePath")

    if os.path.exists(data_files_path):
        sdr_iterations = create_iteration_list(data_files_path, source_name, line)
    else:
        sdr_iterations = []

    log_path = log_path + "SDR/"
    logfile_list = create_log_file_list(log_path, source_name, line)
    result_file_name = result_path + source_name + "_" + get_args("line") + ".json"

    if os.path.isfile(result_file_name):
        pass
    else:
        os.system("touch " + result_file_name)
        result_file = open(result_file_name, "w")
        result_file.write("{ \n" + "\n}")
        result_file.close()

    with open(result_file_name, "r") as result_data:
        result = json.load(result_data)

    stations = list(set(create_station_list(data_files_path, source_name, line)))
    processed_iteration = {station: [] for station in stations}
    processed_iteration2 = {station: [] for station in stations}
    for experiment in result:
        if get_station(experiment) == "IRBENE16":
            station = "ib"
        else:
            station = "ir"

        iteration_in_result = experiment.split("_")[-1]
        if iteration_in_result in sdr_iterations[station] and \
                iteration_in_result not in \
                processed_iteration[station] and result[experiment]["type"] == "SDR":
            processed_iteration[station].append(get_iteration(experiment))

        if iteration_in_result in processed_iteration[station] and \
                result[experiment]["type"] == "SDR" and result[experiment]["flag"]:
            processed_iteration[station].remove(iteration_in_result)

        if iteration_in_result not in processed_iteration2[station] and result[experiment]["type"] == "SDR":
            processed_iteration2[station].append(iteration_in_result)

    for station in stations:
        processed_iteration[station].sort(key=int, reverse=False)
        processed_iteration2[station].sort(key=int, reverse=False)

    for station in stations:
        for iteration in sdr_iterations[station]:
            if iteration not in processed_iteration[station]:
                log_file = find_log_file(logfile_list, iteration, station)
                sdr_fs_parameter = source_name + " " + line + " " + iteration + " " + log_file
                LOGGER.info("Executing python3 " + "sdr_fs.py " + sdr_fs_parameter)
                os.system("python3 " + "sdr_fs.py " + sdr_fs_parameter)

    output_files = os.listdir(output_path + "/" + line + "/" + source_name)
    for output_file in output_files:
        output_file_station = output_file.split("_")[-2].split(".")[0]
        output_file_iteration = output_file.split("_")[-1].split(".")[0]
        if output_file_station == "IRBENE16":
            station = "ib"
        else:
            station = "ir"

        if output_file_iteration not in processed_iteration2[station]:
            if output_file.startswith(source_name):
                with h5py.File(get_configs("paths", "outputFilePath") + get_args("line") +
                               "/" + get_args("source") + "/" + output_file, "r") as input_data_file:
                    input_file_keys = list(input_data_file.keys())
                    input_data_file.close()
                    if "amplitude" in input_file_keys:
                        LOGGER.info("Executing python3 " +
                                   "total_spectrum_analyzer_qt5.py " + output_file + " " + line)
                        os.system("python3 " +
                                  "total_spectrum_analyzer_qt5.py " + output_file + " " + line)


if __name__ == "__main__":
    main()
    sys.exit(0)
