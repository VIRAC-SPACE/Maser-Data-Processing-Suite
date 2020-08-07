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


def find_log_file(log_list, iteration):
    """

    :param log_list: list of log files
    :param iteration: iteration of observations
    :param line: frequency
    :return:
    """
    tmpl = ""
    for log in log_list:
        if iteration in log:
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


def create_iteration_list(path, source, line):
    """

    :param line: frequency
    :param source: source
    :param path: input file path
    :return: None
    """
    iterations = [get_iteration(iteration) for iteration in os.listdir(path) if
                  source in iteration and line in iteration]
    iterations.sort(key=int, reverse=False)
    return iterations


def create_log_file_list(path, source, line):
    """

    :param line: frequency
    :param path: log file path
    :param source: source
    :return: all log files for source
    """
    return [log for log in os.listdir(path) if log.startswith(source) and line in log]


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

    processed_iteration = list()

    for experiment in result:
        if experiment.split("_")[-1] in sdr_iterations and \
                experiment.split("_")[-1] not in \
                processed_iteration and result[experiment]["type"] == "SDR":
            processed_iteration.append(experiment.split("_")[-1])

        if experiment.split("_")[-1] in processed_iteration and \
                experiment.split("_")[-1] in processed_iteration and \
                result[experiment]["type"] == "SDR" and result[experiment]["flag"] == True:
            processed_iteration.remove(experiment.split("_")[-1])

    processed_iteration.sort(key=int, reverse=False)

    for iteration in sdr_iterations:
        if iteration not in processed_iteration:
            log_file = find_log_file(logfile_list, iteration)
            sdr_fs_parameter = source_name + " " + line + " " + iteration + " " + log_file
            LOGGER.info("Executing python3 " + "sdr_fs.py " + sdr_fs_parameter)
            os.system("python3 " + "sdr_fs.py " + sdr_fs_parameter)

    output_files = os.listdir(output_path + "/" + line + "/" + source_name)
    for output_file in output_files:
        if output_file.split(".")[0].split("_")[-1] not in processed_iteration:
            if output_file.startswith(source_name):
                LOGGER.info("Executing python3 " +
                            "total_spectrum_analyzer_qt5.py " + output_file + " " + line)
                os.system("python3 " +
                          "total_spectrum_analyzer_qt5.py " + output_file + " " + line)


if __name__ == "__main__":
    main()
    sys.exit(0)
