#! /usr/bin/python3
# -*- coding: utf-8 -*-
import sys
import argparse

from PyQt5.QtWidgets import (QApplication)

from parsers.configparser_ import ConfigParser

from main_view import MainView


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


class MainGUI(QApplication):
    def __init__(self, sys_argv):
        super(MainGUI, self).__init__(sys_argv)
        source_name = get_args("source")
        line = get_args("line")
        self.main_view = MainView(source_name, line, get_args("config"))
        self.main_view.show()


def main():
    source_name = get_args("source")
    line = get_args("line")
    q_app = QApplication(sys.argv)
    main_view = MainView(source_name, line, get_args("config"))
    main_view.show()
    sys.exit(q_app.exec_())


if __name__ == "__main__":
    main()
    sys.exit(0)
