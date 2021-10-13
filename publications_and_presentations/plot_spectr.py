import sys
import os
import argparse
import h5py
import matplotlib.pyplot as plt

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.append(os.path.normpath(os.path.join( SCRIPT_DIR, PACKAGE_PARENT)))


from parsers.configparser_ import ConfigParser


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''')
    parser.add_argument("outputfile", help="output file", type=str)
    parser.add_argument("table", help="table in output file", type=str)
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="frequency", type=int)
    parser.add_argument("-c", "--config", help="Configuration cfg file",
                        type=str, default="config/config.cfg")
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


def main():
    plt.style.use('config/plot.style')
    table = get_args("table")
    if table not in ["amplitude_corrected", "amplitude_corrected_not_smooht"]:
        print("Wrong table specified. Table not in ['amplitude_corrected', 'amplitude_corrected_not_smooht']")
        sys.exit(1)

    line = get_args("line")
    source = get_args("source")
    output_path = get_configs("paths", "outputFilePath")
    output_file_path = output_path + line + "/" + source + "/" + get_args("outputfile")
    data = h5py.File(output_file_path, "r")[table][()]
    plt.plot(data[:,0], data[:, -1])
    plt.show()
    sys.exit(0)


if __name__ == "__main__":
    main()
