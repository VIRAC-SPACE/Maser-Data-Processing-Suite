import sys
import numpy as np
import h5py
import argparse


def parse_arguments():
    """

    :return: dict with passed args to script
    """
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''')
    parser.add_argument("inputfile", help="input file", type=str)
    parser.add_argument("outputfile", help="output file", type=str)
    parser.add_argument("table", help="table in output file", type=str )
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


def main():
    inputfile = get_args("inputfile")
    outputfile = get_args("outputfile")
    table = get_args("table")

    input_data_file = h5py.File(inputfile, "r")
    if table not in input_data_file.keys():
        print("Wrong table selected possible tables are " + str(list(input_data_file.keys())))
        sys.exit(1)

    data = input_data_file[table][()]
    np.savetxt(outputfile, data)
    sys.exit()


if __name__ == "__main__":
    main()
