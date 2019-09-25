import sys
import os
import argparse
import numpy as np

from parsers._configparser import ConfigParser
from help import *


def parseArguments():
    parser = argparse.ArgumentParser(description='''fix rt 32. ''', epilog="""fix rt 32.""")
    parser.add_argument("source", help="source name", type=str)
    parser.add_argument("line", help="source name", type=str)
    parser.add_argument("factor", help="factor", type=float)
    parser.add_argument("iteration_list", help="numbers of iterations to fix", type=str, nargs='+')
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    args = parser.parse_args()
    return args


def getArgs(key):
    args = parseArguments()
    return str(args.__dict__[key])


def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)


def correct_amplitude(file):
    factor = float(getArgs("factor"))
    data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file), 4))

    x = correctNumpyReadData(data[:, [0]])
    y_u1 = correctNumpyReadData(data[:, [1]]) * factor
    y_u9 = correctNumpyReadData(data[:, [2]]) * factor
    y_avg = correctNumpyReadData(data[:, [3]]) * factor

    totalResults = [x, y_u1, y_u9, y_avg]
    np.savetxt(file, np.transpose(totalResults))

def main():
    iteration_to_fix = getArgs("iteration_list").replace("[", "").replace("]", "").replace("'", "").split(",")
    iteration_to_fix = [i.strip() for i in iteration_to_fix]
    print(iteration_to_fix)
    output_dir = getConfigs("paths", "outputFilePath") + getArgs("source") + "/" + getArgs("line") + "/"

    for file in os.listdir(output_dir):
        if file.startswith(getArgs("source")) and file.split("_")[-1].split(".")[0] in iteration_to_fix and "IRBENE16" not in file:
            correct_amplitude(output_dir + file)


if __name__=="__main__":
    main()
    sys.exit()