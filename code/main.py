#! /usr/bin/python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import json
import coloredlogs, logging

from parsers._configparser import ConfigParser

coloredlogs.install(level='PRODUCTION')
logger = logging.getLogger('Main')


def parseArguments():
    parser = argparse.ArgumentParser(description='''automatically call frequencyShiftingAnalyzer and totalSpectrumAnalyer. ''', epilog="""Main program.""")
    parser.add_argument("source", help="Source Name", type=str)
    parser.add_argument("line", help="frequency", type=int)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-m", "--manual", help="Set manual log data", action='store_true')
    parser.add_argument("-n", "--noGUI", help="Create smoothed and not smothed outputfiles", action='store_true')
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 3.0')
    args = parser.parse_args()
    return args


def getArgs(key):
    return str(parseArguments().__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)


def findLogFile(logList, iteration):
    tmpL = -1
    for l in range(0, len(logList)):
        if logList[l].split("/")[-1].split(".")[0].split("_")[-1] == iteration:
            tmpL = l
            break
    if tmpL == -1:
        logger.warning("Warning " + "log for iteration " + iteration + " do not exist log file " + logList[-1] + " will be used instead!")
    return tmpL


def findLogFileSDR(logList, iteration, line):
    tmpL = -1
    for l in range(0, len(logList)):
        iter = logList[l].split("/")[-1].split(".")[0].split("_")[-1]
        lin = logList[l].split("/")[-1].split(".")[0].split("_")[-2]

        if lin + "_" + iter == "f" + line + "_" +iteration:
            tmpL = l
            break
    if tmpL == -1:
        logger.warning("Warning " + "log for iteration " + iteration + " do not exist log file " + logList[-1] + " will be used instead!")
    return tmpL


def createIterationList(path):
    iterations = [iteration for iteration in os.listdir(path)]
    iterations.sort(key=int, reverse=False)
    return iterations

def createLogFileList(path, source):
    return [log for log in os.listdir(path) if log.startswith(source)]
    
def main():
    args = parseArguments()
    sourceName = getArgs("source")
    dataFilesPath = getConfigs('paths', "dataFilePath")
    resultPath = getConfigs('paths', "resultFilePath")
    logPath = getConfigs('paths', "logPath")

    DBBCpath = dataFilesPath + "DBBC/"
    SDRpath = dataFilesPath + "SDR/"

    if os.path.exists(DBBCpath + sourceName + "/"):
        DBBC_iterations = createIterationList(DBBCpath + sourceName + "/")
    else:
        DBBC_iterations = []

    if os.path.exists(SDRpath + sourceName + "/f" + getArgs("line") + "/"):
        SDR_iterations = createIterationList(SDRpath + sourceName + "/f" + getArgs("line") + "/")
    else:
        SDR_iterations = []

    DBBClogPath = logPath + "DBBC/"
    SDRlOGPath = logPath + "SDR/"

    DBBClogfile_list = createLogFileList(DBBClogPath, sourceName)
    SDRlogfile_list = createLogFileList(SDRlOGPath, sourceName)

    resultFileName = resultPath + sourceName + "_" + getArgs("line") + ".json"
    print(resultFileName)

    if os.path.isfile(resultFileName):
        pass
    else:
        os.system("touch " + resultFileName)
        resultFile = open(resultFileName, "w")
        resultFile.write("{ \n" + "\n}")
        resultFile.close()

    with open(resultFileName, "r") as result_data:
        result = json.load(result_data)

    DBBCprocessed_iteration = list()
    SDRprocessed_iteration = list()

    for experiment in result:
        if experiment.split("_")[-1] in DBBC_iterations and experiment.split("_")[-1] not in DBBCprocessed_iteration and result[experiment]["type"] == "DBBC":
            DBBCprocessed_iteration.append(experiment.split("_")[-1])

    for experiment in result:
        if experiment.split("_")[-1] in SDR_iterations and experiment.split("_")[-1] not in SDRprocessed_iteration and result[experiment]["type"] == "SDR":
            SDRprocessed_iteration.append(experiment.split("_")[-1])

    DBBCprocessed_iteration.sort(key=int, reverse=False)
    SDRprocessed_iteration.sort(key=int, reverse=False)

    #DBBC process
    try:
        if args.manual:
            for i in DBBC_iterations:
                if i not in DBBCprocessed_iteration:
                    frequencyShiftingParametr = sourceName + " " + i + " " + str(DBBClogfile_list[findLogFile(DBBClogfile_list, i)])
                    logger.info("Executing python3 " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr + " -m")
                    os.system("python3  " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr + " -m")
        else:
            for i in DBBC_iterations:
                if i not in DBBCprocessed_iteration:
                    frequencyShiftingParametr = sourceName  + " " + i + " " + str(DBBClogfile_list[findLogFile(DBBClogfile_list, i)])
                    logger.info("Executing python3 " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr)
                    os.system("python3 " + "code/frequencyShiftingAnalyzer_qt5.py " + frequencyShiftingParametr)

        data_files = list()
        for data in os.listdir(DBBCpath):
            if data.startswith(sourceName) and data.endswith(".dat"):
                data_files.append(data)
        for d in data_files:
            if d.split(".")[0].split("_")[-1] not in DBBCprocessed_iteration:
                if args.noGUI:
                    logger.info("Executing python3 " + "code/totalSpectrumAnalyer_qt5.py " + " " + d + " " + getArgs("line") + " -n " + " -t DBBC")
                    os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + " " + d + " " + getArgs("line") + " -n " + " -t DBBC")
                else:
                    logger.info("Executing python3 " + "code/totalSpectrumAnalyer_qt5.py " + d + " " + getArgs("line")  + " -t DBBC")
                    os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + d + " " + getArgs("line")  + " -t DBBC")

    except IOError as e:
        print("IO Error", e)
        sys.exit(1)

    except IndexError as e:
        print("Index Error", e)
        sys.exit(1)

    except ValueError as e:
        print("Cannot crate modified Julian Days", e)

    except TypeError as e:
        print("TypeError", e)
        sys.exit(1)

    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(1)

    #SDR process
    try:
        if args.manual:
            for i in SDR_iterations:
                if i not in SDRprocessed_iteration:
                    frequencyShiftingParametr = sourceName + " " + getArgs("line") + " " + i + " " + str(SDRlogfile_list[findLogFileSDR(SDRlogfile_list, i, getArgs("line"))])
                    logger.info("Executing python3 " + "code/SDR_fs.py " + frequencyShiftingParametr + " -m")
                    os.system("python3  " + "code/SDR_fs.py " + frequencyShiftingParametr + " -m")
        else:
            for i in SDR_iterations:
                if i not in SDRprocessed_iteration:
                    frequencyShiftingParametr = sourceName + " " + getArgs("line") + " " + i + " " + str(SDRlogfile_list[findLogFileSDR(SDRlogfile_list, i, getArgs("line"))])
                    logger.info("Executing python3 " + "code/SDR_fs.py " + frequencyShiftingParametr)
                    os.system("python3 " + "code/SDR_fs.py " + frequencyShiftingParametr)

        data_files = list()
        for data in os.listdir(SDRpath):
            if data.startswith(sourceName) and data.endswith(".dat"):
                data_files.append(data)

        for d in data_files:
            if d.split(".")[0].split("_")[-1] not in SDRprocessed_iteration:
                if args.noGUI:
                    logger.info("Executing python3 " + "code/totalSpectrumAnalyer_qt5.py " + d + " " + getArgs("line") + " " + " -n")
                    os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + d + " " + getArgs("line") + " " + " -n")
                else:
                    logger.info("Executing python3 " + "code/totalSpectrumAnalyer_qt5.py " + d  + " " + getArgs("line"))
                    os.system("python3 " + "code/totalSpectrumAnalyer_qt5.py " + d + " " + getArgs("line"))

    except IOError as e:
        print("IO Error", e)
        sys.exit(1)

    except IndexError as e:
        print("Index Error", e)
        sys.exit(1)

    except ValueError as e:
        print("Cannot crate modified Julian Days", e)

    except TypeError as e:
        print("TypeError", e)
        sys.exit(1)

    except:
        print("Unexpected error:", sys.exc_info()[0])
        sys.exit(1)


if __name__=="__main__":
    main()
    sys.exit(0)
