#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QDesktopWidget, QPushButton, QInputDialog)
from PyQt5.QtGui import QIcon
from PyQt5 import QtCore
import argparse
import re
import numpy as np

from parsers._configparser import ConfigParser
from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes
from ploting_qt5 import Plot
from help import *

def parseArguments():
    parser = argparse.ArgumentParser(description='''Creates input file for plotting tool. ''', epilog="""PRE PLOTTER.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("iteration_number", help="iteration number ", type=int)
    parser.add_argument("logFile", help="Experiment log file name", type=str)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v", "--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    return args

def getArgs(key):
    return str(parseArguments().__dict__[key])

def getConfigs(key, value):
    configFilePath = getArgs("config")
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

class Analyzer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.setWindowTitle("SDR")
        self.center()
        self.DataDir = getConfigs("paths", "dataFilePath")  + "SDR/" + getArgs("source") + "/" + getArgs("iteration_number") + "/"
        self.DataFiles = os.listdir(self.DataDir)
        self.ScanPairs = self.createScanPairs()
        self.index = 0
        self.logs = LogReaderFactory.getLogReader(LogTypes.SDR, getConfigs("paths", "logPath") + getArgs("logFile"), getConfigs("paths", "prettyLogsPath") + getArgs("source") + "_" + getArgs("iteration_number")).getLogs()
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)

        self.__UI__()

    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def __getScanName__(self, dataFileName):
        return dataFileName.split(".")[0].split("_")[3][2:len(dataFileName.split(".")[0].split("_")[3])].lstrip('0')

    def __getData(self, dataFileName):
        data = np.fromfile(dataFileName, dtype="float64", count=-1, sep=" ").reshape((file_len(dataFileName), 3))
        frequency = correctNumpyReadData(data[:, [0]])
        polarizationU1 = correctNumpyReadData(data[:, [1]])
        polarizationU9 = correctNumpyReadData(data[:, [2]])
        return(frequency, polarizationU1, polarizationU9)

    def __getDataFileForScan__(self, scanName):
        fileName = ""

        for file in self.DataFiles:
            if  self.__getScanName__(file) == scanName:
                fileName = file
                break

        return fileName

    def createScanPairs(self):
        scanNames = [self.__getScanName__(file) for file in self.DataFiles]
        scansNumbers = list(set([int(re.findall("[0-9]+", s)[0]) for s in scanNames]))
        scansNumbers = sorted(scansNumbers)
        scanPairs = []

        for scan in scansNumbers:
            scanPairs.append( ((str(scan) + "r" + "0", str(scan) + "s" + "0"), (str(scan) + "r" + "1", str(scan) + "s" + "1")) )

        return scanPairs

    def plotPair(self, index):
        pair = self.ScanPairs[index]
        file1 = self.DataDir + self.__getDataFileForScan__(pair[0][0]) #r0
        file2 = self.DataDir + self.__getDataFileForScan__(pair[0][1]) #s0
        file3 = self.DataDir + self.__getDataFileForScan__(pair[1][0]) #r1
        file4 = self.DataDir + self.__getDataFileForScan__(pair[1][1]) #s1

        frequencyA = self.__getData(file1)[0] #r0
        polarizationU1A = self.__getData(file1)[1] #r0
        polarizationU9A = self.__getData(file1)[2] #r0

        frequencyB = self.__getData(file2)[0] #s0
        polarizationU1B = self.__getData(file2)[1] #s0
        polarizationU9B = self.__getData(file2)[2] #s0

        frequencyC = self.__getData(file3)[0] #r1
        polarizationU1C = self.__getData(file3)[1] #r1
        polarizationU9C = self.__getData(file3)[2] #r1

        frequencyD = self.__getData(file4)[0] #s1
        polarizationU1D = self.__getData(file4)[1] #s1
        polarizationU9D = self.__getData(file4)[2] #s1

        # fft shift
        polarizationU1A = np.fft.fftshift(polarizationU1A) #r0
        polarizationU9A = np.fft.fftshift(polarizationU9A) #r0
        polarizationU1B = np.fft.fftshift(polarizationU1B) #s0
        polarizationU9B = np.fft.fftshift(polarizationU9B) #s0
        polarizationU1C = np.fft.fftshift(polarizationU1C) #r1
        polarizationU9C = np.fft.fftshift(polarizationU9C) #r1
        polarizationU1D = np.fft.fftshift(polarizationU1D) #s1
        polarizationU9D = np.fft.fftshift(polarizationU9D) #s1

        #plot1
        self.plot_start_u1A = Plot()
        self.plot_start_u1A.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "u1 Polarization", (1, 0))
        self.plot_start_u1A.plot(frequencyA, polarizationU1A, 'b', label=str(index+1) + "r0")
        self.plot_start_u1A.plot(frequencyB, polarizationU1B, 'g', label=str(index+1) + "s0")
        self.plot_start_u1A.plot(frequencyC, polarizationU1C, 'r', label=str(index + 1) + "r1")
        self.plot_start_u1A.plot(frequencyD, polarizationU1D, 'y', label=str(index + 1) + "s1")
        self.grid.addWidget(self.plot_start_u1A, 0, 0)

        #plot2
        self.plot_start_u9B = Plot()
        self.plot_start_u9B.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "u9 Polarization", (1, 1))
        self.plot_start_u9B.plot(frequencyA, polarizationU9A, 'b', label=str(index+1) + "r0")
        self.plot_start_u9B.plot(frequencyB, polarizationU9B, 'g', label=str(index+1) + "s0")
        self.plot_start_u9B.plot(frequencyC, polarizationU9C, 'r', label=str(index + 1) + "r1")
        self.plot_start_u9B.plot(frequencyD, polarizationU9D, 'y', label=str(index + 1) + "s1")
        self.grid.addWidget(self.plot_start_u9B, 0, 1)

        df_div = 16
        BW = 1.5625
        f_shift = BW / df_div
        l_spec = len(frequencyA)
        f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
        n_shift = int(f_shift / f_step)
        avg_interval = 0.5 # inner 50%
        si = int(l_spec / 2 - l_spec * avg_interval / 2)
        ei = int(l_spec / 2 + l_spec * avg_interval / 2)

        Tsys_off_1U1 = float(self.logs["header"]["Tcal"][0]) * ((polarizationU1D + polarizationU1B) - np.mean(polarizationU1D[si:ei] - polarizationU1B[si:ei])) / (2 * np.mean(polarizationU1D[si:ei] - polarizationU1B[si:ei]))
        Tsys_off_2U1 = float(self.logs["header"]["Tcal"][1]) * ((polarizationU1C + polarizationU1A) - np.mean(polarizationU1C[si:ei] - polarizationU1A[si:ei])) / (2 * np.mean(polarizationU1C[si:ei] - polarizationU1A[si:ei]))

        Tsys_off_1U9 = float(self.logs["header"]["Tcal"][0]) * ((polarizationU9D + polarizationU9B) - np.mean(polarizationU9D[si:ei] - polarizationU9B[si:ei])) / (2 * np.mean(polarizationU9D[si:ei] - polarizationU9B[si:ei]))
        Tsys_off_2U9 = float(self.logs["header"]["Tcal"][1]) * ((polarizationU9C + polarizationU9A) - np.mean(polarizationU9C[si:ei] - polarizationU9A[si:ei])) / (2 * np.mean(polarizationU9C[si:ei] - polarizationU9A[si:ei]))

        Ta_1_caloffU1 = Tsys_off_1U1 * (polarizationU1A - polarizationU1B) / polarizationU1B # non-cal phase
        Ta_1_caloffU9 = Tsys_off_1U9 * (polarizationU9A - polarizationU9B) / polarizationU9B  # non-cal phase

        Ta_1_calonU1 = (Tsys_off_1U1 + float(self.logs["header"]["Tcal"][0])) * (polarizationU1C - polarizationU1D) / polarizationU1D  # cal phase
        Ta_1_calonU9 = (Tsys_off_1U9 + float(self.logs["header"]["Tcal"][1])) * (polarizationU9C - polarizationU9D) / polarizationU9D  # cal phase

        Ta_sigU1 = (Ta_1_caloffU1 + Ta_1_calonU1) / 2
        Ta_sigU9 = (Ta_1_caloffU9 + Ta_1_calonU9) / 2

        Ta_2_caloffU1 = Tsys_off_2U1 * (polarizationU1A - polarizationU1A) / polarizationU1A  # non-cal phase
        Ta_2_caloffU9 = Tsys_off_2U9 * (polarizationU9A - polarizationU1A) / polarizationU1A  # non-cal phase

        Ta_2_calonU1 = (Tsys_off_2U1 + float(self.logs["header"]["Tcal"][0])) * (polarizationU1D - polarizationU1C) / polarizationU1C  # cal phase
        Ta_2_calonU9 = (Tsys_off_2U9 + float(self.logs["header"]["Tcal"][1])) * (polarizationU9D - polarizationU9C) / polarizationU9C  # cal phase

        Ta_refU1 = (Ta_2_caloffU1 + Ta_2_calonU1) / 2
        Ta_refU9 = (Ta_2_caloffU9 + Ta_2_calonU9) / 2;

        Ta_sigU1 = np.roll(Ta_sigU1, +n_shift)
        Ta_sigU9 = np.roll(Ta_sigU9, +n_shift)

        Ta_refU1 = np.roll(Ta_refU1, -n_shift)
        Ta_refU9 = np.roll(Ta_refU9, -n_shift)

        TaU1 = (Ta_sigU1 + Ta_refU1) / 2
        TaU9 = (Ta_sigU9 + Ta_refU9) / 2

        El = 80
        G_El = [-0.0000333143, 0.0033676682, 0.9144626256]
        SfU1 = TaU1 / float(self.logs["header"]["DPFU"][0]) * np.polyval(G_El, El)
        SfU9 = TaU9 / float(self.logs["header"]["DPFU"][1]) * np.polyval(G_El, El)

        # plot3
        self.total_u1 = Plot()
        self.total_u1.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "", (4, 0))
        self.total_u1.plot(frequencyA[si:ei], SfU1[si:ei], 'b', label=str(index + 1))
        self.grid.addWidget(self.total_u1, 3, 0)

        # plot4
        self.total__u9 = Plot()
        self.total__u9.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "", (4, 1))
        self.total__u9.plot(frequencyA[si:ei], SfU9[si:ei], 'b', label=str(index + 1))
        self.grid.addWidget(self.total__u9, 3, 1)

        if index == len(self.ScanPairs) - 1:
            self.nextPairButton.setText('Move to total results')
            self.nextPairButton.clicked.connect(self.plotTotalResults)

    def plotTotalResults(self):
        pass

    def nextPair(self):
        if self.index == len(self.ScanPairs)- 1:
            pass

        else:
            self.plot_start_u1A.removePolt()
            self.plot_start_u9B.removePolt()
            self.index = self.index + 1

            self.plotPair(self.index)

    def __UI__(self):
        if self.index != len(self.ScanPairs) - 1:  # cheking if there is not one pair
            self.nextPairButton = QPushButton("Next pair", self)
            self.nextPairButton.clicked.connect(self.nextPair)
            self.grid.addWidget(self.nextPairButton, 4, 2)

        self.plotPair(self.index)

def main():
    qApp = QApplication(sys.argv)
    a = Analyzer()
    a.show()
    sys.exit(qApp.exec_())
    sys.exit(0)

if __name__ == "__main__":
    main()

