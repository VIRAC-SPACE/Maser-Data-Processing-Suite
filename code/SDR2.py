#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QDesktopWidget, QPushButton, QInputDialog)
from PyQt5.QtGui import QIcon
import argparse
import re
import numpy as np
import scipy.constants
from astropy.time import Time
import datetime
import pickle

from parsers._configparser import ConfigParser
from ExperimentsLogReader.experimentsLogReader import LogReaderFactory, LogTypes
from ploting_qt5 import Plot
from vlsr import lsr
from help import *

def parseArguments():
    parser = argparse.ArgumentParser(description='''Creates input file for plotting tool. ''', epilog="""PRE PLOTTER.""")
    parser.add_argument("source", help="Experiment source", type=str, default="")
    parser.add_argument("line", help="frequency", type=str)
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

def dopler(ObservedFrequency, velocityReceiver, f0):
    c = scipy.constants.speed_of_light
    velocitySoure = (-((ObservedFrequency / f0) - 1) * c + (velocityReceiver * 1000)) / 1000
    return velocitySoure

class Result(object):
    __slots__ = ('matrix', 'specie')
    def __init__(self, matrix, specie):
        self.matrix = matrix
        self.specie = specie

    def getMatrix(self):
        return self.matrix

    def getSpecie(self):
        return self.specie

class Analyzer(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.setWindowTitle("SDR")
        self.center()
        self.DataDir = getConfigs("paths", "dataFilePath") + "SDR/" + getArgs("source") + "/f" + getArgs("line") + "/" + getArgs("iteration_number") + "/"
        self.DataFiles = os.listdir(self.DataDir)
        self.ScanPairs = self.createScanPairs()
        self.index = 0
        self.Sf_left = list()
        self.Sf_right = list()
        self.logs = LogReaderFactory.getLogReader(LogTypes.SDR, getConfigs("paths", "logPath") + "SDR/"+ getArgs("logFile"), getConfigs("paths", "prettyLogsPath") + getArgs("source") + "_" + getArgs("iteration_number")).getLogs()
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
        polarization_left = correctNumpyReadData(data[:, [1]])
        polarization_right = correctNumpyReadData(data[:, [2]])
        return(frequency, polarization_left, polarization_right)

    def __getDataFileForScan__(self, scanName):
        fileName = ""

        for file in self.DataFiles:
            if self.__getScanName__(file) == scanName:
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

        print("data files", file1, file2, file3, file4)

        frequencyA = self.__getData(file1)[0] #r0
        p_sig_left = self.__getData(file1)[1] #r0
        p_sig_right = self.__getData(file1)[2] #r0

        frequencyB = self.__getData(file2)[0] #s0
        p_ref_left = self.__getData(file2)[1] #s0
        p_ref_right = self.__getData(file2)[2] #s0

        frequencyC = self.__getData(file3)[0] #r1
        p_sig_on_left = self.__getData(file3)[1] #r1
        p_sig_on_right = self.__getData(file3)[2] #r1

        frequencyD = self.__getData(file4)[0] #s1
        p_ref_on_left = self.__getData(file4)[1] #s1
        p_ref_on_right = self.__getData(file4)[2] #s1

        # fft shift
        p_sig_left = np.fft.fftshift(p_sig_left) #r0
        p_sig_right  = np.fft.fftshift(p_sig_right) #r0
        p_ref_left = np.fft.fftshift(p_ref_left) #s0
        p_ref_right = np.fft.fftshift(p_ref_right) #s0
        p_sig_on_left = np.fft.fftshift(p_sig_on_left) #r1
        p_sig_on_right = np.fft.fftshift(p_sig_on_right) #r1
        p_ref_on_left = np.fft.fftshift(p_ref_on_left) #s1
        p_ref_on_right = np.fft.fftshift(p_ref_on_right) #s1

        #plot1
        self.plot_start__leftA = Plot()
        self.plot_start__leftA.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "Left Polarization", (1, 0), "linear")
        self.plot_start__leftA.plot(frequencyA, p_sig_left, 'b', label=str(index+1) + "r0")
        self.plot_start__leftA.plot(frequencyB, p_ref_left, 'g', label=str(index+1) + "s0")
        self.plot_start__leftA.plot(frequencyC, p_sig_on_left, 'r', label=str(index+1) + "r1")
        self.plot_start__leftA.plot(frequencyD, p_ref_on_left, 'y', label=str(index+1) + "s1")
        self.grid.addWidget(self.plot_start__leftA, 0, 0)

        #plot2
        self.plot_start__rightB = Plot()
        self.plot_start__rightB.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "Right Polarization", (1, 1), "linear")
        self.plot_start__rightB.plot(frequencyA, p_sig_right , 'b', label=str(index+1) + "r0")
        self.plot_start__rightB.plot(frequencyB, p_ref_right, 'g', label=str(index+1) + "s0")
        self.plot_start__rightB.plot(frequencyC, p_sig_on_right, 'r', label=str(index+1) + "r1")
        self.plot_start__rightB.plot(frequencyD, p_ref_on_right, 'y', label=str(index+1) + "s1")
        self.grid.addWidget(self.plot_start__rightB, 0, 1)

        df_div = float(self.logs["header"]["df_div,df"][0])
        BW = float(self.logs["header"]["Fs,Ns,RBW"][0])
        f_shift = BW / df_div
        l_spec = len(frequencyA)
        f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
        n_shift = int(np.rint(f_shift / f_step))
        avg_interval = 0.5 # inner 50%
        si = int(l_spec / 2 - l_spec * avg_interval / 2)
        ei = int(l_spec / 2 + l_spec * avg_interval / 2)

        Tsys_off_1_left = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_left + p_ref_left) - np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei])) / (2 * np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei]))
        Tsys_off_2_left = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_left + p_sig_left) - np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei])) / (2 * np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei]))

        Tsys_off_1_right = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_right + p_ref_right) - np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei])) / (2 * np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei]))
        Tsys_off_2_right = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_right + p_sig_right ) - np.mean(p_sig_on_right[si:ei] - p_sig_right [si:ei])) / (2 * np.mean(p_sig_on_right[si:ei] - p_sig_right [si:ei]))

        Ta_1_caloff_left = Tsys_off_1_left * (p_sig_left - p_ref_left) / p_ref_left # non-cal phase
        Ta_1_caloff_right = Tsys_off_1_right * (p_sig_right - p_ref_right) / p_ref_right  # non-cal phase

        Ta_1_calon_left = (Tsys_off_1_left + float(self.logs["header"]["Tcal"][0])) * (p_sig_on_left - p_ref_on_left) / p_ref_on_left  # cal phase
        Ta_1_calon_right = (Tsys_off_1_right + float(self.logs["header"]["Tcal"][1])) * (p_sig_on_right - p_ref_on_right) / p_ref_on_right  # cal phase

        Ta_sig_left = (Ta_1_caloff_left + Ta_1_calon_left) / 2
        Ta_sig_right = (Ta_1_caloff_right + Ta_1_calon_right) / 2

        Ta_2_caloff_left = Tsys_off_2_left * (p_ref_left - p_sig_left) / p_sig_left  # non-cal phase
        Ta_2_caloff_right = Tsys_off_2_right * (p_ref_right - p_sig_right ) / p_sig_right   # non-cal phase

        Ta_2_calon_left = (Tsys_off_2_left + float(self.logs["header"]["Tcal"][0])) * (p_ref_on_left - p_sig_on_left) / p_sig_on_left  # cal phase
        Ta_2_calon_right = (Tsys_off_2_right + float(self.logs["header"]["Tcal"][1])) * (p_ref_on_right - p_sig_on_right) / p_sig_on_right  # cal phase

        Ta_ref_left = (Ta_2_caloff_left + Ta_2_calon_left) / 2
        Ta_ref_right = (Ta_2_caloff_right + Ta_2_calon_right) / 2

        Ta_sig_left = np.roll(Ta_sig_left, +n_shift)
        Ta_sig_right = np.roll(Ta_sig_right, +n_shift)

        Ta_ref_left = np.roll(Ta_ref_left, -n_shift)
        Ta_ref_right = np.roll(Ta_ref_right, -n_shift)

        Ta_left = (Ta_sig_left + Ta_ref_left) / 2
        Ta_right = (Ta_sig_right + Ta_ref_right) / 2

        El = (float(self.logs[pair[0][0]]["AzEl"][1]) + float(self.logs[pair[0][1]]["AzEl"][1]) + float(self.logs[pair[1][0]]["AzEl"][1]) + float(self.logs[pair[1][1]]["AzEl"][1]))/ 4
        G_El = self.logs["header"]["Elev_poly"]
        G_El = [float(gel) for gel in G_El]
        G_ELtmp = [0,0,0]
        G_ELtmp[0] = G_El[2]
        G_ELtmp[1] = G_El[1]
        G_ELtmp[2] = G_El[0]
        G_El = G_ELtmp

        Sf_leftscan = Ta_left / (((-1) * float(self.logs["header"]["DPFU"][0])) * np.polyval(G_El, El))
        Sf_rightscan = Ta_right / (((-1) * float(self.logs["header"]["DPFU"][1])) * np.polyval(G_El, El))

        self.si = si
        self.ei = ei

        self.Sf_left.append(Sf_leftscan[self.si:self.ei])
        self.Sf_right.append(Sf_rightscan[self.si:self.ei])

        # plot3
        self.total__left = Plot()
        self.total__left.creatPlot(self.grid, 'Frequency Mhz', 'Flux density (Jy)', "", (4, 0), "linear")
        self.x = frequencyA[self.si:self.ei]
        self.total__left.plot(frequencyA[self.si:self.ei], Sf_leftscan[self.si:self.ei], 'b', label=str(index + 1))
        self.grid.addWidget(self.total__left, 3, 0)

        # plot4
        self.total__right = Plot()
        self.total__right.creatPlot(self.grid, 'Frequency Mhz', 'Flux density (Jy)', "", (4, 1), "linear")
        self.total__right.plot(frequencyA[self.si:self.ei], Sf_rightscan[self.si:self.ei], 'b', label=str(index + 1))
        self.grid.addWidget(self.total__right, 3, 1)

        if index == len(self.ScanPairs) - 1:
            self.nextPairButton.setText('Move to total results')
            self.nextPairButton.clicked.connect(self.plotTotalResults)

    def plotTotalResults(self):
        self.grid.removeWidget(self.plot_start__leftA)
        self.grid.removeWidget(self.plot_start__rightB)
        self.grid.removeWidget(self.total__left)
        self.grid.removeWidget(self.total__right)

        self.plot_start__leftA.hide()
        self.plot_start__rightB.hide()
        self.total__left.hide()
        self.total__right.hide()

        self.plot_start__leftA.close()
        self.plot_start__rightB.close()
        self.total__left.close()
        self.total__right.close()

        self.plot_start__leftA.removePolt()
        self.plot_start__rightB.removePolt()
        self.total__left.removePolt()
        self.total__right.removePolt()

        del self.plot_start__leftA
        del self.plot_start__rightB
        del self.total__left
        del self.total__right

        self.grid.removeWidget(self.nextPairButton)
        self.nextPairButton.hide()
        self.nextPairButton.close()
        del self.nextPairButton

        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()

        velocitys_avg = np.zeros(len(self.x))
        y__left_avg = np.zeros(len(self.x))
        y__right_avg = np.zeros(len(self.x))

        station = self.logs["header"]["station"]
        if station == "RT-32":
            station = "IRBENE"
        else:
            station = "IRBENE16"

        stationCordinations = getConfigs("stations", station)
        stationCordinations = stationCordinations.replace(" ", "").split(",")
        x = np.float64(stationCordinations[0])
        y = np.float64(stationCordinations[1])
        z = np.float64(stationCordinations[2])

        maximums_left = []
        minimums_left = []
        maximums_right = []
        minimums_right = []

        for p in range(0, len(self.ScanPairs)):
            scanNumber = self.ScanPairs[p][0][0]
            scan_1 = self.logs[str(scanNumber)]
            stringTime = scan_1["date"].replace("T", " ")
            t = datetime.datetime.strptime(scan_1["date"], '%Y-%m-%dT%H:%M:%S')
            time = t.isoformat()
            date = Time(time, format='isot', scale='utc')

            print("stationCordinations", stationCordinations)

            sourceCordinations = getConfigs("sources",  getArgs("source")).split(",")
            sourceCordinations = [sc.strip() for sc in sourceCordinations]
            RA = sourceCordinations[0]
            DEC = sourceCordinations[1]

            ra = list()
            dec = list()
            ra.append(RA[0:2])
            ra.append(RA[2:4])
            ra.append(RA[4:len(RA)])

            if DEC[0] == "-":
                dec.append(DEC[0:3])
                dec.append(DEC[3:5])
                dec.append(DEC[5:len(DEC)])
            else:
                dec.append(DEC[0:2])
                dec.append(DEC[2:4])
                dec.append(DEC[4:len(DEC)])

            RaStr = ra[0] + "h" + ra[1] + "m" + ra[2] + "s"
            if int(dec[0]) > 0:
                DecStr = "+" + dec[0] + "d" + dec[1] + "m" + dec[2] + "s"
            else:
                DecStr = dec[0] + "d" + dec[1] + "m" + dec[2] + "s"

            print("Vel Total params", RaStr, DecStr, date, stringTime, x, y, z)
            VelTotal = lsr(RaStr, DecStr, date, stringTime, x, y, z)
            print("VelTotal", VelTotal)

            self.max_y_left_index = self.Sf_left[p].argmax(axis=0)
            self.max_y_right_index = self.Sf_right[p].argmax(axis=0)

            line = getConfigs('base_frequencies_SDR', "f" + getArgs("line")).replace(" ", "").split(",")
            lineF = float(line[0]) * (10**9)
            lineS = line[1]
            specie = lineS

            print("specie", specie, "\n")
            LO = float(self.logs["header"]["f_obs,LO,IF"][1])
            velocitys = dopler((self.x + LO) * (10 ** 6), VelTotal, lineF)

            maximums_left.append(np.max(self.Sf_left[p]))
            minimums_left.append(np.min(self.Sf_left[p]))
            maximums_right.append(np.max(self.Sf_right[p]))
            minimums_right.append(np.min(self.Sf_right[p]))

            y__left_avg = y__left_avg + self.Sf_left[p]
            y__right_avg = y__right_avg + self.Sf_right[p]
            velocitys_avg = velocitys_avg + velocitys

        print("maximums_left", maximums_left, "minimums_left", minimums_left, "maximums_right", maximums_right, "minimums_right", minimums_right, "\n")
        maximums_left = np.min(maximums_left)
        minimums_left = np.max(minimums_left)
        maximums_right = np.min(maximums_right)
        minimums_right = np.max(minimums_right)
        print("maximums_left", maximums_left, "minimums_left", minimums_left, "maximums_right", maximums_right, "minimums_right", minimums_right, "\n")

        velocitys_avg = velocitys_avg / len(self.ScanPairs)
        y__left_avg = y__left_avg / len(self.ScanPairs)
        y__right_avg = y__right_avg / len(self.ScanPairs)

        index_min_left = (np.abs(y__left_avg - minimums_left)).argmin()
        index_max_left = (np.abs(y__left_avg - maximums_left)).argmin()
        index_min_right = (np.abs(y__right_avg - minimums_right)).argmin()
        index_max_right = (np.abs(y__right_avg - maximums_right)).argmin()

        print("index", index_min_left, index_max_left, index_min_right, index_max_right)

        velocitys_avg_left = velocitys_avg[index_min_left:index_max_left]
        y__left_avg = y__left_avg[index_min_left:index_max_left]

        velocitys_avg_right = velocitys_avg[index_min_right:index_max_right]
        y__right_avg = y__right_avg[index_min_right:index_max_right]

        self.plot_velocity__left = Plot()
        self.plot_velocity__left.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Left Polarization",(1, 0), "linear")
        self.plot_velocity__left.plot(velocitys_avg_left, y__left_avg, 'b')

        self.plot_velocity__right = Plot()
        self.plot_velocity__right.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Right Polarization",(1, 1), "linear")
        self.plot_velocity__right.plot(velocitys_avg_right, y__right_avg, 'b')

        self.grid.addWidget(self.plot_velocity__left, 0, 0)
        self.grid.addWidget(self.plot_velocity__right, 0, 1)

        totalResults = np.transpose(np.array([np.transpose(velocitys_avg), np.transpose(y__left_avg), np.transpose(y__right_avg)]))

        # source_day_Month_year_houre:minute:righte_station_iteration.dat
        day = scan_1["date"].split("-")[2][0:2]
        month = scan_1["date"].split("-")[1]
        months = {"Jan": "1", "Feb": "2", "Mar": "3", "Apr": "4", "May": "5", "Jun": "6", "Jul": "7", "Aug": "8", "Sep": "9", "Oct": "10", "Nov": "11", "Dec": "12"}
        month = list(months.keys())[int(month) -1]
        year = scan_1["date"].split("-")[0]
        houre = scan_1["date"].split("T")[1].split(":")[0]
        minute = scan_1["date"].split("T")[1].split(":")[1]
        righte = scan_1["date"].split("T")[1].split(":")[2]
        if not os.path.exists(getConfigs("paths", "dataFilePath")  + "SDR/" + getArgs("line") + "/"):
            os.makedirs(getConfigs("paths", "dataFilePath")  + "SDR/" + getArgs("line") + "/")

        output_file_name =  getConfigs("paths", "dataFilePath")  + "SDR/" + getArgs("line") + "/" + getArgs("source") + "_" + day + "_" + month + "_" + year + "_"  + houre + ":" + minute + ":" + righte + "_" + station + "_" + getArgs("iteration_number") + ".dat"
        print("output_file_name", output_file_name)
        output_file_name = output_file_name.replace(" ", "")
        result = Result(totalResults, specie)
        pickle.dump(result, open(output_file_name, 'wb'))

    def nextPair(self):
        if self.index == len(self.ScanPairs)- 1:
            pass

        else:
            self.plot_start__leftA.removePolt()
            self.plot_start__rightB.removePolt()
            self.index = self.index + 1
            self.plotPair(self.index)

    def skipAll(self):
        while self.index < len(self.ScanPairs):
            pair = self.ScanPairs[self.index]
            file1 = self.DataDir + self.__getDataFileForScan__(pair[0][0]) #r0
            file2 = self.DataDir + self.__getDataFileForScan__(pair[0][1]) #s0
            file3 = self.DataDir + self.__getDataFileForScan__(pair[1][0]) #r1
            file4 = self.DataDir + self.__getDataFileForScan__(pair[1][1]) #s1

            frequencyA = self.__getData(file1)[0] #r0
            p_sig_left = self.__getData(file1)[1] #r0
            p_sig_right  = self.__getData(file1)[2] #r0

            frequencyB = self.__getData(file2)[0] #s0
            p_ref_left = self.__getData(file2)[1] #s0
            p_ref_right = self.__getData(file2)[2] #s0

            frequencyC = self.__getData(file3)[0] #r1
            p_sig_on_left = self.__getData(file3)[1] #r1
            p_sig_on_right = self.__getData(file3)[2] #r1

            frequencyD = self.__getData(file4)[0] #s1
            p_ref_on_left = self.__getData(file4)[1] #s1
            p_ref_on_right = self.__getData(file4)[2] #s1

            # fft shift
            p_sig_left = np.fft.fftshift(p_sig_left) #r0
            p_sig_right = np.fft.fftshift(p_sig_right) #r0
            p_ref_left = np.fft.fftshift(p_ref_left) #s0
            p_ref_right = np.fft.fftshift(p_ref_right) #s0
            p_sig_on_left = np.fft.fftshift(p_sig_on_left) #r1
            p_sig_on_right = np.fft.fftshift(p_sig_on_right) #r1
            p_ref_on_left = np.fft.fftshift(p_ref_on_left) #s1
            p_ref_on_right = np.fft.fftshift(p_ref_on_right) #s1

            df_div = float(self.logs["header"]["df_div,df"][0])
            BW = float(self.logs["header"]["Fs,Ns,RBW"][0])
            f_shift = BW / df_div
            l_spec = len(frequencyA)
            f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
            n_shift = int(np.rint(f_shift/f_step))
            avg_interval = 0.5 # inner 50%
            si = int(l_spec / 2 - l_spec * avg_interval / 2)
            ei = int(l_spec / 2 + l_spec * avg_interval / 2)

            Tsys_off_1_left = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_left + p_ref_left) - np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei])) / (2 * np.mean(p_ref_on_left[si:ei] - p_ref_left[si:ei]))
            Tsys_off_2_left = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_left + p_sig_left) - np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei])) / (2 * np.mean(p_sig_on_left[si:ei] - p_sig_left[si:ei]))

            Tsys_off_1_right = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_right + p_ref_right) - np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei])) / (2 * np.mean(p_ref_on_right[si:ei] - p_ref_right[si:ei]))
            Tsys_off_2_right = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_right + p_sig_right ) - np.mean(p_sig_on_right[si:ei] - p_sig_right [si:ei])) / (2 * np.mean(p_sig_on_right[si:ei] - p_sig_right [si:ei]))

            Ta_1_caloff_left = Tsys_off_1_left * (p_sig_left - p_ref_left) / p_ref_left # non-cal phase
            Ta_1_caloff_right = Tsys_off_1_right * (p_sig_right - p_ref_right) / p_ref_right  # non-cal phase

            Ta_1_calon_left = (Tsys_off_1_left + float(self.logs["header"]["Tcal"][0])) * (p_sig_on_left - p_ref_on_left) / p_ref_on_left  # cal phase
            Ta_1_calon_right = (Tsys_off_1_right + float(self.logs["header"]["Tcal"][1])) * (p_sig_on_right - p_ref_on_right) / p_ref_on_right  # cal phase

            Ta_sig_left = (Ta_1_caloff_left + Ta_1_calon_left) / 2
            Ta_sig_right = (Ta_1_caloff_right + Ta_1_calon_right) / 2

            Ta_2_caloff_left = Tsys_off_2_left * (p_sig_left - p_sig_left) / p_sig_left  # non-cal phase
            Ta_2_caloff_right = Tsys_off_2_right * (p_sig_right - p_sig_right ) / p_sig_right   # non-cal phase

            Ta_2_calon_left = (Tsys_off_2_left + float(self.logs["header"]["Tcal"][0])) * (p_ref_on_left - p_sig_on_left) / p_sig_on_left  # cal phase
            Ta_2_calon_right = (Tsys_off_2_right + float(self.logs["header"]["Tcal"][1])) * (p_ref_on_right - p_sig_on_right) / p_sig_on_right  # cal phase

            Ta_ref_left = (Ta_2_caloff_left + Ta_2_calon_left) / 2
            Ta_ref_right = (Ta_2_caloff_right + Ta_2_calon_right) / 2

            Ta_sig_left = np.roll(Ta_sig_left, +n_shift)
            Ta_sig_right = np.roll(Ta_sig_right, +n_shift)

            Ta_ref_left = np.roll(Ta_ref_left, -n_shift)
            Ta_ref_right = np.roll(Ta_ref_right, -n_shift)

            Ta_left = (Ta_sig_left + Ta_ref_left) / 2
            Ta_right = (Ta_sig_right + Ta_ref_right) / 2

            El = (float(self.logs[pair[0][0]]["AzEl"][1]) + float(self.logs[pair[0][1]]["AzEl"][1]) + float(self.logs[pair[1][0]]["AzEl"][1]) + float(self.logs[pair[1][1]]["AzEl"][1]))/ 4
            G_El = self.logs["header"]["Elev_poly"]
            G_El = [float(gel) for gel in G_El]
            G_ELtmp = [0,0,0]
            G_ELtmp[0] = G_El[2]
            G_ELtmp[1] = G_El[1]
            G_ELtmp[2] = G_El[0]
            G_El = G_ELtmp

            Sf_leftscan = Ta_left / (((-1) * (float(self.logs["header"]["DPFU"][0])) ) * np.polyval(G_El, El))
            Sf_rightscan = Ta_right / (((-1) * (float(self.logs["header"]["DPFU"][1])) ) * np.polyval(G_El, El))

            self.Sf_left.append(Sf_leftscan[self.si:self.ei])
            self.Sf_right.append(Sf_rightscan[self.si:self.ei])
            self.x = frequencyA[self.si:self.ei]

            self.si = si
            self.ei = ei

            if self.index == len(self.ScanPairs) - 1:
                self.nextPairButton.setText('Move to total results')
                self.nextPairButton.clicked.connect(self.plotTotalResults)

            self.index +=1

        self.plotTotalResults()

    def __UI__(self):
        if self.index != len(self.ScanPairs) - 1:  # cheking if there is not one pair
            self.nextPairButton = QPushButton("Next pair", self)
            self.nextPairButton.clicked.connect(self.nextPair)
            self.grid.addWidget(self.nextPairButton, 4, 2)

        self.skipAllButton = QPushButton("Skip to end", self)
        self.skipAllButton.clicked.connect(self.skipAll)
        self.grid.addWidget(self.skipAllButton, 5, 2)

        self.plotPair(self.index)

def main():
    qApp = QApplication(sys.argv)
    a = Analyzer()
    a.show()
    sys.exit(qApp.exec_())
    sys.exit(0)

if __name__ == "__main__":
    main()
