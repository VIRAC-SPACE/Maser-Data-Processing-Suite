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
        self.DataDir = getConfigs("paths", "dataFilePath")  + "SDR/" + getArgs("source") + "/f" + getArgs("line")  + "/" + getArgs("iteration_number") + "/"
        self.DataFiles = os.listdir(self.DataDir)
        self.ScanPairs = self.createScanPairs()
        self.index = 0
        self.Sf_first = list()
        self.Sf_second = list()
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
        polarization_first = correctNumpyReadData(data[:, [1]])
        polarization_second = correctNumpyReadData(data[:, [2]])
        return(frequency, polarization_first, polarization_second)

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
        p_sig_first = self.__getData(file1)[1] #r0
        p_sig_second  = self.__getData(file1)[2] #r0

        frequencyB = self.__getData(file2)[0] #s0
        p_ref_first = self.__getData(file2)[1] #s0
        p_ref_second = self.__getData(file2)[2] #s0

        frequencyC = self.__getData(file3)[0] #r1
        p_sig_on_first = self.__getData(file3)[1] #r1
        p_sig_on_second = self.__getData(file3)[2] #r1

        frequencyD = self.__getData(file4)[0] #s1
        p_ref_on_first = self.__getData(file4)[1] #s1
        p_ref_on_second = self.__getData(file4)[2] #s1

        # fft shift
        p_sig_first = np.fft.fftshift(p_sig_first) #r0
        p_sig_second  = np.fft.fftshift(p_sig_second) #r0
        p_ref_first = np.fft.fftshift(p_ref_first) #s0
        p_ref_second = np.fft.fftshift(p_ref_second) #s0
        p_sig_on_first = np.fft.fftshift(p_sig_on_first) #r1
        p_sig_on_second = np.fft.fftshift(p_sig_on_second) #r1
        p_ref_on_first = np.fft.fftshift(p_ref_on_first) #s1
        p_ref_on_second = np.fft.fftshift(p_ref_on_second) #s1

        #plot1
        self.plot_start__firstA = Plot()
        self.plot_start__firstA.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "First Polarization", (1, 0), "linear")
        self.plot_start__firstA.plot(frequencyA, p_sig_first, 'b', label=str(index+1) + "r0")
        self.plot_start__firstA.plot(frequencyB, p_ref_first, 'g', label=str(index+1) + "s0")
        self.plot_start__firstA.plot(frequencyC, p_sig_on_first, 'r', label=str(index+1) + "r1")
        self.plot_start__firstA.plot(frequencyD, p_ref_on_first, 'y', label=str(index+1) + "s1")
        self.grid.addWidget(self.plot_start__firstA, 0, 0)

        #plot2
        self.plot_start__secondB = Plot()
        self.plot_start__secondB.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "Second Polarization", (1, 1), "linear")
        self.plot_start__secondB.plot(frequencyA, p_sig_second , 'b', label=str(index+1) + "r0")
        self.plot_start__secondB.plot(frequencyB, p_ref_second, 'g', label=str(index+1) + "s0")
        self.plot_start__secondB.plot(frequencyC, p_sig_on_second, 'r', label=str(index+1) + "r1")
        self.plot_start__secondB.plot(frequencyD, p_ref_on_second, 'y', label=str(index+1) + "s1")
        self.grid.addWidget(self.plot_start__secondB, 0, 1)

        df_div = float(self.logs["header"]["df_div,df"][0])
        BW = float(self.logs["header"]["Fs,Ns,RBW"][0])
        f_shift = BW / df_div
        l_spec = len(frequencyA)
        f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
        n_shift = int(np.rint(f_shift / f_step))
        avg_interval = 0.5 # inner 50%
        si = int(l_spec / 2 - l_spec * avg_interval / 2)
        ei = int(l_spec / 2 + l_spec * avg_interval / 2)

        Tsys_off_1_first = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_first + p_ref_first) - np.mean(p_ref_on_first[si:ei] - p_ref_first[si:ei])) / (2 * np.mean(p_ref_on_first[si:ei] - p_ref_first[si:ei]))
        Tsys_off_2_first = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_first + p_sig_first) - np.mean(p_sig_on_first[si:ei] - p_sig_first[si:ei])) / (2 * np.mean(p_sig_on_first[si:ei] - p_sig_first[si:ei]))

        Tsys_off_1_second = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_second + p_ref_second) - np.mean(p_ref_on_second[si:ei] - p_ref_second[si:ei])) / (2 * np.mean(p_ref_on_second[si:ei] - p_ref_second[si:ei]))
        Tsys_off_2_second = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_second + p_sig_second ) - np.mean(p_sig_on_second[si:ei] - p_sig_second [si:ei])) / (2 * np.mean(p_sig_on_second[si:ei] - p_sig_second [si:ei]))

        Ta_1_caloff_first = Tsys_off_1_first * (p_sig_first - p_ref_first) / p_ref_first # non-cal phase
        Ta_1_caloff_second = Tsys_off_1_second * (p_sig_second - p_ref_second) / p_ref_second  # non-cal phase

        Ta_1_calon_first = (Tsys_off_1_first + float(self.logs["header"]["Tcal"][0])) * (p_sig_on_first - p_ref_on_first) / p_ref_on_first  # cal phase
        Ta_1_calon_second = (Tsys_off_1_second + float(self.logs["header"]["Tcal"][1])) * (p_sig_on_second - p_ref_on_second) / p_ref_on_second  # cal phase

        Ta_sig_first = (Ta_1_caloff_first + Ta_1_calon_first) / 2
        Ta_sig_second = (Ta_1_caloff_second + Ta_1_calon_second) / 2

        Ta_2_caloff_first = Tsys_off_2_first * (p_ref_first - p_sig_first) / p_sig_first  # non-cal phase
        Ta_2_caloff_second = Tsys_off_2_second * (p_ref_second - p_sig_second ) / p_sig_second   # non-cal phase

        Ta_2_calon_first = (Tsys_off_2_first + float(self.logs["header"]["Tcal"][0])) * (p_ref_on_first - p_sig_on_first) / p_sig_on_first  # cal phase
        Ta_2_calon_second = (Tsys_off_2_second + float(self.logs["header"]["Tcal"][1])) * (p_ref_on_second - p_sig_on_second) / p_sig_on_second  # cal phase

        Ta_ref_first = (Ta_2_caloff_first + Ta_2_calon_first) / 2
        Ta_ref_second = (Ta_2_caloff_second + Ta_2_calon_second) / 2

        Ta_sig_first = np.roll(Ta_sig_first, +n_shift)
        Ta_sig_second = np.roll(Ta_sig_second, +n_shift)

        Ta_ref_first = np.roll(Ta_ref_first, -n_shift)
        Ta_ref_second = np.roll(Ta_ref_second, -n_shift)

        Ta_first = (Ta_sig_first + Ta_ref_first) / 2
        Ta_second = (Ta_sig_second + Ta_ref_second) / 2

        El = (float(self.logs[pair[0][0]]["AzEl"][1]) + float(self.logs[pair[0][1]]["AzEl"][1]) + float(self.logs[pair[1][0]]["AzEl"][1]) + float(self.logs[pair[1][1]]["AzEl"][1]))/ 4
        G_El = self.logs["header"]["Elev_poly"]
        G_El = [float(gel) for gel in G_El]
        G_ELtmp = [0,0,0]
        G_ELtmp[0] = G_El[2]
        G_ELtmp[1] = G_El[1]
        G_ELtmp[2] = G_El[0]
        G_El = G_ELtmp

        Sf_firstscan = Ta_first / (((-1) * float(self.logs["header"]["DPFU"][0])) * np.polyval(G_El, El))
        Sf_secondscan = Ta_second / (((-1) * float(self.logs["header"]["DPFU"][1])) * np.polyval(G_El, El))

        self.si = si
        self.ei = ei

        self.Sf_first.append(Sf_firstscan[self.si:self.ei])
        self.Sf_second.append(Sf_secondscan[self.si:self.ei])

        # plot3
        self.total__first = Plot()
        self.total__first.creatPlot(self.grid, 'Frequency Mhz', 'Amplitude', "", (4, 0), "linear")
        self.x = frequencyA[self.si:self.ei]
        self.total__first.plot(frequencyA[self.si:self.ei], Sf_firstscan[self.si:self.ei], 'b', label=str(index + 1))
        self.grid.addWidget(self.total__first, 3, 0)

        # plot4
        self.total__second = Plot()
        self.total__second.creatPlot(self.grid, 'Frequency Mhz', 'Flux density (Jy)', "", (4, 1), "linear")
        self.total__second.plot(frequencyA[self.si:self.ei], Sf_secondscan[self.si:self.ei], 'b', label=str(index + 1))
        self.grid.addWidget(self.total__second, 3, 1)

        if index == len(self.ScanPairs) - 1:
            self.nextPairButton.setText('Move to total results')
            self.nextPairButton.clicked.connect(self.plotTotalResults)

    def plotTotalResults(self):
        self.grid.removeWidget(self.plot_start__firstA)
        self.grid.removeWidget(self.plot_start__secondB)
        self.grid.removeWidget(self.total__first)
        self.grid.removeWidget(self.total__second)

        self.plot_start__firstA.hide()
        self.plot_start__secondB.hide()
        self.total__first.hide()
        self.total__second.hide()

        self.plot_start__firstA.close()
        self.plot_start__secondB.close()
        self.total__first.close()
        self.total__second.close()

        self.plot_start__firstA.removePolt()
        self.plot_start__secondB.removePolt()
        self.total__first.removePolt()
        self.total__second.removePolt()

        del self.plot_start__firstA
        del self.plot_start__secondB
        del self.total__first
        del self.total__second

        self.grid.removeWidget(self.nextPairButton)
        self.nextPairButton.hide()
        self.nextPairButton.close()
        del self.nextPairButton

        for i in reversed(range(self.grid.count())):
            self.grid.itemAt(i).widget().deleteLater()

        velocitys_avg = np.zeros(len(self.x))
        y__first_avg = np.zeros(len(self.x))
        y__second_avg = np.zeros(len(self.x))

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

            self.max_y_first_index = self.Sf_first[p].argmax(axis=0)
            self.max_y_second_index = self.Sf_second[p].argmax(axis=0)

            line = getConfigs('base_frequencies_SDR', "f" + getArgs("line")).replace(" ", "").split(",")
            lineF = float(line[0]) * (10**9)
            lineS = line[1]
            specie = lineS

            print("specie", specie, "\n")
            LO = float(self.logs["header"]["f_obs,LO,IF"][1])
            velocitys = dopler((self.x + LO) * (10 ** 6), VelTotal, lineF)
            y__first_avg = y__first_avg + self.Sf_first[p]
            y__second_avg = y__second_avg + self.Sf_second[p]
            velocitys_avg = velocitys_avg + velocitys

        velocitys_avg = velocitys_avg / len(self.ScanPairs)
        y__first_avg = y__first_avg / len(self.ScanPairs)
        y__second_avg = y__second_avg / len(self.ScanPairs)

        self.plot_velocity__first = Plot()
        self.plot_velocity__first.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "First Polarization",(1, 0), "linear")
        self.plot_velocity__first.plot(velocitys_avg, y__first_avg, 'b')

        self.plot_velocity__second = Plot()
        self.plot_velocity__second.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Second Polarization",(1, 1), "linear")
        self.plot_velocity__second.plot(velocitys_avg, y__second_avg, 'b')

        self.grid.addWidget(self.plot_velocity__first, 0, 0)
        self.grid.addWidget(self.plot_velocity__second, 0, 1)

        totalResults = np.transpose(np.array([np.transpose(velocitys_avg), np.transpose(y__first_avg), np.transpose(y__second_avg)]))

        # source_day_Month_year_houre:minute:seconde_station_iteration.dat
        day = scan_1["date"].split("-")[2][0:2]
        month = scan_1["date"].split("-")[1]
        months = {"Jan": "1", "Feb": "2", "Mar": "3", "Apr": "4", "May": "5", "Jun": "6", "Jul": "7", "Aug": "8", "Sep": "9", "Oct": "10", "Nov": "11", "Dec": "12"}
        month = list(months.keys())[int(month) -1]
        year = scan_1["date"].split("-")[0]
        houre = scan_1["date"].split("T")[1].split(":")[0]
        minute = scan_1["date"].split("T")[1].split(":")[1]
        seconde = scan_1["date"].split("T")[1].split(":")[2]
        output_file_name =  getConfigs("paths", "dataFilePath")  + "SDR/" + getArgs("source") + "_" + day + "_" + month + "_" + year + "_"  + houre + ":" + minute + ":" + seconde + "_" + station + "_" + getArgs("iteration_number") + ".dat"
        print("output_file_name", output_file_name)
        output_file_name = output_file_name.replace(" ", "")
        result = Result(totalResults, specie)
        pickle.dump(result, open(output_file_name, 'wb'))

    def nextPair(self):
        if self.index == len(self.ScanPairs)- 1:
            pass

        else:
            self.plot_start__firstA.removePolt()
            self.plot_start__secondB.removePolt()
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
            p_sig_first = self.__getData(file1)[1] #r0
            p_sig_second  = self.__getData(file1)[2] #r0

            frequencyB = self.__getData(file2)[0] #s0
            p_ref_first = self.__getData(file2)[1] #s0
            p_ref_second = self.__getData(file2)[2] #s0

            frequencyC = self.__getData(file3)[0] #r1
            p_sig_on_first = self.__getData(file3)[1] #r1
            p_sig_on_second = self.__getData(file3)[2] #r1

            frequencyD = self.__getData(file4)[0] #s1
            p_ref_on_first = self.__getData(file4)[1] #s1
            p_ref_on_second = self.__getData(file4)[2] #s1

            # fft shift
            p_sig_first = np.fft.fftshift(p_sig_first) #r0
            p_sig_second = np.fft.fftshift(p_sig_second) #r0
            p_ref_first = np.fft.fftshift(p_ref_first) #s0
            p_ref_second = np.fft.fftshift(p_ref_second) #s0
            p_sig_on_first = np.fft.fftshift(p_sig_on_first) #r1
            p_sig_on_second = np.fft.fftshift(p_sig_on_second) #r1
            p_ref_on_first = np.fft.fftshift(p_ref_on_first) #s1
            p_ref_on_second = np.fft.fftshift(p_ref_on_second) #s1

            df_div = float(self.logs["header"]["df_div,df"][0])
            BW = float(self.logs["header"]["Fs,Ns,RBW"][0])
            f_shift = BW / df_div
            l_spec = len(frequencyA)
            f_step = (frequencyA[l_spec - 1] - frequencyA[0]) / (l_spec - 1)
            n_shift = int(np.rint(f_shift/f_step))
            avg_interval = 0.5 # inner 50%
            si = int(l_spec / 2 - l_spec * avg_interval / 2)
            ei = int(l_spec / 2 + l_spec * avg_interval / 2)

            Tsys_off_1_first = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_first + p_ref_first) - np.mean(p_ref_on_first[si:ei] - p_ref_first[si:ei])) / (2 * np.mean(p_ref_on_first[si:ei] - p_ref_first[si:ei]))
            Tsys_off_2_first = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_first + p_sig_first) - np.mean(p_sig_on_first[si:ei] - p_sig_first[si:ei])) / (2 * np.mean(p_sig_on_first[si:ei] - p_sig_first[si:ei]))

            Tsys_off_1_second = float(self.logs["header"]["Tcal"][0]) * ((p_ref_on_second + p_ref_second) - np.mean(p_ref_on_second[si:ei] - p_ref_second[si:ei])) / (2 * np.mean(p_ref_on_second[si:ei] - p_ref_second[si:ei]))
            Tsys_off_2_second = float(self.logs["header"]["Tcal"][1]) * ((p_sig_on_second + p_sig_second ) - np.mean(p_sig_on_second[si:ei] - p_sig_second [si:ei])) / (2 * np.mean(p_sig_on_second[si:ei] - p_sig_second [si:ei]))

            Ta_1_caloff_first = Tsys_off_1_first * (p_sig_first - p_ref_first) / p_ref_first # non-cal phase
            Ta_1_caloff_second = Tsys_off_1_second * (p_sig_second - p_ref_second) / p_ref_second  # non-cal phase

            Ta_1_calon_first = (Tsys_off_1_first + float(self.logs["header"]["Tcal"][0])) * (p_sig_on_first - p_ref_on_first) / p_ref_on_first  # cal phase
            Ta_1_calon_second = (Tsys_off_1_second + float(self.logs["header"]["Tcal"][1])) * (p_sig_on_second - p_ref_on_second) / p_ref_on_second  # cal phase

            Ta_sig_first = (Ta_1_caloff_first + Ta_1_calon_first) / 2
            Ta_sig_second = (Ta_1_caloff_second + Ta_1_calon_second) / 2

            Ta_2_caloff_first = Tsys_off_2_first * (p_sig_first - p_sig_first) / p_sig_first  # non-cal phase
            Ta_2_caloff_second = Tsys_off_2_second * (p_sig_second - p_sig_second ) / p_sig_second   # non-cal phase

            Ta_2_calon_first = (Tsys_off_2_first + float(self.logs["header"]["Tcal"][0])) * (p_ref_on_first - p_sig_on_first) / p_sig_on_first  # cal phase
            Ta_2_calon_second = (Tsys_off_2_second + float(self.logs["header"]["Tcal"][1])) * (p_ref_on_second - p_sig_on_second) / p_sig_on_second  # cal phase

            Ta_ref_first = (Ta_2_caloff_first + Ta_2_calon_first) / 2
            Ta_ref_second = (Ta_2_caloff_second + Ta_2_calon_second) / 2

            Ta_sig_first = np.roll(Ta_sig_first, +n_shift)
            Ta_sig_second = np.roll(Ta_sig_second, +n_shift)

            Ta_ref_first = np.roll(Ta_ref_first, -n_shift)
            Ta_ref_second = np.roll(Ta_ref_second, -n_shift)

            Ta_first = (Ta_sig_first + Ta_ref_first) / 2
            Ta_second = (Ta_sig_second + Ta_ref_second) / 2

            El = (float(self.logs[pair[0][0]]["AzEl"][1]) + float(self.logs[pair[0][1]]["AzEl"][1]) + float(self.logs[pair[1][0]]["AzEl"][1]) + float(self.logs[pair[1][1]]["AzEl"][1]))/ 4
            G_El = self.logs["header"]["Elev_poly"]
            G_El = [float(gel) for gel in G_El]
            G_ELtmp = [0,0,0]
            G_ELtmp[0] = G_El[2]
            G_ELtmp[1] = G_El[1]
            G_ELtmp[2] = G_El[0]
            G_El = G_ELtmp

            Sf_firstscan = Ta_first / (((-1) * (float(self.logs["header"]["DPFU"][0])) ) * np.polyval(G_El, El))
            Sf_secondscan = Ta_second / (((-1) * (float(self.logs["header"]["DPFU"][1])) ) * np.polyval(G_El, El))

            self.Sf_first.append(Sf_firstscan[self.si:self.ei])
            self.Sf_second.append(Sf_secondscan[self.si:self.ei])
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
