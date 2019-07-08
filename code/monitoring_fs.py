#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import json
import numpy as np
from numpy import genfromtxt
from scipy.interpolate import griddata
import matplotlib
from matplotlib.ticker import MaxNLocator
from astropy.time import Time
from astropy.stats import LombScargle
from astropy.modeling import models, fitting
from astropy.modeling.fitting import LevMarLSQFitter
import datetime
from functools import reduce
from operator import itemgetter
from multiprocessing import Process
import math

from PyQt5.QtWidgets import QWidget, QGridLayout, QApplication, QPushButton, QLabel, QLineEdit, QDesktopWidget, QComboBox, QGroupBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

from ploting_qt5 import  Plot
from parsers._configparser import ConfigParser
from result import  Result
from monitor.months import Months
from monitor.monitoringViewHelper import MonitoringViewHelper
from help import *

matplotlib.use('Qt5Agg')


def parseArguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')
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


class PlottingView(QWidget):
    __slots__ = ['grid']
    def __init__(self):

        super().__init__()
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
                
    def _addWidget(self, widget, row, collon):
        self.grid.addWidget(widget, row, collon)
            
    def getGrid(self):
        return self.grid

    
class Spectre_View(PlottingView):
        def __init__(self):
            super().__init__()
            self.setWindowTitle(" ")


class Gauss_View2(PlottingView):
    __slots__ = ['AreaList', 'source', 'GaussDatePointsList', 'gaussLocationList', 'velocitie_dict_componet', "date_list"]
    def __init__(self, AreaList, source, GaussDatePointsList, gaussLocationList, gaussIterationList, velocitie_dict_componet, date_list):

        super().__init__()
        self.setWindowTitle(" ")
        self.AreaList = AreaList
        self.source = source
        self.GaussDatePointsList = GaussDatePointsList
        self.gaussLocationList = gaussLocationList
        self.gaussIterationList = gaussIterationList
        self.velocitie_dict_componet = velocitie_dict_componet
        self.date_list = date_list

    def plot(self):
        gaussLines = getConfigs("gauss_lines", self.source).replace(" ", "").split(",")
        gaussLineDict = dict()

        for l in gaussLines:
            gaussLineDict[l] = list()

        for areas in self.AreaList:
            i = 0
            for g in gaussLines:
                gaussLineDict[g].append(areas[i])
                i += 1

        self.gaussAreaPlot = Plot()
        self.gaussAreaPlot.creatPlot(self.getGrid(), "Time", "Gauss Amplitude", None, (1, 0), "linear")

        i = 0
        colors = []
        r = 0.25
        g = 0.2
        b = 0.15
        a = 1

        for index in range(0, len(gaussLines)):
            r += 0.08
            g += 0.01
            b += 0.02
            colors.append((r,g,b,a))

        labels2 = ["Date is " + d.strftime('%d %m %Y') for d in self.GaussDatePointsList]

        for gl in (gaussLines):
            self.gaussAreaPlot.plot(self.GaussDatePointsList, gaussLineDict[gl], ".", color=colors[i], fontsize=8, visible=True, picker=5, label=gl)

            i += 1

        self.gaussAreaPlot.addPickEvent(self.chooseGaussFitt)
        self._addWidget(self.gaussAreaPlot, 0, 0)
        self.gaussAreaPlot.addCursor(labels2)

        self.testView = PlottingView()
        self.testPlot = Plot()
        self.testPlot.creatPlot(self.testView.getGrid(), "Time", "Gauss Amplitude", None, (1, 0), "linear")
        self.testPlot.plot(self.date_list, self.velocitie_dict_componet, "r.", label="Orginal data")
        self.testPlot.plot(self.GaussDatePointsList, gaussLineDict["-1.77"],  "g.", label="Gauss model")
        self.testView._addWidget(self.testPlot, 0, 0)
        self.testView.show()

    def chooseGaussFitt(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ind = event.ind
        index = [ind][0]
        fittFile = getConfigs("paths", "notSmoohtFilePath") + self.source + "/" + self.source + "_" + MonitoringViewHelper.formatDate(xdata, index) + "_" + MonitoringViewHelper.getLocation(self.gaussLocationList, int(index[0])) + "_"  + MonitoringViewHelper.getIteration(self.gaussIterationList, int(index[0])) + ".dat"
        gaussLines = getConfigs("gauss_lines", self.source).replace(" ", "").split(",")
        data = np.fromfile(fittFile, dtype="float64", count=-1, sep=" ").reshape((file_len(fittFile), 4))
        velocity = correctNumpyReadData(data[:, [0]])
        ampvid = correctNumpyReadData(data[:, [3]])
        indexies = [(np.abs(velocity - float(line))).argmin() for line in gaussLines]
        mons = [max(ampvid[index - 5:index + 5]) for index in indexies]
        gaussian = [models.Gaussian1D(mons[index], gaussLines[index], 0.05, bounds={'stddev': (None, 0.15)}) for index in range(0, len(mons))]

        def sum(a, b):
            return a + b

        gg_init = reduce(sum, gaussian)
        fit = LevMarLSQFitter()
        gg_fit = fit(gg_init, velocity, ampvid)
        sts = [models.Gaussian1D(gg_fit[index].amplitude, gg_fit[index].mean, gg_fit[index].stddev) for index in range(0, len(gaussLines))]

        self.fitPlotview = PlottingView()
        self.fitPlot = Plot()
        self.fitPlot.creatPlot(self.fitPlotview .getGrid(), "Velocity", "Gauss Fits", MonitoringViewHelper.formatDate(xdata, index), (1, 0), "linear")
        self.fitPlot.plot(velocity, ampvid, 'C0+', label="data")
        self.fitPlot.plot(velocity, gg_fit(velocity), "C1-", label="total fit", linewidth=4)

        colors = []

        for index in range(2, len(gaussLines) + 2):
            if index % 2:
                index -= 1
            else:
                index += 1
            colors.append("C" + str(index))

        c = 0
        for st in sts:
            self.fitPlot.plot(velocity, st(velocity), colors[c], label=gaussLines[c])
            c += 1

        self.fitPlot.plot(velocity, (gg_fit(velocity) - ampvid) * 1 - 35, "r")
        self.fitPlot.plot(velocity, velocity * 0 - 35, "g")
        self.fitPlotview._addWidget(self.fitPlot, 0, 0)
        self.fitPlotview.show()


class Gauss_View(PlottingView):
    __slots__ = ['AreaList', 'source', 'GaussDatePointsList', 'gaussLocationList']

    def __init__(self, AreaList, source, GaussDatePointsList, gaussLocationList, gaussIterationList):
        super().__init__()
        self.setWindowTitle(" ")
        self.AreaList = AreaList
        self.source = source
        self.GaussDatePointsList = GaussDatePointsList
        self.gaussLocationList = gaussLocationList
        self.gaussIterationList = gaussIterationList

    def plot(self):
        gaussLines = getConfigs("gauss_lines", self.source).replace(" ", "").split(",")
        gaussLineDict = dict()

        for l in gaussLines:
            gaussLineDict[l] = list()

        for areas in self.AreaList:
            i = 0
            for g in gaussLines:
                gaussLineDict[g].append(areas[i])
                i += 1

        self.gaussAreaPlot = Plot()
        self.gaussAreaPlot.creatPlot(self.getGrid(), "Time", "Gauss Areas", None, (1, 0), "linear")

        i = 0
        colors = []
        r = 0.25
        g = 0.2
        b = 0.15
        a = 1

        for index in range(0, len(gaussLines)):
            r += 0.08
            g += 0.01
            b += 0.02
            colors.append((r, g, b, a))

        labels2 = ["Date is " + d.strftime('%d %m %Y') for d in self.GaussDatePointsList]
        for gl in (gaussLines):
            self.gaussAreaPlot.plot(self.GaussDatePointsList, gaussLineDict[gl], ".", fontsize=8, color=colors[i], visible=True, picker=5, label=gl)
            i += 1

        self.gaussAreaPlot.addPickEvent(self.chooseGaussFitt)
        self._addWidget(self.gaussAreaPlot, 0, 0)
        self.gaussAreaPlot.addCursor(labels2)

    def chooseGaussFitt(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ind = event.ind
        index = [ind][0]
        gaussLines = getConfigs("gauss_lines", self.source).replace(" ", "").split(",")
        fittFile = getConfigs("paths", "notSmoohtFilePath") + self.source + "/" + self.source + "_" + MonitoringViewHelper.formatDate(xdata, index) + "_" + MonitoringViewHelper.getLocation(self.gaussLocationList, int(index[0])) + "_"  + MonitoringViewHelper.getIteration(self.gaussIterationList, int(index[0])) + ".dat"

        gaussLines = getConfigs("gauss_lines", self.source).replace(" ", "").split(",")

        data = np.fromfile(fittFile, dtype="float64", count=-1, sep=" ").reshape((file_len(fittFile), 4))
        velocity = correctNumpyReadData(data[:, [0]])
        ampvid = correctNumpyReadData(data[:, [3]])
        indexies = [(np.abs(velocity - float(line))).argmin() for line in gaussLines]
        mons = [max(ampvid[index - 5:index + 5]) for index in indexies]
        gaussian = [models.Gaussian1D(mons[index], gaussLines[index], 0.05, bounds={'stddev': (None, 0.15)}) for index in range(0, len(mons))]

        def sum(a, b):
            return a + b

        gg_init = reduce(sum, gaussian)
        fit = LevMarLSQFitter()
        gg_fit = fit(gg_init, velocity, ampvid)
        sts = [models.Gaussian1D(gg_fit[index].amplitude, gg_fit[index].mean, gg_fit[index].stddev) for index in range(0, len(gaussLines))]

        self.fitPlotview = PlottingView()
        self.fitPlot = Plot()
        self.fitPlot.creatPlot(self.fitPlotview .getGrid(), "Velocity", "Gauss Fits", MonitoringViewHelper.formatDate(xdata, index), (1, 0), "linear")
        self.fitPlot.plot(velocity, ampvid, 'C0+', label="data")
        self.fitPlot.plot(velocity, gg_fit(velocity), "C1-", label="total fit", linewidth=4)

        colors = []

        for index in range(2, len(gaussLines) + 2):
            if index % 2:
                index -= 1
            else:
                index += 1
            colors.append("C" + str(index))

        c = 0
        for st in sts:
            self.fitPlot.plot(velocity, st(velocity), colors[c], label=gaussLines[c])
            c += 1

        self.fitPlot.plot(velocity, (gg_fit(velocity) - ampvid) * 1 - 35, "r")
        self.fitPlot.plot(velocity, velocity * 0 - 35, "g")
        self.fitPlotview._addWidget(self.fitPlot, 0, 0)
        self.fitPlotview.show()

class Maps_View(PlottingView):
        __slots__ = ['source', 'label']            
        def __init__(self, source):
            super().__init__()
            self.setWindowTitle(" ")
            self.source = source
            
        def plotMaps(self, output_path):
            parametrs = json.load(open("python.txt", "r"))
            fin = open(output_path + self.source + "/" + "3d.txt","r")
            max_flux_limit = 9999
            JD_shift = 0
            jd_min = 0
            jd_max = float(parametrs["days"])
            vmin = float(parametrs["min_v"])
            vmax = float(parametrs["max_v"]) 
            jd_first = float(parametrs["jd_first"])
            vrange = vmax-vmin
            sources_vrange = genfromtxt('DB_vrange.csv', names=True, delimiter=',', dtype=None, encoding=None)
            
            i = 0
            found_vrange = False
            found_ind = -1
            for source in sources_vrange:
                if(str(source[0]) == self.source):
                    print('We found vrange for source '+source[0]+': from '+str(source[1])+' to '+str(source[2]))
                    found_vrange = True
                    found_ind = i
                    vmin = source[1]
                    vmax = source[2]
                i = i + 1

            velocity = []
            observed_flux = []
            observed_time = []
            
            print('Reading data...\n')
            
            while True:
                line = fin.readline()
                if(len(line)) == 0:
                    break
                if (line.find("#")==-1):
                    line = line.split()
                    velocity.append(float(line[0]))
                    if (float(line[1]) > max_flux_limit):
                        observed_flux.append( max_flux_limit)
                    else:
                        if float(line[1]) < 0:
                            observed_flux.append( 0.001 )
                        else:
                            observed_flux.append( float(line[1])) 
                
                    observed_time.append( float(line[2])+JD_shift )
            #        if (float(line[2]) == 316) and ( float(line[0])== -9.41):
            #            print float(line[1])    
                        
            fin.close()
            velocity = np.array(velocity)
            observed_flux = np.array(observed_flux)
            observed_time = np.array(observed_time)
            
            print('Initialize arrays...\n')
            
            x = np.arange(vmin, vmax, 0.01)
            y = np.arange(0, jd_max, 0.5)
            X, Y = np.meshgrid(x,y)
            
            print('Regridding data...\n')
            #Z = matplotlib.mlab.griddata(velocity,observed_time,observed_flux, X, Y, interp='linear')
            Z = griddata((velocity, observed_time), observed_flux, (X[:], Y[:]), method='linear')
            self.mapPlot = Plot()
            self.mapPlot.creatPlot(self.getGrid(), "Velocity (km/s)", "JD (days) -  " + str(jd_first), None, (1,0), "linear")
            CS = self.mapPlot.contourf(X, Y, Z)
    
            cbar = self.mapPlot.colorbar(CS)
            cbar.set_clim(vmin=0) # ,vmax=max_flux_limit
            cbar.ax.set_ylabel(r'$Flux~(\mathrm{Jy})$')
            cbar.locator = MaxNLocator(nbins = 50)
            
            text_file = open(output_path + self.source + "/" + "labels.txt", "r")
            days = text_file.readline().rstrip().split(',')
            days = [float(i) for i in days]
            dates = text_file.readline().rstrip().split(',')
            dates_days = text_file.readline().rstrip().split(',')
            dates_days = [float(i) for i in dates_days]
            text_file.close()
            
            print ("adding date to plot....")
            for iday in days:
                CS = self.mapPlot.plot([vmax-math.fabs(vmax-vmin)/8.0, vmax], [iday+0.0, iday+0.0], 'w', linewidth=0.5, alpha=0.9)
                
            for i in range(0,len(dates_days)-1):
                iday = dates_days[i]
                self.mapPlot.setAxiesText(vmin+0.2, iday, dates[i], fontsize=5, color='white', horizontalalignment='left', alpha= 1.0)

            self.mapPlot.setAxiesText(0.01,1.01, dates[len(dates_days)-1], fontsize=5, color='black', horizontalalignment='left', alpha= 1.0)
            print ("Done")
            self._addWidget(self.mapPlot, 0, 0)


class Period_View(PlottingView):
        def __init__(self):
            super().__init__()
            self.setWindowTitle("Periods in days ")
            
        def convertDatetimeObjectToMJD(self, time):
            time=time.isoformat()
            t=Time(time, format='isot')
            return t.mjd 
        
        def PlotPeriods(self, dateList, velocitie_dict, source_velocities, velocityIndex, plotSimbol):
            x = dateList
            t = np.array([self.convertDatetimeObjectToMJD(i) for i in x], dtype="float64")
            y =  velocitie_dict["avg"][source_velocities[velocityIndex]]
            error = np.array(y) * 0.1
            ls  = LombScargle(t,  y, error, fit_mean=True)
            
            def dateDelta(d1, d2):
                return abs(d1 -d2)
            
            def getMaxDateDelta():
                maxDelta = list()
                for x in range(0, len(t) -1):
                    maxDelta.append(dateDelta(t[x], t[x + 1]))
    
                return np.max(maxDelta)
            
            def getMinDateDelta():
                minDelta = list()
                for x in range(0, len(t) -1):
                    minDelta.append(dateDelta(t[x], t[x + 1]))
    
                return np.min(minDelta)
            
            nyquist_factor = 2 * getMaxDateDelta()
            maximum_frequency = 2  * getMinDateDelta()
            minimum_frequency = 1/dateDelta(t[0], t[-1])
                
            print("nyquist_factor", nyquist_factor)
            print("minimum_frequency", minimum_frequency)
            print("maximum_frequency", maximum_frequency)
            
            frequency, power = ls.autopower(method='fastchi2', normalization='model', nyquist_factor=nyquist_factor, minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency, samples_per_peak=20)
            false_alarm = ls.false_alarm_probability(power.max(), method="bootstrap", nyquist_factor=nyquist_factor, minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency, samples_per_peak=20)
            
            print("max power", power.max(), "false_alarm", false_alarm)
            
            period_days = 1. / frequency
            best_period = period_days[np.argmax(power)]
            print("Best period: {0:.2f} hours".format(24 * best_period))
            
            self.periodPlot = Plot()
            self.periodPlot.creatPlot(self.getGrid(), "Period (days)", "Power", None, (1,0), "linear")
            self.periodPlot.plot(period_days, power, plotSimbol, label="polarization AVG " + "Velocity " + source_velocities[velocityIndex], rasterized=True)
            self._addWidget(self.periodPlot, 0, 0)
            self.show()


class Monitoring_View(PlottingView):
        def __init__(self, iteration_list, location_list, source, output_path, source_velocities, date_list, velocitie_dict, AreaList, GaussDatePointsList, gaussLocationList, gaussIterationList, gauss_ampList, velocitie_dict_componet):
            __slots__ = ['grid', 'polarization', 'labels', 'lines', 'iteration_list', 'location_list', 'source', 'output_path', 'new_spectre', 'spectrumSet', 'plotList', 'months', 'dateList', 'velocitie_dict', 'periodPlotSet', 'AreaList', 'GaussDatePointsList', 'gaussLocationList', "gaussIterationList", "gauss_ampList", "velocitie_dict_componet", "bad_points", "selected_points", "multiple_spectre"]
            super().__init__()
            self.setWindowTitle("Monitoring")
            self._addWidget(self.createControlGroup(), 1, 1)
             
            self.polarization = "polarization AVG"
            self.labels = list()
            self.lines = list()
            self.iteration_list = iteration_list
            self.location_list = location_list
            self.source = source
            self.output_path = output_path
            self.new_spectre = True
            self.spectrumSet = set()
            self.plotList = list()
            self.months = Months()
            self.source_velocities = source_velocities
            self.dateList = date_list
            self.velocitie_dict = velocitie_dict
            self.periodPlotSet = set()
            self.AreaList = AreaList
            self.GaussDatePointsList = GaussDatePointsList
            self.gaussLocationList = gaussLocationList
            self.gaussIterationList = gaussIterationList
            self.gauss_ampList = gauss_ampList
            self.velocitie_dict_componet = velocitie_dict_componet
            self.bad_points = list()
            self.selected_points = list()
            self.multiple_spectre = False
            
        def setLineDict(self, lineDict):
            self.lineDict = lineDict
              
        def createPeriodView(self):
            symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
            colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
            velocity = self.componentInput.text()
            if velocity != "":
                velocityIndex = self.source_velocities.index(velocity)
                plotSimbol = symbols[velocityIndex] + colors[velocityIndex]
                self.period_view = Period_View()
                self.period_view.PlotPeriods(self.dateList, self.velocitie_dict, self.source_velocities, velocityIndex, plotSimbol)
                self.periodPlotSet.add(self.period_view)
            else:
                print("No velocity choosed")

        def createGaussAreaView(self):
            self.gauss_view = Gauss_View(self.AreaList, self.source, self.GaussDatePointsList, self.gaussLocationList, self.gaussIterationList)
            self.gauss_view.plot()
            self.gauss_view.show()

        def createGaussAmplitudesView(self):
            self.gauss_view2 = Gauss_View2(self.gauss_ampList, self.source, self.GaussDatePointsList, self.gaussLocationList, self.gaussIterationList, self.velocitie_dict_componet, self.dateList)
            self.gauss_view2.plot()
            self.gauss_view2.show()

        def createMapVew(self):
            os.system("perl " + "code/find_multiple.pl " + self.output_path + " " + self.source)
            self.maps_view = Maps_View(self.source)
            self.maps_view.plotMaps(self.output_path)
            self.maps_view.show()
        
        def createControlGroup(self):
            groupBox = QGroupBox("")
            groupBox.setFixedWidth(120)
            groupBox.setFixedHeight(120)
            
            comboBox = QComboBox(self)
            comboBox.addItem("polarization AVG")
            comboBox.addItem("polarization U1")
            comboBox.addItem("polarization U9")
            comboBox.addItem("ALL")
            comboBox.activated[str].connect(self.getPolarization)    
                  
            controlGrid = QGridLayout()
            controlGrid.addWidget(comboBox, 5, 0)
            plotPeriodsbutton = QPushButton('Plot periods', self)
            plotPeriodsbutton.clicked.connect(self.createPeriodView)
            controlGrid.addWidget(plotPeriodsbutton, 1, 0)
            plotMapsbutton = QPushButton('Plot maps', self)
            plotMapsbutton.clicked.connect(self.createMapVew)

            plotGaussAreasbutton = QPushButton('Plot Gauss areas', self)
            plotGaussAreasbutton.clicked.connect(self.createGaussAreaView)
            controlGrid.addWidget(plotGaussAreasbutton, 3, 0)

            plotGaussAmplitudesbutton = QPushButton('Plot Gauss amplitudes', self)
            plotGaussAmplitudesbutton.clicked.connect(self.createGaussAmplitudesView)
            controlGrid.addWidget(plotGaussAmplitudesbutton, 4, 0)

            controlGrid.addWidget(plotMapsbutton, 2, 0)
            self.componentInput = QLineEdit()
            self.componentInput.setFixedWidth(100)
            controlGrid.addWidget(self.componentInput, 0, 0)
            groupBox.setLayout(controlGrid)
            return groupBox
        
        def setPolarization(self, polarization):
            self.polarization = polarization
            
        def setLabels(self, labels):
            self.labels = labels
             
        def setLines(self, lines):
            self.lines = lines
        
        def getIndexiesOfPolarization(self):
            polarization = self.polarization.split(" ")[-1]
            
            if polarization == "ALL":
                for line in self.lines:
                    line.set_picker(5)
                    line.set_visible(True)            
                 
            else:
                if polarization == "AVG":
                    for line in self.lineDict["avg"]:
                        line.set_picker(5)
                        line.set_visible(True)
                        
                    for line in self.lineDict["u1"]:
                        line.set_picker(False)
                        line.set_visible(False)
                        
                    for line in self.lineDict["u9"]:
                        line.set_picker(False)
                        line.set_visible(False)
                        
                elif polarization == "U1":
                    for line in self.lineDict["u1"]:
                        line.set_picker(5)
                        line.set_visible(True)
                        
                    for line in self.lineDict["avg"]:
                        line.set_picker(False)
                        line.set_visible(False)
                        
                    for line in self.lineDict["u9"]:
                        line.set_picker(False)
                        line.set_visible(False)
                
                elif polarization == "U9":
                    for line in self.lineDict["u9"]:
                        line.set_picker(5)
                        line.set_visible(True)
                        
                    for line in self.lineDict["u1"]:
                        line.set_picker(False)
                        line.set_visible(False)
                        
                    for line in self.lineDict["avg"]:
                        line.set_picker(False)
                        line.set_visible(False)
                                            
        def getPolarization(self, polarization):
            self.setPolarization(polarization)
            self.getIndexiesOfPolarization()
            
        def keyPressEvent(self, e):
            if e.key() == Qt.Key_Shift:
                self.new_spectre = True

            elif e.key() == Qt.Key_Alt:
                if self.multiple_spectre:
                    self.multiple_spectre = False
                else:
                    self.multiple_spectre = True


        def setMonitoringPlot(self, plot):
            self.monitoringPlot = plot
                
        def chooseSpectrum(self, event):

            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            index = [ind][0]
            selected_point = (xdata[ind], ydata[ind])
            self.selected_points.append(selected_point)

            iteration = MonitoringViewHelper.getIteration(self.iteration_list, int(index[0]))
            resultFileName = self.source + ".json"

            with open(getConfigs("paths", "resultFilePath") + resultFileName) as result_data:
                results = json.load(result_data)

            p = tuple(zip(xdata[ind], ydata[ind]))

            if event.mouseevent.button == 1:
                if self.multiple_spectre:

                    date1 = MonitoringViewHelper.formatDate(xdata, [int(index) - 1])
                    location1 = MonitoringViewHelper.getLocation(self.location_list, int(index[0]) - 1)
                    iteration1 = str(int(MonitoringViewHelper.getIteration(self.iteration_list, int(index[0]) -1)))

                    date2 = MonitoringViewHelper.formatDate(xdata, [int(index)])
                    location2 = MonitoringViewHelper.getLocation(self.location_list, int(index[0]))
                    iteration2 = str(int(MonitoringViewHelper.getIteration(self.iteration_list, int(index[0]))))

                    date3 = MonitoringViewHelper.formatDate(xdata, [int(index) + 1])
                    location3 = MonitoringViewHelper.getLocation(self.location_list, int(index[0]) + 1)
                    iteration3 = str(int(MonitoringViewHelper.getIteration(self.iteration_list, int(index[0]) + 1)))

                    source_path = self.source + "/" + self.source + "_"

                    spectra_file_name1 = source_path + date1 + "_" + location1 + "_" + iteration1 + ".dat"
                    spectra_file_name2 = source_path + date2 + "_" + location2 + "_" + iteration2 + ".dat"
                    spectra_file_name3 = source_path + date3 + "_" + location3 + "_" + iteration3 + ".dat"
                    self.plotSpecter([spectra_file_name1, spectra_file_name2, spectra_file_name3], self.polarization)

                else:
                    spectraFileName = self.source + "/" + self.source + "_" + MonitoringViewHelper.formatDate(xdata, index) + "_" + MonitoringViewHelper.getLocation(self.location_list, int(index[0])) + "_" + MonitoringViewHelper.getIteration(self.iteration_list,int(index[0])) + ".dat"
                    self.plotSpecter(spectraFileName, self.polarization)

            elif event.mouseevent.button == 2:
                good_index = self.selected_points.index(selected_point)
                self.bad_points[good_index].set_marker(None)
                self.monitoringPlot.canvasShow()

                for experiment in results:
                    if experiment.endswith("_" + iteration):
                        results[experiment]["flag"] = False

                with open(getConfigs("paths", "resultFilePath") + resultFileName, "w") as result_data:
                    result_data.write(json.dumps(results, indent=2))

            elif event.mouseevent.button == 3:
                line = self.monitoringPlot.plot(p[0][0], p[0][1], "rx", markersize=10)
                self.bad_points.append(line[0])
                self.monitoringPlot.canvasShow()

                for experiment in results:
                    if experiment.endswith("_" + iteration):
                        results[experiment]["flag"] = True

                with open(getConfigs("paths", "resultFilePath") + resultFileName, "w") as result_data:
                    result_data.write(json.dumps(results, indent=2))

        def plotSpecter(self, spectraFileName, polarization):

            amplitude_colon = 3
            if polarization == "U1":
                amplitude_colon = 1
            elif polarization == "U9":
                amplitude_colon = 2
            elif polarization == "AVG" or polarization=="ALL":
                amplitude_colon = 3

            if self.new_spectre:
                self.Spectre_View = Spectre_View()
                self.spectrumSet.add(self.Spectre_View)
                self.spectrPlot = Plot()
                if isinstance(spectraFileName, list):
                    source_name = getConfigs("Full_source_name", spectraFileName[0].split("/")[-1].split("_")[0])
                    self.spectrPlot.creatPlot(self.Spectre_View.getGrid(), "Velocity (km sec$^{-1}$)", "Flux density (Jy)", source_name, (1, 0), "linear")
                else:
                    source_name = getConfigs("Full_source_name", spectraFileName.split("/")[-1].split("_")[0])
                    self.spectrPlot.creatPlot(self.Spectre_View.getGrid(), "Velocity (km sec$^{-1}$)","Flux density (Jy)", source_name, (1, 0), "linear")
                self.Spectre_View._addWidget(self.spectrPlot, 0, 0)
                self.Spectre_View.show()
                self.new_spectre = False
                self.plotList.clear()

            if self.multiple_spectre:
                for sf in spectraFileName:
                    file_name = self.output_path + sf
                    data = np.fromfile(file_name, dtype="float64", count=-1, sep=" ").reshape((file_len(file_name), 4))
                    tmpDate = sf.split("/")[-1].split("_")
                    tmpDate[-4] = self.months.getMonthNumber([tmpDate[-4]][0])
                    plot_name = datetime.datetime.strptime(" ".join(tmpDate[1:-2]), "%H %M %S %d %m %Y")

                    x = data[:, [0]]
                    y = data[:, [amplitude_colon]]

                    if plot_name not in self.plotList:
                        self.plotList.append(plot_name)
                        self.spectrPlot.plot(x, y, "-", label=plot_name)
                        self.spectrPlot.canvasShow()

            else:
                spectraFileName = self.output_path + spectraFileName
                data = np.fromfile(spectraFileName, dtype="float64", count=-1, sep=" ") .reshape((file_len(spectraFileName),4))
                tmpDate = spectraFileName.split("/")[-1].split("_")
                tmpDate[-4] = self.months.getMonthNumber([tmpDate[-4]][0])
                plot_name = datetime.datetime.strptime(" ".join(tmpDate[1:-2]), "%H %M %S %d %m %Y")

                x = data[:, [0]]
                y = data[:, [amplitude_colon]]

                if plot_name not in self.plotList:
                    self.plotList.append(plot_name)
                    self.spectrPlot.plot(x, y, "-", label=plot_name)
                    self.spectrPlot.canvasShow()


class MonitoringApp(QWidget):
    __slots__ = ['configFilePath', 'grid', 'months', 'new_spectre', 'flag']

    def __init__(self, configFilePath):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        self.configFilePath = configFilePath
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)
        self.chooseSource()
        self.new_spectre = True
        self.months = Months()
        self.flag = "All"
        
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    def __addWidget(self, widget, row, colomn):
        self.grid.addWidget(widget, row, colomn)
        
    def chooseSource(self):
        self.setWindowTitle("Choose Source")
        chooseLabel = QLabel("Choose source")
        self.__addWidget(chooseLabel, 0, 0)
        self.sourceInput = QLineEdit()
        self.__addWidget(self.sourceInput, 1, 0)
        chooseSourceButton = QPushButton("Ok", self)
        chooseSourceButton.clicked.connect(self.plot)
        chooseSourceButton.setStyleSheet("background-color: green")
        self.__addWidget(chooseSourceButton, 1, 2)

        comboBox2 = QComboBox(self)
        comboBox2.addItem("All")
        comboBox2.addItem("Not Flag")
        comboBox2.activated[str].connect(self.setFlags)
        self.__addWidget(comboBox2, 1, 1)

    def setFlags(self, flag):
        self.flag = flag
        
    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Return:
            if len (self.sourceInput.text()) > 1:
                self.plot()
                                                                                                                                                        
    def plotMonitoring(self, resultDir, source_velocities, source):
        result_list = list()
        velocitie_dict = {"u1": dict(), "u9": dict(), "avg": dict()}
        iteration_list = list()
        date_list = list()
        
        for velocitie in velocitie_dict:
            for vel in source_velocities:
                velocitie_dict[velocitie][vel] = list()
                                
        resultFileName = source + ".json"
        
        with open(resultDir + resultFileName) as result_data:    
                results = json.load(result_data)
        
        self.location_list = list()
        labels2 = list()
        badIteration = list()
        for experiment in results:
                scanData = results[experiment]
                date = scanData["Date"]
                location = scanData["location"]
                amplitudes_for_u1 = scanData["polarizationU1"] # Got poitns for all experiments for polarization u1
                amplitudes_for_u9 = scanData["polarizationU9"] # Got poitns for all experiments for polarization u9
                amplitudes_for_uAVG = scanData["polarizationAVG"] # Got poitns for all experiments for polarization uAVG
                iter_number = scanData["Iteration_number"]
                specie = scanData["specie"]
                type = scanData["type"]
                modifiedJulianDays = scanData["modifiedJulianDays"]
                flag = scanData["flag"]

                if "areas" in scanData.keys():
                    areas = scanData["areas"]
                else:
                    print("No Gauss areas in result file for iteration", iter_number)

                if "gauss_amp" in scanData.keys():
                    gauss_amp = scanData["gauss_amp"]
                else:
                    badIteration.append(iter_number)
                    gauss_amp = np.zeros(len(getConfigs("gauss_lines", self.source).replace(" ", "").split(",")))
                    print("No Gauss gauss_amp in result file for iteration", iter_number)

                dates = date.split("_")
                monthsNumber = dates[1]
                dates[1] = self.months.getMonthNumber([monthsNumber][0])
                date = scanData["time"].replace(":", " ") + " " + " ".join(dates)

                if self.flag == "Not Flag":
                    if flag == False:
                        result = Result(location, datetime.datetime.strptime(date, '%H %M %S %d %m %Y'), amplitudes_for_u1, amplitudes_for_u9, amplitudes_for_uAVG, iter_number, specie, type, modifiedJulianDays, areas, gauss_amp)
                        result_list.append(dict(result))

                elif self.flag == "All":
                    result = Result(location, datetime.datetime.strptime(date, '%H %M %S %d %m %Y'), amplitudes_for_u1,amplitudes_for_u9, amplitudes_for_uAVG, iter_number, specie, type, modifiedJulianDays, areas, gauss_amp)
                    result_list.append(dict(result))

        print("badIteration", badIteration)
        result_list = sorted(result_list, key=itemgetter('date'), reverse=False)

        self.AreaList = list()
        self.GaussDatePointsList = list()
        self.gaussLocationList = list()
        self.gaussIterationList = list()
        self.gauss_ampList = list()
        modifiedJulianDaysList = list()

        for experiment in result_list:
            u1 = experiment["polarizationU1"]
            u9 = experiment["polarizationU9"]
            avg = experiment["polarizationUAVG"]
            iteration_list.append(experiment["iteration_number"])
            date_list.append(experiment["date"])
            location = experiment["location"]
            specie = experiment["specie"]
            areas = experiment["areas"]
            gauss_amp = experiment["gauss_amp"]
            gauss_amp = [float(g) for g in gauss_amp]

            gaussLines = getConfigs("gauss_lines", self.source).replace(" ", "").split(",")
            if len(areas) == len(gaussLines):
                self.AreaList.append(areas)
                self.GaussDatePointsList.append(experiment["date"])
                self.gaussLocationList.append(location)
                self.gaussIterationList.append(experiment["iteration_number"])
                self.gauss_ampList.append(gauss_amp)
            self.location_list.append(location)
            
            for i in u1:
                for vel in source_velocities:
                    if float(vel) == float(i[0]):
                        velocitie_dict["u1"][vel].append(i[1]) 
            
            for j in u9:
                for vel in source_velocities:
                    if float(vel) == float(j[0]):
                        velocitie_dict["u9"][vel].append(j[1]) 
                        
            for k in avg:
                for vel in source_velocities:
                    if float(vel) == float(k[0]):
                        velocitie_dict["avg"][vel].append(k[1]) 

            type = experiment["type"]
            label = "Station is " + location + "\n" + "Date is " + experiment["date"].strftime('%d %m %Y') + "\n " + "Iteration number " + str(experiment["iteration_number"]) + "\n " + "Specie " + str(specie) + "\n " + "Type " + type
            labels2.append(label)
            modifiedJulianDaysList.append(experiment["modifiedJulianDays"])
            
        self.iteration_list = iteration_list
        Symbols = ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        
        self.Monitoring_View = Monitoring_View(self.iteration_list, self.location_list, self.source, self.output_path, source_velocities, date_list, velocitie_dict, self.AreaList, self.GaussDatePointsList, self.gaussLocationList, self.gaussIterationList, self.gauss_ampList, velocitie_dict["avg"][source_velocities[0]])
        self.monitoringPlot = Plot()
        self.monitoringPlot.creatPlot(self.Monitoring_View.getGrid(), "Time", "Flux density (Jy)", getConfigs("Full_source_name", self.source), (1,0), "log")
        self.Monitoring_View.setMonitoringPlot(self.monitoringPlot)
        
        def convertDatetimeObjectToMJD(time):
            time=time.isoformat()
            t=Time(time, format='isot')
            return t.mjd
        
        #.strftime("%H %M %d %m %Y")
        #x = [date.strftime("%H %M %d %m %Y") for date in date_list]

        # MODIFICETS
        
        lineDict = {"u1":list(), "u9":list(), "avg":list()}
        lines = list()
        monitoringResults = [[convertDatetimeObjectToMJD(d) for d in date_list]]

        for i in range(0, len(source_velocities)):
            l1, = self.monitoringPlot.plot(date_list, velocitie_dict["u1"][source_velocities[i]], Symbols[i]+colors[i], fontsize=8, visible=False, picker=False)
            l2, = self.monitoringPlot.plot(date_list, velocitie_dict["u9"][source_velocities[i]], Symbols[i]+colors[i], fontsize=8, visible=False, picker=False)
            l3, = self.monitoringPlot.plot(date_list, velocitie_dict["avg"][source_velocities[i]], Symbols[i]+colors[i], fontsize=8, label="Velocity " + source_velocities[i], visible=True, picker=5)
            monitoringResults.append(velocitie_dict["avg"][source_velocities[i]])
            lines.append(l1)
            lines.append(l2)
            lines.append(l3)
            lineDict["u1"].append(l1)
            lineDict["u9"].append(l2)
            lineDict["avg"].append(l3)

        np.savetxt("monitoring/" + self.source + ".txt", np.transpose(monitoringResults))
        self.Monitoring_View._addWidget(self.monitoringPlot, 0, 0)
        #self.monitoringPlot.setXtics(date_list, [convertDatetimeObjectToMJD(date) for date in  date_list], '30')
        
        self.monitoringPlot.addCursor(labels2)
        labels = [str(line.get_label()) for line in lines]
        #self.monitoringPlot.addSecondAxis2(modifiedJulianDaysList, label="Modified Julian Days", axiss=None)
        self.Monitoring_View.setLineDict(lineDict)
        self.Monitoring_View.setLabels(labels)
        self.Monitoring_View.setLines(lines)
        self.monitoringPlot.addPickEvent(self.Monitoring_View.chooseSpectrum)
        self.Monitoring_View.showMaximized()
        self.Monitoring_View.show()
        
    def plot(self):
        config = ConfigParser.getInstance()
        config.CreateConfig(self.configFilePath)
        resultDir = config.getConfig('paths', "resultFilePath")
        self.output_path = config.getConfig('paths', "outputFilePath")
        self.source = self.sourceInput.text()
        source_velocities = config.getConfig('velocities', self.source).split(",")
        source_velocities = [x.strip() for x in source_velocities]
        self.plotMonitoring(resultDir, source_velocities, self.source)


class Main(object):
    
    def __parse_arguments(self):
        self.args = parseArguments()
        self.configFilePath = str(self.args.__dict__["config"])
        
    def __create_app(self):
        qApp = QApplication(sys.argv)
        aw = MonitoringApp(self.configFilePath)
        aw.show()
        sys.exit(qApp.exec_())
        sys.exit(0)
    
    @classmethod
    def run(cls):
        Main.__parse_arguments(cls)
        Main.__create_app(cls)


def main():
    Main().run()


if __name__ == "__main__":
    p = Process(target=main,)
    p.start()
    p.join()
    sys.exit(0)
