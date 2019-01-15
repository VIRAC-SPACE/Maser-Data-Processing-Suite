#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import _thread
import argparse
import json
import numpy as np
import matplotlib
import datetime
from operator import itemgetter

from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QPushButton, QLabel, QLineEdit, QDesktopWidget, QComboBox, QGroupBox)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

from ploting_qt5 import  Plot
from parsers._configparser import ConfigParser
from result import  Result
from monitor.months import Months
from monitor.monitoringViewHelper import MonitoringViewHelper

matplotlib.use('Qt5Agg')

def parseArguments():
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''', epilog="""Monitor.""")
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')
    args = parser.parse_args()
    return args

class PlotingView(QWidget):
    def __init__(self):
        super().__init__()
        self.grid = QGridLayout()
        self.grid.setSpacing(10)
        self.setLayout(self.grid)
        self.show()
        
    def _addWidget(self, widget, row, collon):
        self.grid.addWidget(widget, row, collon)
            
    def getGrid(self):
        return self.grid

class TimeView(PlotingView):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Monitoring")
        self.showMaximized()
        
        def createPeriodView():
            #self.Period_View = Period_View()
            #self.Period_View.PlotPeriods()
            print("qwerty")
        
        def createControlGroup():
            groupBox = QGroupBox("")
                
            comboBox = QComboBox(self)
            comboBox.addItem("polarization AVG")
            comboBox.addItem("polarization U1")
            comboBox.addItem("polarization U9")
            comboBox.addItem("ALL")
            comboBox.activated[str].connect(self.getPolarization)
                
            controlGrid = QGridLayout()
            controlGrid.addWidget(comboBox, 2, 0)
            plotPeriodsbutton = QPushButton('Plot periods', self)
            plotPeriodsbutton.clicked.connect(createPeriodView)
            controlGrid.addWidget(plotPeriodsbutton, 1, 0)
            componentInput = QLineEdit()
            controlGrid.addWidget(componentInput, 0, 0)
            groupBox.setLayout(controlGrid)
                
            return groupBox
        
        self._addWidget(createControlGroup(), 1, 2)
        
    def getPolarization(self, polarization):
            self.setPolarization(polarization)
            self.getIndexiesOfPolarization(self.labels)
        
    def createPlot(self, x, source_velocities, velocitie_dict, cursorLabels):
        lines = list()
        Symbols =  ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        monitoringPlot = Plot()
        monitoringPlot.creatPlot(self.getGrid(), "Time", "Flux density (Jy)", None, (1,0))
       
        for i in range(0, len(source_velocities)):
            l1, = monitoringPlot.plot(x, velocitie_dict["u1"][source_velocities[i]], Symbols[i]+"r", label="polarization U1 " + "Velocity " + source_velocities[i], visible=False, picker=False)
            l2, = monitoringPlot.plot(x, velocitie_dict["u9"][source_velocities[i]], Symbols[i]+"g", label="polarization U9 " + "Velocity " + source_velocities[i], visible=False, picker=False)
            l3, = monitoringPlot.plot(x, velocitie_dict["avg"][source_velocities[i]], Symbols[i]+"b", label="polarization AVG " + "Velocity " + source_velocities[i], visible=True, picker=5)
            
            lines.append(l1)
            lines.append(l2)
            lines.append(l3)
            
        monitoringPlot.setXtics(x, [date.strftime("%H %M %d %m %Y") for date in  x], '30')
        monitoringPlot.addCursor(cursorLabels)
        self._addWidget(monitoringPlot, 0, 0)
                    
class SpectralView(PlotingView):
    def __init__(self):
        super().__init__()

class MonitoringApp(QWidget):
    def __init__(self, configFilePath):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)
        self.configFilePath = configFilePath
        self.chooseSource()
        self.source = ""
        self.months = Months()
        
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
        chooseSourceButton.clicked.connect(self.plotTimes)
        chooseSourceButton.setStyleSheet("background-color: green")
        self.__addWidget(chooseSourceButton, 1, 1)
        
    def keyPressEvent(self, e):
        if e.key() == Qt.Key_Return:
            if len (self.sourceInput.text()) > 1:
                self.plotTimes()
                       
    def plotTimes(self):
        self.source = self.sourceInput.text()
        self.timeView = TimeView()
        self.timeView.createPlot(self.__createDateList(), self.getSourceVelocities(), self.__createVelocitie_dict(), self.__getCursorLabel())
        
    def getSource(self):
        return self.source
        
    def getConfig(self, key, value):
        config = ConfigParser.getInstance()
        config.CreateConfig(self.configFilePath)
        return config.getConfig(key, value)
            
    def getSourceVelocities(self):
        source_velocities = self.getConfig('velocities', self.getSource()).replace(" ", "").split(",")
        return source_velocities
    
    def getResultDir(self):
        resultDir = self.getConfig('paths', 'resultFilePath')
        return resultDir
    
    def getResults(self):
        resultFileName = self.getSource() + ".json"
        
        with open(self.getResultDir() + resultFileName) as result_data:
            results = json.load(result_data)
            
        return results
    
    def __createResultList(self):
        result_list = list()
        
        for experiment in self.getResults():
                scanData = self.getResults()[experiment]
                date = scanData["Date"]
                location = scanData["location"]
                amplitudes_for_u1 = scanData["polarizationU1"] # Got poitns for all experiments for polarization u1
                amplitudes_for_u9 = scanData["polarizationU9"] # Got poitns for all experiments for polarization u9
                amplitudes_for_uAVG = scanData["polarizationAVG"] # Got poitns for all experiments for polarization uAVG
                iter_number = scanData["Iteration_number"]
                specie = scanData["specie"]
                dates = date.split("_")
                monthsNumber = dates[1]
                dates[1] = self.months.getMonthNumber([monthsNumber][0])
                date = scanData["time"].replace(":", " ") + " " +  " ".join(dates) 
                result = Result(location, datetime.datetime.strptime(date, '%H %M %S %d %m %Y'), amplitudes_for_u1, amplitudes_for_u9, amplitudes_for_uAVG, iter_number, specie)
                result_list.append(dict(result))
              
        result_list = sorted(result_list, key=itemgetter('date'), reverse=False)
        return result_list
    
    def __createDateList(self):
        date_list = list()
        for experiment in self.__createResultList():
            date_list.append(experiment["date"])
            
        return date_list
    
    def __getCursorLabel(self):
        labels = list()
        for experiment in self.__createResultList():
            location = experiment["location"]
            specie = experiment["specie"]
            label = "Station is " + location + "\n" + "Date is " + experiment["date"].strftime('%d %m %Y') + "\n " + "Iteration number " + str(experiment["iteration_number"]) + "\n " + "Specie " + str(specie)
            labels.append(label)
        return labels
            
    def __createVelocitie_dict(self):
        velocitie_dict = {"u1":dict(), "u9":dict(), "avg":dict()}
        
        for velocitie in velocitie_dict:
            for vel in self.getSourceVelocities():
                velocitie_dict[velocitie][vel] = list()
                   
        for experiment in self.__createResultList():
            u1 = experiment["polarizationU1"]
            u9 = experiment["polarizationU9"]
            avg = experiment["polarizationUAVG"]
            
            for i in u1:
                for vel in self.getSourceVelocities():
                    if float(vel) == float(i[0]):
                        velocitie_dict["u1"][vel].append(i[1]) 
            
            for j in u9:
                for vel in self.getSourceVelocities():
                    if float(vel) == float(j[0]):
                        velocitie_dict["u9"][vel].append(j[1]) 
                        
            for k in avg:
                for vel in self.getSourceVelocities():
                    if float(vel) == float(k[0]):
                        velocitie_dict["avg"][vel].append(k[1])
                        
        return velocitie_dict

class Main():
    
    def __parseArguments(self):
        self.args = parseArguments()
        self.configFilePath = str(self.args.__dict__["config"])
        
    def __CreateApp(self): 
        qApp = QApplication(sys.argv)
        aw = MonitoringApp(self.configFilePath)
        aw.show()
        sys.exit(qApp.exec_())
        sys.exit(0)
    
    def run(self):
        self.__parseArguments()
        self.__CreateApp()
        
def main():
    app = Main()
    app.run()
    
if __name__=="__main__":
    main() 
