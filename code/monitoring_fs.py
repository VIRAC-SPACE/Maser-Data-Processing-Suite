#! /usr/bin/python3
# -*- coding: utf-8 -*-
import sys
import matplotlib.pyplot  as plt
from matplotlib.widgets import CheckButtons
import mplcursors
import datetime
import json
import argparse
import configparser
from operator import itemgetter
import numpy as np

from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QPushButton, QLabel, QLineEdit, QDesktopWidget)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

#from ploting_qt5 import  Plot
from result import  Result

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Monitoring velocity amplitudes in time. ''',
    epilog="""Monitor.""")

    # Positional mandatory arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 2.0')

    # Parse arguments
    args = parser.parse_args()

    return args

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class Monitoring(QWidget):
    def __init__(self, configFilePath):
        super().__init__()
        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        
        self.configFilePath = configFilePath
        
        self.grid = QGridLayout()
        self.setLayout(self.grid)
        self.grid.setSpacing(10)
        
        self.chooseSource()
        
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
    def chooseSource(self):
        self.setWindowTitle("Choose Source")
        
        self.chooseLabel = QLabel("Choose source")
        self.grid.addWidget(self.chooseLabel, 0, 0)
        
        self.sourceInput = QLineEdit()
        self.grid.addWidget(self.sourceInput, 1, 0)
        
        self.chooseSource = QPushButton("Ok", self)
        self.chooseSource.clicked.connect(self.plot)
        self.chooseSource.setStyleSheet("background-color: green")
        self.grid.addWidget(self.chooseSource, 1, 1)
        
    def keyPressEvent(self, e):
        
        if e.key() == Qt.Key_Return:
            if len (self.sourceInput.text()) > 1:
                self.plot()
            
    def chooseSpectrum(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ind = event.ind
        index = [ind][0]
        polarization = thisline.get_label().split()[1]
        date = xdata[index][0].strftime("%H %M %S %d %m %Y").split()
        month = datetime.date(1900, int(date[-2]) , 1).strftime('%B')[0:3].title().replace("ū", "u").replace("i", "y").replace("k", "c")
        date[-2] = month
        date = "_".join(date)
        iteration = self.iteration_list[int(index)]
        location = self.location_list[int(index)]
        spectraFileName = self.source + "_" + date + "_" + location + "_"  + str(iteration) + ".dat"
        self.plotSpecter(spectraFileName, polarization)
        
    def plotSpecter(self, spectraFileName, polarization):
        if polarization == "U1":
            amplitude_colon = 1
        elif polarization == "U9":
            amplitude_colon = 2
        elif polarization == "AVG":
            amplitude_colon = 3
            
        spectraFileName = self.output_path + spectraFileName
        
        data = np.fromfile(spectraFileName, dtype="float64", count=-1, sep=" ") .reshape((file_len(spectraFileName),4))
        plot_name =  " ".join(spectraFileName.split("/")[-1].split("_")[0:-2])
        
        x = data[:, [0]]
        y = data[:, [amplitude_colon]]
        
        fig = plt.figure("spectrum  for " + plot_name + " Polarization " + polarization)
        graph = fig.add_subplot(111)
        
        graph.plot(x,y)
        
        plt.show()

    def plotMonitoring(self, resultDir, source_velocities, source):
        result_list = list()
        velocitie_dict = {"u1":dict(), "u9":dict(), "avg":dict()}
        iteration_list = list()
        date_list = list()
        
        for velocitie in velocitie_dict:
            for vel in source_velocities:
                velocitie_dict[velocitie][vel] = list()
                                
        resultFileName = source + ".json"
        
        months = {"Jan":"1", "Feb":"2", "Mar":"3", "Apr":"4", "May":"5", "Jun":"6", "Jul":"7", "Aug":"8", "Sep":"9", "Oct":"10", "Nov":"11", "Dec":"12"}
        
        with open(resultDir + resultFileName) as result_data:    
                results = json.load(result_data)
        
        self.location_list = list()
        labels2 = list()
        for experiment in results:
                scanData = results[experiment]
                date = scanData["Date"]
                location = scanData["location"]
                amplitudes_for_u1 = scanData["polarizationU1"] # Got poitns for all experiments for polarization u1
                amplitudes_for_u9 = scanData["polarizationU9"] # Got poitns for all experiments for polarization u9
                amplitudes_for_uAVG = scanData["polarizationAVG"] # Got poitns for all experiments for polarization uAVG
                iter_number = scanData["Iteration_number"]
                dates = date.split("_")
                monthsNumber = dates[1]
                dates[1] = months[monthsNumber]
                date = scanData["time"].replace(":", " ") + " " +  " ".join(dates) 
            
                result = Result(location, datetime.datetime.strptime(date, '%H %M %S %d %m %Y'), amplitudes_for_u1, amplitudes_for_u9, amplitudes_for_uAVG, iter_number)
                
                result_list.append(dict(result))
              
        result_list = sorted(result_list, key=itemgetter('date'), reverse=False)
         
        for experiment in result_list:
            u1 = experiment["polarizationU1"]
            u9 = experiment["polarizationU9"]
            avg = experiment["polarizationUAVG"]
            iteration_list.append(experiment["iteration_number"])
            date_list.append(experiment["date"])
            location = experiment["location"]
            
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
                         
            label = "Station is " + location + "\n" + "Date is " + experiment["date"].strftime('%d %m %Y') + "\n " + "iteration number " + str(experiment["iteration_number"])
            labels2.append(label)
            
        self.iteration_list = iteration_list
       
        x = list()
        for a in range(0, len(date_list)):
            x.append(date_list[a])
        
        Symbols =  ["*", "o", "v", "^", "<", ">", "1", "2", "3", "4"]
        
        lines = list()
        
        fig = plt.figure("Monitoring for " + source)
        graph = fig.add_subplot(111)
        for i in range(0, len(source_velocities)):
            l1, = graph.plot(x, velocitie_dict["u1"][source_velocities[i]], Symbols[i]+"r", label="polarization U1 " + "Velocity " + source_velocities[i], visible=False, picker=False)
            l2, = graph.plot(x, velocitie_dict["u9"][source_velocities[i]], Symbols[i]+"g", label="polarization U9 " + "Velocity " + source_velocities[i], visible=False, picker=False)
            l3, = graph.plot(x, velocitie_dict["avg"][source_velocities[i]], Symbols[i]+"b", label="polarization AVG " + "Velocity " + source_velocities[i], visible=True, picker=5)
            
            lines.append(l1)
            lines.append(l2)
            lines.append(l3)
        
        plt.xticks(x, [date.strftime("%H %M %S %d %m %Y") for date in  date_list], rotation='30')
        
        plt.legend()
        
        cursor =  mplcursors.cursor(hover=True, highlight=True)
        cursor.connect("add", lambda sel: sel.annotation.set_text(labels2[sel.target.index]))
        
        labels = list()
        labels = [str(line.get_label()) for line in lines]
        visibility = [line.get_visible() for line in lines]
        
        check = CheckButtons(plt.axes([0.8, 0.1, 0.1, 0.7]),  labels, visibility)
        
        def func(label):
            index = labels.index(label)
            lines[index].set_visible(not lines[index].get_visible())
            lines[index].set_picker = 5
            plt.draw()
            
        check.on_clicked(func)
        
        cid = fig.canvas.mpl_connect('pick_event', self.chooseSpectrum)
       
        plt.show()
        
    def plot(self):
        #Creating config parametrs
        config = configparser.RawConfigParser()
        config.read(self.configFilePath)
        resultDir = config.get('paths', "resultFilePath")
        self.output_path = config.get('paths', "output")
        self.source = self.sourceInput.text()
        source_velocities = config.get('velocities', self.source).split(",")
            
        self.plotMonitoring(resultDir, source_velocities, self.source)
        
def main():
    # Parse the arguments
    args = parseArguments()
    configFilePath = str(args.__dict__["config"])
     
    #Create App
    qApp = QApplication(sys.argv)
    aw = Monitoring(configFilePath)
    aw.show()
    sys.exit(qApp.exec_())
    
    sys.exit(0)
    
    sys.exit(0)
    
if __name__=="__main__":
    main()
    