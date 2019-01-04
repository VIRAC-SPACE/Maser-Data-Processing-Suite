#! /usr/bin/python3
# -*- coding: utf-8 -*-

import sys
import os
import argparse
import configparser
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.time import Time
from datetime import datetime
import peakutils
import json
import pickle
from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QPushButton, QMessageBox, QLabel, QLineEdit, QSlider, QDesktopWidget, QLCDNumber)
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon
from PyQt5.QtGui import QColor

from ploting_qt5 import  Plot

def parseArguments():
    parser = argparse.ArgumentParser(description='''plotting tool. ''', epilog="""PRE PLOTTER.""")
    parser.add_argument("datafile", help="Experiment correlation file name", type=str)
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")
    parser.add_argument("-r", "--rawdata", help="Use raw data, skip smoothing", action='store_true')
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 1.0')
    args = parser.parse_args()
    
    return args

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def indexies(array, value):
    indexs = list()
    for i in range(0, len(array)-1):
        if array[i] == value:
            indexs.append(i)
    return indexs

class Result():
    def __init__(self, matrix, specie):
        self.matrix = matrix
        self.specie = specie
        
    def getMatrix(self):
        return self.matrix
    
    def getSpecie(self):
        return self.specie

class Analyzer(QWidget):
    def __init__(self, datafile, resultFilePath, source_velocities, cuts, output, index_range_for_local_maxima, skipsmooth):
        super().__init__()
        self.skipsmoothing = skipsmooth;

        self.setWindowIcon(QIcon('viraclogo.png'))
        self.center()
        
        self.polynomialOrder = 3
        self.source = datafile.split("/")[-1].split(".")[0].split("_")[0]
        self.expername = datafile.split("/")[-1].split(".")[0]
        self.location = datafile.split("/")[-1].split(".")[0].split("_")[-2]
        self.date = "_".join([datafile.split("/")[-1].split(".")[0].split("_")[1], datafile.split("/")[-1].split(".")[0].split("_")[2], datafile.split("/")[-1].split(".")[0].split("_")[3]])
        self.time = datafile.split("/")[-1].split(".")[0].split("_")[-3]
        self.iteration_number = datafile.split("/")[-1].split(".")[0].split("_")[-1]
        self.resultFilePath = resultFilePath
        self.source_velocities = source_velocities
        self.cuts = cuts
        self.output = output
        self.index_range_for_local_maxima = index_range_for_local_maxima
        
        self.infoSet = set()
        self.infoSet_2 = list()
        
        self.changeParms = False
        
        try:
           
            result = pickle.load(open(datafile, "rb"))
            self.specie = result.getSpecie()
            data = result.getMatrix()
        
        except IOError as e:
            print ("IO Error",  e)
            sys.exit(1)
            
        except TypeError as e:
            print ("TypeError",  e)
            sys.exit(1)
            
        except AttributeError as e:
            print ("AttributeError",  e)
            sys.exit(1)
                
        except:
            print("Unexpected error:", sys.exc_info()[0])
            sys.exit(1)
                 
        else:
            self.dataPoints = data.shape[0] 
            self.m = 0
            self.n = self.dataPoints
            
            self.xdata = data[:, [0]]
            self.y_u1 = data[:, [1]] 
            self.y_u9 = data[:, [2]]
            
            #Making sure that data is numpy array
            self.xarray = np.zeros(self.dataPoints)
            self.y1array = np.zeros(self.dataPoints)
            self.y2array = np.zeros(self.dataPoints)
            
            for i in range (0, self.dataPoints):
                self.xarray[i] = self.xdata[i]
                self.y1array[i] = self.y_u1[i]
                self.y2array[i] = self.y_u9[i]
                
            self.xarray =  np.flip(self.xarray,0)
            self.y1array =  np.flip(self.y1array,0)
            self.y2array =  np.flip(self.y2array,0)
            
            self.grid = QGridLayout()
            self.setLayout(self.grid)
            self.grid.setSpacing(10)
            
            self.plotShortSpectr()
            
    def center(self):
        qr = self.frameGeometry()
        cp = QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
    
    def onpickU1(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_3.plot(p[0][0], p[0][1], 'ro',  markersize=1,  picker=5)
        self.points_1u.append(p[0])
        #self.points_9u.append(p[0])
        self.plot_3.canvasShow()
        
    def onpickU9(self, event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_4.plot(p[0][0], p[0][1], 'ro', markersize=1, picker=5)
        self.points_9u.append(p[0])
        #self.points_1u.append(p[0])
        self.plot_4.canvasShow()
                
    def plotInitData(self):
        self.setWindowTitle("Change params")
        
        self.changeParms = True
        
        self.changeParams.hide()
        self.changeParams.close()
        self.grid.removeWidget(self.changeParams)
        del self.changeParams
        
        self.plotPoly.hide()
        self.plotPoly.close()
        self.grid.removeWidget(self.plotPoly)
        del self.plotPoly
        
        self.changeDataButton = QPushButton("Change Data", self)
        self.changeDataButton.clicked.connect(self.changeData)
        self.changeDataButton.setStyleSheet("background-color: blue")
        self.grid.addWidget(self.changeDataButton, 5, 3)
        
        self.plotSmoothDataButton = QPushButton("Create Shorter specter", self)
        self.plotSmoothDataButton.clicked.connect(self.plotShortSpectr)
        self.plotSmoothDataButton.setStyleSheet("background-color: green")
        self.grid.addWidget(self.plotSmoothDataButton, 5, 4)
         
        self.plot_1 = Plot()
        self.plot_1.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "u1 Polarization", (1, 0))
        self.plot_1.plot(self.xarray, self.y1array, 'ko', label='Data Points', markersize=1)
            
        self.plot_2 = Plot()
        self.plot_2.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "u9 Polarization", (1, 1))
        self.plot_2.plot(self.xarray, self.y2array, 'ko', label='Data Points', markersize=1)
            
        self.grid.addWidget(self.plot_1, 0, 0)
        self.grid.addWidget(self.plot_2, 0, 1)
        
        infoPanelLabelsText = ["Polynomial order"]
        infoPanelEntryText = [{"defaultValue":str(self.polynomialOrder), "addEntry":True}]
        
        for i in range(0, len(infoPanelLabelsText)):
            
            self.infoLabel = QLabel(infoPanelLabelsText[i])
            self.grid.addWidget(self.infoLabel, i +  3, 3)
            self.infoSet.add(self.infoLabel)
            
            if  infoPanelEntryText[i]["addEntry"]:
                self.infoInputField = QLineEdit()
                self.infoInputField.setText(str(infoPanelEntryText[i]["defaultValue"]))
                self.grid.addWidget(self.infoInputField, i + 3, 4)
                self.infoSet_2.append(self.infoInputField)
        
        self.a = ((np.abs(self.xarray-float(self.cuts[0][0]))).argmin())
        self.b = ((np.abs(self.xarray-float(self.cuts[-1][1]))).argmin())
        
        print ("self.a, self.b", self.a, self.b)
               
        #sliders
        self.previousM = self.m
        self.previousN = self.n -1
        
        self.m_slider = QSlider(Qt.Horizontal, self)
        self.n_slider = QSlider(Qt.Horizontal, self)
        self.m_slider.setFocusPolicy(Qt.NoFocus)
        self.n_slider.setFocusPolicy(Qt.NoFocus)
        self.m_slider.setTickInterval(20)
        self.m_slider.setSingleStep(20)
        self.m_slider.setMinimum(0) 
        self.m_slider.setMaximum(self.a-1)
        self.n_slider.setMinimum(self.b-1)
        self.n_slider.setMaximum(self.n)
        self.n_slider.setValue(self.n)
        self.m_slider.setMinimumSize(500, 0)
        self.m_slider.setMinimumSize(500, 0)
        self.m_slider.valueChanged[int].connect(self.change_M)
        self.n_slider.valueChanged[int].connect(self.change_N)
        
        self.m_lcd = QLCDNumber(self)
        self.n_lcd = QLCDNumber(self)
        
        self.m_lcd.setSegmentStyle(QLCDNumber.Flat)
        self.n_lcd.setSegmentStyle(QLCDNumber.Flat)
        mpalette = self.m_lcd.palette()
        npalette = self.n_lcd.palette()
        mpalette.setColor(mpalette.Dark, QColor(0, 255, 0))
        npalette.setColor(npalette.Dark, QColor(0, 255, 0))
        self.m_lcd.setPalette(mpalette)
        self.n_lcd.setPalette(npalette)
        
        self.mLabel = QLabel('M', self)
        self.nLabel = QLabel('N', self)
        
        self.grid.addWidget(self.mLabel, 1,3)
        self.grid.addWidget(self.nLabel, 2,3)
        
        self.grid.addWidget(self.m_slider, 1,4)
        self.grid.addWidget(self.n_slider, 2,4)
        
        self.grid.addWidget(self.m_lcd, 1,5)
        self.grid.addWidget(self.n_lcd, 2,5)
        
        self.m_slider.valueChanged.connect(self.m_lcd.display)
        self.n_slider.valueChanged.connect(self.n_lcd.display)
        
        self.m = self.m_slider.value()
        self.n = self.n_slider.value()
            
    def changeData (self):
        newValues = list()
        
        for value in self.infoSet_2:
            newValues.append(value.text())
            
        self.polynomialOrder = float(newValues[0])
        
        QMessageBox.information(self, "Info", "Data was changed")
        
    def change_M(self, value):
        self.plot_1.plot(self.xarray[int(self.previousM)], self.y1array[int(self.previousM)], 'ko', markersize=1)
        self.plot_2.plot(self.xarray[int(self.previousM)], self.y2array[int(self.previousM)], 'ko', markersize=1)
        
        self.plot_1.annotation(self.xarray[int(self.previousM)], self.y1array[int(self.previousM)], " ")
        self.plot_2.annotation(self.xarray[int(self.previousM)], self.y2array[int(self.previousM)], " ")
        
        self.plot_1.remannotation()
        self.plot_2.remannotation()

        self.plot_1.annotation(self.xarray[int(value)], self.y1array[int(value)], "M")
        self.plot_2.annotation(self.xarray[int(value)], self.y2array[int(value)], "M")
        
        self.plot_1.plot(self.xarray[int(value)], self.y1array[int(value)], 'ro', markersize=1)
        self.plot_2.plot(self.xarray[int(value)], self.y2array[int(value)], 'ro', markersize=1)
         
        self.plot_1.canvasShow()
        self.plot_2.canvasShow()
 
        self.previousM = value
    
    def change_N(self, value):
        self.plot_1.plot(self.xarray[int(self.previousN-1)], self.y1array[int(self.previousN-1)], 'ko', markersize=1)
        self.plot_2.plot(self.xarray[int(self.previousN-1)], self.y2array[int(self.previousN-1)], 'ko', markersize=1)
        
        self.plot_1.annotation(self.xarray[int(self.previousN-1)], self.y1array[int(self.previousN-1)], " ")
        self.plot_2.annotation(self.xarray[int(self.previousN-1)], self.y2array[int(self.previousN-1)], " ")

        self.plot_1.remannotation()
        self.plot_2.remannotation()
        
        self.plot_1.annotation(self.xarray[int(value-1)], self.y1array[int(value-1)], "N")
        self.plot_2.annotation(self.xarray[int(value-1)], self.y2array[int(value-1)], "N")
        
        self.plot_1.plot(self.xarray[int(value-1)], self.y1array[int(value-1)], 'ro', markersize=1)
        self.plot_2.plot(self.xarray[int(value-1)], self.y2array[int(value-1)], 'ro', markersize=1)
        
        self.plot_1.canvasShow()
        self.plot_2.canvasShow()
        
        self.previousN = value
    
    def plotShortSpectr(self):
        self.setWindowTitle("Spectrum")
        
        if self.changeParms:
            self.m = self.m_slider.value()
            self.n = self.n_slider.value()
            
            self.plot_1.hide()
            self.plot_2.close()
            self.plot_1.hide()
            self.plot_2.close()
            self.grid.removeWidget(self.plot_1)
            self.grid.removeWidget(self.plot_2)
            self.plot_1.removePolt()
            self.plot_2.removePolt()
            del self.plot_1
            del self.plot_2
            
            self.m_slider.hide()
            self.m_slider.close()
            self.n_slider.hide()
            self.n_slider.close()
            self.grid.removeWidget(self.m_slider)
            self.grid.removeWidget(self.n_slider)
            del self.m_slider
            del self.n_slider
            
            self.plotSmoothDataButton.hide()
            self.plotSmoothDataButton.close()
            self.grid.removeWidget(self.plotSmoothDataButton)
            del self.plotSmoothDataButton
            
            self.m_lcd.hide()
            self.n_lcd.hide()
            self.m_lcd.close()
            self.n_lcd.close()
            self.grid.removeWidget(self.m_lcd)
            self.grid.removeWidget(self.n_lcd)
            del self.m_lcd
            del self.n_lcd
            
            self.mLabel.hide()
            self.mLabel.close()
            self.grid.removeWidget(self.mLabel)
            del self.mLabel
            
            self.nLabel.hide()
            self.nLabel.close()
            self.grid.removeWidget(self.nLabel)
            del self.nLabel
            
            while len(self.infoSet) != 0:
                info_item = self.infoSet.pop()
                info_item.hide()
                info_item.close()
                self.grid.removeWidget(info_item)
                del info_item 
                
            while len(self.infoSet_2) != 0:   
                info_item = self.infoSet_2.pop()
                info_item.hide()
                info_item.close()
                self.grid.removeWidget(info_item)
                del info_item 
                
            del self.infoSet
            del self.infoSet_2
            
            self.changeDataButton.hide()
            self.changeDataButton.close()
            self.grid.removeWidget(self.changeDataButton)
            del self.changeDataButton
                    
            if self.m < 0:
                self.m = 0
        else:
            self.m = 0
            self.n = self.dataPoints
            self.changeParams = QPushButton("Change Params", self)
            self.grid.addWidget(self.changeParams, 4, 3)
            self.changeParams.clicked.connect(self.plotInitData)
            self.changeParams.setStyleSheet("background-color: blue")
            
        self.xarray = self.xarray[self.m:self.n]
        self.y1array = self.y1array[self.m:self.n]
        self.y2array = self.y2array[self.m:self.n]
        
        #u1 plot
        self.plot_10 = Plot()
        self.plot_10.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization", (1, 0))
        self.plot_10.plot(self.xarray, self.y1array, 'ko', label='Data Points',  markersize=1)
        
        #u9 plot
        self.plot_11 = Plot()
        self.plot_11.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization", (1, 1))
        self.plot_11.plot(self.xarray, self.y2array, 'ko', label='Data Points',  markersize=1)
        
        self.grid.addWidget(self.plot_10, 0, 0)
        self.grid.addWidget(self.plot_11, 0, 1)
        
        self.plotPoly = QPushButton("Create Polynomial", self)
        self.grid.addWidget(self.plotPoly, 3, 3)
        self.plotPoly.clicked.connect(self.removeCuts)
        self.plotPoly.setStyleSheet("background-color: green")
        
    def removeCuts(self):
        if not self.changeParms:
            self.changeParams.hide()
            self.changeParams.close()
            self.grid.removeWidget(self.changeParams)
            del self.changeParams
        
        cutsIndex = list()
        cutsIndex.append(0)

        for cut in self.cuts:
            cutsIndex.append((np.abs(self.xarray-float(cut[0]))).argmin()) 
            cutsIndex.append((np.abs(self.xarray-float(cut[1]))).argmin())
            
        cutsIndex.append(self.n)
        
        polyArray_x = list()
        polyArray_u1 = list()
        polyArray_u9 = list()
        
        i = 0
        j = 1
        
        while i != len(cutsIndex):
            polyArray_x.append(self.xarray[cutsIndex[i] : cutsIndex[j]])
            polyArray_u1.append(self.y1array[cutsIndex[i] : cutsIndex[j]])
            polyArray_u9.append(self.y2array[cutsIndex[i] : cutsIndex[j]])
            i = i + 2 
            j = j + 2
            
        poly_x = list()
        poly_u1 = list()
        poly_u9 = list()
        
        for p in polyArray_x:
            for p1 in p:
                poly_x.append(p1)
                
        for p in polyArray_u1:
            for p1 in p:
                poly_u1.append(p1)
        
        for p in polyArray_u9:
            for p1 in p:
                poly_u9.append(p1)
            
        self.polyx = np.array(poly_x)
        self.polyu1 = np.array(poly_u1)
        self.polyu9 = np.array(poly_u9)
        
        z_u1 = np.polyfit(self.polyx, self.polyu1, self.polynomialOrder)
        self.p_u1 = np.poly1d(z_u1)
        
        z_u9 = np.polyfit(self.polyx, self.polyu9, self.polynomialOrder)
        self.p_u9 = np.poly1d(z_u9)
        
        self.plot_10.hide()
        self.plot_11.close()
        self.plot_10.hide()
        self.plot_11.close()
        self.grid.removeWidget(self.plot_10)
        self.grid.removeWidget(self.plot_11)
        self.plot_10.removePolt()
        self.plot_11.removePolt()
        del self.plot_10
        del self.plot_11
        
        self.plotPoly.hide()
        self.plotPoly.close()
        self.grid.removeWidget(self.plotPoly)
        del self.plotPoly
        
        #u1 plot
        self.plot_5 = Plot()
        self.plot_5.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization", (1, 0))
        self.plot_5.plot(self.polyx, self.polyu1, 'ko', label='Data Points',  markersize=1)
        #self.plot_5.plot(self.xarray[self.m:self.n], self.ceb_1(self.xarray[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        self.plot_5.plot(self.xarray, self.p_u1(self.xarray), 'b', label='Numpy polyfit', markersize=1)
        
        #u9 plot
        self.plot_6 = Plot()
        self.plot_6.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization", (1, 1))
        self.plot_6.plot(self.polyx, self.polyu9, 'ko', label='Data Points',  markersize=1)
        #self.plot_6.plot(self.xarray[self.m:self.n], self.ceb_2(self.xarray[self.m:self.n]), 'r', label='Chebyshev polynomial', markersize=1)
        self.plot_6.plot(self.xarray, self.p_u9(self.xarray), 'b', label='Numpy polyfit', markersize=1)
        
        self.grid.addWidget(self.plot_5, 0, 0)
        self.grid.addWidget(self.plot_6, 0, 1)
        
        self.plotPoly = QPushButton("Create local maximums", self)
        self.grid.addWidget(self.plotPoly, 3, 3)
        self.plotPoly.clicked.connect(self.plotLocalMaximum)
        self.plotPoly.setStyleSheet("background-color: green")
    
    def plotLocalMaximum(self):
        self.setWindowTitle("Local maximums")
        
        self.plot_5.hide()
        self.plot_5.close()
        self.plot_6.hide()
        self.plot_6.close()
        self.grid.removeWidget(self.plot_5)
        self.grid.removeWidget(self.plot_6)
        self.plot_5.removePolt()
        self.plot_6.removePolt()
        del self.plot_5
        del self.plot_6
        
        self.plotPoly.hide()
        self.plotPoly.close()
        self.grid.removeWidget(self.plotPoly)
        del self.plotPoly
        
        self.monitoringButton = QPushButton("Add points to monitoring", self)
        self.grid.addWidget(self.monitoringButton, 3, 3)
        self.monitoringButton.clicked.connect(self.createResult)
        self.monitoringButton.setStyleSheet("background-color: green")
        
        self.localMax_Array_u1 = self.y1array - self.p_u1(self.xarray)
        self.localMax_Array_u9 = self.y2array - self.p_u9(self.xarray)
        
        # Smoohting
        if self.skipsmoothing:
            self.z1 = self.localMax_Array_u1
            self.z2 = self.localMax_Array_u9
        else:
            g1 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
            g2 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)

            self.z1 = convolve(self.localMax_Array_u1, g1, boundary='extend')
            self.z2 = convolve(self.localMax_Array_u9, g2, boundary='extend')

        three_sigma_u1 = 3 * np.std(self.polyu1)
        three_sigma_u9 = 3 * np.std(self.polyu9)
        polyuAVG = (self.polyu1 + self.polyu9)/2
        self.avg_y = (self.z1 + self.z2) / 2
        three_sigma_uAVG = 3 * np.std(polyuAVG)

        smart_tres_u1=2.5*three_sigma_u1/np.max(self.z1)
        smart_tres_u9=2.5*three_sigma_u9/np.max(self.z2)
        smart_tres_uAVG=2.5*three_sigma_uAVG/np.max(self.avg_y)

        thres=0.3
           
        #indexsu apreikinasana
        indexes_for_ceb = peakutils.indexes(self.z1, thres=smart_tres_u1, min_dist=3)
        indexes_for_ceb2 = peakutils.indexes(self.z2, thres=smart_tres_u9, min_dist=3)
        indexes_for_avg = peakutils.indexes(self.avg_y, thres=smart_tres_uAVG, min_dist=3)
        
        #u1
        self.plot_7 = Plot()
        self.plot_7.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "1u Polarization", (1, 0))
        self.plot_7.plot(self.xarray, self.z1, 'b', label='Signal - polynomial', markersize=1)
        self.plot_7.plot(self.xarray[indexes_for_ceb], self.z1[indexes_for_ceb], 'dr', label="Local Maximums for signal", markersize=2)
        self.plot_7.annotations(self.xarray[indexes_for_ceb], self.z1[indexes_for_ceb])
        
        #u9
        self.plot_8 = Plot()
        self.plot_8.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "9u Polarization", (1, 1))
        self.plot_8.plot(self.xarray, self.z2, 'b', label='Signal - polynomial', markersize=1)
        self.plot_8.plot(self.xarray[indexes_for_ceb2], self.z2[indexes_for_ceb2], 'dr', label="Local Maximums for signal", markersize=2)
        self.plot_8.annotations(self.xarray[indexes_for_ceb2], self.z2[indexes_for_ceb2])
        
        #uAVG
        self.plot_9 = Plot()
        self.plot_9.creatPlot(self.grid, 'Velocity (km sec$^{-1}$)', 'Flux density (Jy)', "Average Polarization", (1, 2))
        self.plot_9.plot(self.xarray, self.avg_y, 'b', label='Signal - polynomial', markersize=1)
        self.plot_9.plot(self.xarray[indexes_for_avg], self.avg_y[indexes_for_avg], 'dr', label="Local Maximums for signal", markersize=2)
        self.plot_9.annotations(self.xarray[indexes_for_avg], self.avg_y[indexes_for_avg])
        
        self.grid.addWidget(self.plot_7, 0, 0)
        self.grid.addWidget(self.plot_8, 0, 1)
        self.grid.addWidget(self.plot_9, 0, 2)
        
        self.maxU1 = list()
        self.maxU9 = list()
        self.avgMax = list()
        self.maxu1_index = list()
        self.maxu9_index = list()
        self.maxavg_index = list()
        
    def createResult(self):
        resultFileName = self.source + ".json"
        
        if os.path.isfile(self.resultFilePath + resultFileName):
            pass
        else:
            os.system("touch " + self.resultFilePath +  resultFileName)
            
            resultFile = open (self.resultFilePath +  resultFileName, "w")
            resultFile.write("{ \n" + "\n}")
            resultFile.close()
        
        with open(self.resultFilePath + resultFileName) as result_data:    
            result = json.load(result_data)
        
        if self.expername not in result:
            result[self.expername] = dict()
        
        indexies_for_source_velocities = [0] * len(self.source_velocities)
        
        for index in range (0, len(self.source_velocities)):
            indexies_for_source_velocities[index] =  (np.abs(self.xarray-float(self.source_velocities[index]))).argmin()
        print (indexies_for_source_velocities)
        
        max_amplitude_list_u1 = list()
        max_amplitude_list_u9 = list()
        max_amplitude_list_uavg = list()
        for index in indexies_for_source_velocities:
            max_amplitude_list_tmp_u1 = list()
            max_amplitude_list_tmp_u9 = list()
            max_amplitude_list_tmp_uavg = list()
            for i in range (index - self.index_range_for_local_maxima, index + self.index_range_for_local_maxima):
                max_amplitude_list_tmp_u1.append(self.z1[i])
                max_amplitude_list_tmp_u9.append(self.z2[i])
                max_amplitude_list_tmp_uavg.append(self.avg_y[i])
                
            max_amplitude_list_u1.append(max_amplitude_list_tmp_u1)
            max_amplitude_list_u9.append(max_amplitude_list_tmp_u9)
            max_amplitude_list_uavg.append(max_amplitude_list_tmp_uavg)
        
        max_apmlitudes_u1 = [np.max(value) for value  in max_amplitude_list_u1]
        max_apmlitudes_u9 = [np.max(value) for value  in max_amplitude_list_u9]
        max_apmlitudes_uavg = [np.max(value) for value  in max_amplitude_list_uavg]
        
        for max in range(0,len(max_apmlitudes_u1)):
            max_apmlitudes_u1[max] = [self.source_velocities[max], max_apmlitudes_u1[max]]
            max_apmlitudes_u9[max] = [self.source_velocities[max], max_apmlitudes_u9[max]]
            max_apmlitudes_uavg[max] = [self.source_velocities[max], max_apmlitudes_uavg[max]]
            
        time = self.time + self.date.replace(" ", "_")
        
        try:
            time=datetime.strptime(time, "%H:%M:%S%d_%b_%Y" ).isoformat()
            t=Time(time, format='isot')
            MJD = t.mjd
            result[self.expername]["modifiedJulianDays"] = MJD
        except ValueError as e:
            print ("Cannot crate modified Julian Days",  e)
        except:
            print("Unexpected error:", sys.exc_info()[0])                                       
            
        result[self.expername]["location"] = self.location
        result[self.expername]["Date"] = self.date
        result[self.expername]["Iteration_number"] = int(self.iteration_number)
        result[self.expername]["time"] = self.time
        result[self.expername]["specie"] = self.specie
        
        result[self.expername]["polarizationU1"] =  max_apmlitudes_u1
        result[self.expername]["polarizationU9"] = max_apmlitudes_u9
        result[self.expername]["polarizationAVG"] = max_apmlitudes_uavg
                
        resultFile = open (self.resultFilePath +  resultFileName, "w")
        resultFile.write(json.dumps(result, indent=2))
        resultFile.close()
        
        totalResults = [self.xarray,  self.z1,  self.z2,  self.avg_y]
        output_file_name = self.output + self.source + "_" + self.time.replace(":", "_") + "_" + self.date.replace(" ", "_") + "_" + self.location + "_" + str(self.iteration_number) + ".dat"
        output_file_name = output_file_name.replace(" ", "")
        np.savetxt(output_file_name, np.transpose(totalResults)) 
        
        self._quit()
    
    @QtCore.pyqtSlot()   
    def _quit(self):
        for i in reversed(range(self.grid.count())): 
            self.grid.itemAt(i).widget().deleteLater()
        
        self.hide()
        self.close()
        del self
              
def main():
    args = parseArguments()
   
    datafile = str(args.__dict__["datafile"])
    configFilePath = str(args.__dict__["config"])

    if args.rawdata:
        skipsmooth = True
    else:
        skipsmooth = False
    
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    dataFilesPath =  config.get('paths', "dataFilePath")
    resultFilePath =  config.get('paths', "resultFilePath")
    output =  config.get('paths', "outputFilePath")
    source = datafile.split("/")[-1].split(".")[0].split("_")[0]
    cuts = config.get('cuts', source).split(";")
    cuts = [c.split(",") for c in  cuts]
    print ("source", source)
    source_velocities = config.get('velocities', source).replace(" ", "").split(",")
    index_range_for_local_maxima = int(config.get('parameters', "index_range_for_local_maxima"))
    
    #Create App
    qApp = QApplication(sys.argv)

    aw = Analyzer(dataFilesPath + datafile, resultFilePath, source_velocities, cuts, output, index_range_for_local_maxima, skipsmooth)
    aw.show()
    aw.showMaximized() 
    sys.exit(qApp.exec_())
    
    sys.exit(0)
                
if __name__=="__main__":
    main()