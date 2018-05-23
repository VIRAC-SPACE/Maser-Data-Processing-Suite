from __future__ import unicode_literals
import sys
import os
import matplotlib

matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Time New Roman']
rcParams['font.size'] = 12

class Plot(FigureCanvas):
    
    def __init__(self, parent=None, width=5, height=5, dpi=100):
        self.parent = parent
        fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, fig)
        self.setParent(self.parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def plot(self, x, y, line, **options):
        self.graph.plot(x,y, line, **options)
        self.graph.legend()
        
    def creatPlot(self, sides, x_label, y_label, title):
        self.graph = self.figure.add_subplot(111)
        
        self.x_label = x_label
        self.y_label = y_label
        self.title = title
        
        if title != None:
            self.graph.set_title(title,  y=1.08) 
                
        #self.toolbar = tkagg.NavigationToolbar2TkAgg(self.canvas, self.window)
        #self.toolbar.update()
        #self.canvas._tkcanvas.pack(side=LEFT, expand=YES)
        self.graph.grid(True)
        self.graph.set_xlabel(x_label)
        self.graph.set_ylabel (y_label)
    
    def removePolt(self):
        self.figure.clf()
        del self.graph
        #FigureCanvas.destroy()
        #self.toolbar.destroy()
        
        
