from __future__ import unicode_literals
import sys
import os
import matplotlib

matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.backends.backend_qt5agg as qt5agg
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Time New Roman']
rcParams['font.size'] = 12

class Plot(FigureCanvas):
    
    def __init__(self, parent=None, width=5, height=5, dpi=100):
        self.parent = parent
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        FigureCanvas.__init__(self, self.fig)
        self.setParent(self.parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    
    def plot(self, x, y, line, **options):
        self.graph.plot(x,y, line, **options)
        self.graph.legend()
        
    def creatPlot(self, grid, x_label, y_label, title,toolbarpos):
        self.graph = self.figure.add_subplot(111)
        self.grid = grid
        
        self.x_label = x_label
        self.y_label = y_label
        self.title = title
        
        if title != None:
            self.graph.set_title(title,  y=1.08) 
        
        self.toolbar = qt5agg.NavigationToolbar2QT(self, self.parent)
        self.toolbar.update()       
        self.grid.addWidget(self.toolbar, toolbarpos[0], toolbarpos[1])
        
        self.graph.grid(True)
        self.graph.set_xlabel(x_label)
        self.graph.set_ylabel (y_label)
    
    def removePolt(self):
        self.figure.clf()
        self.fig.clf()
        del self.graph
        #FigureCanvas.destroy()
        self.grid.removeWidget(self.toolbar)
        self.toolbar.hide()
        self.toolbar.close()
        self.toolbar.destroy()
