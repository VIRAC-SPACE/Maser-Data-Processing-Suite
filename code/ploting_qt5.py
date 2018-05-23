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
    
    def __init__(self, parent=None, width=7, height=7):
        self.parent = parent
        self.fig = Figure(figsize=(width, height))
        FigureCanvas.__init__(self, self.fig)
        self.setParent(self.parent)
        self.canvas = FigureCanvas(self.fig)
        ##self.canvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        #self.canvas.updateGeometry(self)
    
    def plot(self, x, y, line, **options):
        self.graph.plot(x,y, line, **options)
        
        self.graph.legend()
        
    def creatPlot(self, grid, x_label, y_label, title,toolbarpos):
        self.graph = self.fig.add_subplot(111)
        self.grid = grid
        
        self.x_label = x_label
        self.y_label = y_label
        self.title = title
        
        if self.title != None:
            self.graph.set_title(title,  y=1.08) 
        
        self.toolbar = qt5agg.NavigationToolbar2QT(self, self.parent)
        self.toolbar.update()       
        self.grid.addWidget(self.toolbar, toolbarpos[0], toolbarpos[1])
        
        self.graph.grid(True)
        self.graph.set_xlabel(x_label)
        self.graph.set_ylabel(y_label)
        
    def annotations(self, xvalues, yvalues):
        ax = self.figure.add_subplot(111)
        for xy in zip(xvalues, yvalues):                        
            ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
            
    def annotation(self, xvalue, yvalue, text):
        self.ax = self.figure.add_subplot(111)
        self.annotate = self.ax.annotate(text, xy=(xvalue, yvalue),  textcoords='data')
    
    def remannotation(self):
        self.annotate.remove()
        del self.annotate
        
    def canvasShow(self):
        self.canvas.draw()
    
    def addPickEvent(self, callback):
        self.cid = self.canvas.mpl_connect('pick_event', callback)
    
    def removePickEvent(self):
        self.figure.canvas.mpl_disconnect(self.cid)
        
    def addSecondAxis(self, axiss, label, start, stop, step):
        self.second_x_axis = self.graph.twiny()
        self.second_x_axis.set_xlabel(label)
        self.graph.tick_params(axis=axiss)
        self.second_x_axis.set_xticks(range(start, stop, step))
    
    def removePolt(self):
        self.fig.clf()
        del self.graph
        self.grid.removeWidget(self.toolbar)
        self.toolbar.hide()
        self.toolbar.close()
        self.toolbar.destroy()
        self.canvas.destroy()
        
    def __del__(self):
        del self.fig
        
        del self.canvas
        del self.toolbar
        del self
