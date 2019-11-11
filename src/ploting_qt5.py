from __future__ import unicode_literals
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.backends.backend_qt5agg as qt5agg
from matplotlib.widgets import Slider
from matplotlib.figure import Figure
from matplotlib import rcParams

from parsers._configparser import ConfigParser

def getConfigsItems():
    configFilePath = "config/plot.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getItems("main")

configItems = getConfigsItems()
for key, value in configItems.items():
    rcParams[key] = value

import mplcursors
import numpy as np
from PyQt5 import QtWidgets

class Plot(FigureCanvas):
    
    def __init__(self, parent=None, width=16, height=9):
        self.parent = parent
        self.fig = Figure(figsize=(width, height))
        FigureCanvas.__init__(self, self.fig)
        self.setParent(self.parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def plot(self, x, y, line, fontsize=12, **options):
        line = self.graph.plot(x,y, line, **options)
        #self.graph.set_ylim(np.min(y), np.max(y))
        #self.graph.set_ylim(0, 1)

        if self.type == "log":
            self.graph.set_yscale("log")
            yTicks = []
            y_min = np.min(y)
            if y_min < 0:
                ymin = 10

            ytick = 1
            n = 1
            while ytick < np.max(y) + 100:
                if ytick > y_min:
                    yTicks.append(ytick)
                n +=1
                ytick = 10**n

            yTicks.append(np.max(y))
            self.graph.set_yticks(yTicks)

        elif self.type == "linear":
            self.graph.set_yscale("linear")

        else:
            print("Wrong plot type")

        if "label" in options.keys():
            self.legend = self.graph.legend(prop={'size': fontsize})
            self.legend.set_draggable(True, update='loc')
        return line

    def errorbar(self,  x, y, error, line,  **options):
        self.graph.errorbar(x, y, yerr=error, fmt=line, **options)

    def set_tick_params(self, axis, direction, which, length, width, labelsize, rotation):
        self.graph.tick_params(axis=axis, direction=direction, which=which, length=length, width=width, labelsize=labelsize, rotation=rotation)

    def save_fig(self, fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None, metadata=None):
        self.fig.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None, metadata=None)

    def get_xlim(self):
        return self.graph.get_xlim()

    def set_xlim(self, new_x_lim):
        self.graph.set_xlim(new_x_lim)

    def get_ylim(self):
        return self.graph.get_ylim()

    def set_ylim(self, new_y_lim):
        self.graph.set_ylim(new_y_lim)
    
    def contourf(self, x, y, z, **options):
        cs = self.graph.contourf(x, y, z, 500, cmap='jet', **options)
        return cs

    def remove_markers(self):
        for line in self.graph.lines:
            line.set_marker(None)
    
    def colorbar(self, cs):
        cbar = self.fig.colorbar(cs)
        return cbar

    def creatPlot(self, grid, x_label, y_label, title, toolbarpos, type):
        self.graph = self.fig.add_subplot(111)
        #self.graph2 = self.fig.add_subplot(111, sharex=self.graph, frameon=False)

        self.type = type

        '''
        self.graph2.yaxis.tick_right()
        self.graph2.yaxis.set_label_position("right")
        self.graph2.set_yscale("linear")
        self.graph2.set_yticks([])
        mn, mx = self.graph.get_ylim()
        self.graph2.set_ylim(mn, mx)
        '''

        self.grid = grid
        
        self.x_label = x_label
        self.y_label = y_label
        self.title = title
        
        if self.title != None:
            self.graph.set_title(title,  y=1.08) 
        
        self.toolbar = qt5agg.NavigationToolbar2QT(self, self.parent)
        self.toolbar.update()
        self.grid.addWidget(self.toolbar, toolbarpos[0], toolbarpos[1])
        
        #self.graph.grid(True, which='major', color='k', linestyle='-', linewidth=0.5)
        #self.graph.grid(False, which='minor')

        self.graph.yaxis.set_ticks_position('both')
        self.graph.xaxis.set_ticks_position('both')
        self.fig.tight_layout(pad=2.3, h_pad=None, w_pad=None, rect=None)

        #self.graph.autoscale()

        if self.x_label != None:
            self.graph.set_xlabel(self.x_label)
        if self.y_label !=None:
            self.graph.set_ylabel(self.y_label)
        
    def get_label(self):
        return self.graph.get_label()

    def setAxiesText(self, x, y, s,  **options):
        self.graph.axes.text(x, y, s,  **options)

    def get_axes(self):
        return self.graph.axes
    
    def get_visible(self):
        return self.graph.get_visible()
        
    def setXtics(self, x, y, rotation, **options):
        self.graph.set_xticks(x, y, **options)
        self.graph.set_xticklabels(y, rotation=rotation, **options)
       
    def annotations(self, xvalues, yvalues):
        ax = self.figure.add_subplot(111)
        for xy in zip(xvalues, yvalues):                        
            ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
            
    def addCursor(self, labels):
        cursor =  mplcursors.cursor(self.graph, hover=True, highlight=True)
        cursor.connect("add", lambda sel: sel.annotation.set_text(labels[sel.target.index]))
            
    def annotation(self, xvalue, yvalue, text):
        self.ax = self.figure.add_subplot(111)
        self.annotate = self.ax.annotate(text, xy=(xvalue, yvalue),  textcoords='data')
         
    def remannotation(self):
        self.annotate.remove()
        del self.annotate
        
    def canvasShow(self):
        FigureCanvas.draw(self)

    def addZoomEvent(self, callback):
        self.zoom = self.graph.callbacks.connect('xlim_changed', callback)
        self.zoom = self.graph.callbacks.connect('ylim_changed', callback)
    
    def addPickEvent(self, callback):
        self.cidPick = FigureCanvas.mpl_connect(self, 'pick_event', callback)

    def addClickEvent(self, callback):
        self.cidClick = FigureCanvas.mpl_connect(self, 'button_press_event', callback)

    def removePickEvent(self):
        FigureCanvas.mpl_disconnect(self, self.cidPick)

    def removeClickEvent(self):
        FigureCanvas.mpl_disconnect(self, self.cidClick)
        
    def addKeyPressEvent(self, callback):
        self.cidKeyPress = FigureCanvas.mpl_connect(self, 'key_press_event', callback)
    
    def removeKeyPressEvent(self):
        FigureCanvas.mpl_disconnect(self, self.cidKeyPress)
        
    def addSecondAxis(self, axiss, label, start, stop, step):
        self.second_x_axis = self.graph.twiny()
        self.second_x_axis.set_xlabel(label)
        self.graph.tick_params(axis=axiss)
        self.second_x_axis.set_xticks(range(start, stop, step))

    def addSecondAxis2(self, values, label, axiss):
        self.second_x_axis = self.graph.twiny()
        self.second_x_axis.set_xlabel(label)
        self.graph.tick_params(axis=axiss)
        self.second_x_axis.set_xticklabels(values)
        
    def addSlider(self, cords, label, start, stop, init, callback):
        self.figure.subplots_adjust(bottom=0.25)
        axcolor = 'lightgoldenrodyellow'
        Axes = self.figure.add_axes(cords, axisbg=axcolor)
        slider=Slider(Axes, label,  start, stop, valinit=init)
        slider.on_changed(callback)
        
    def removePolt(self):
        self.fig.clf()
        del self.graph
        self.grid.removeWidget(self.toolbar)
        self.toolbar.hide()
        self.toolbar.close()
        self.toolbar.destroy()
        
    def __del__(self):
        del self.fig
        del self.toolbar
        del self
