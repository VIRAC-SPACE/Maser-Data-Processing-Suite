from __future__ import unicode_literals
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.backends.backend_qt5agg as qt5agg
from matplotlib.widgets import Slider
from matplotlib.figure import Figure
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Time New Roman']
rcParams['font.size'] = 10
import mplcursors
from PyQt5 import QtWidgets

class Plot(FigureCanvas):
    
    def __init__(self, parent=None, width=7, height=7):
        self.parent = parent
        self.fig = Figure(figsize=(width, height))
        FigureCanvas.__init__(self, self.fig)
        self.setParent(self.parent)
        FigureCanvas.setSizePolicy(self, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    def plot(self, x, y, line, **options):
        line = self.graph.plot(x,y, line, **options)
        self.graph.legend()
        return line


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
        
    def get_label(self):
        return self.graph.get_label()
    
    def get_visible(self):
        return self.graph.get_visible()
        
    def setXtics(self, x, y, rotation, **options):
        self.graph.set_xticks(x, y, **options)
        self.graph.set_xticklabels(x, rotation=rotation, **options)
       
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
    
    def addPickEvent(self, callback):
        self.cidPick = FigureCanvas.mpl_connect(self, 'pick_event', callback)
    
    def removePickEvent(self):
        FigureCanvas.mpl_disconnect(self, self.cidPick)
        
    def addKeyPressEvent(self, callback):
        self.cidKeyPress = FigureCanvas.mpl_connect(self, 'key_press_event', callback)
    
    def removeKeyPressEvent(self):
        FigureCanvas.mpl_disconnect(self, self.cidKeyPress)
        
    def addSecondAxis(self, axiss, label, start, stop, step):
        self.second_x_axis = self.graph.twiny()
        self.second_x_axis.set_xlabel(label)
        self.graph.tick_params(axis=axiss)
        self.second_x_axis.set_xticks(range(start, stop, step))
        
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
        
        #del self.canvas
        del self.toolbar
        del self
