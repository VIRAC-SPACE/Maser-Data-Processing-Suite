import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.backends.backend_tkagg as tkagg
from matplotlib.widgets import Slider

try:
    from tkinter import *
    import tkinter as tk
except:
    from Tkinter import *
    import Tkinter as tk

class Plot():
    
    def __init__(self, fig_size_x, fig_size_y, window, frame):
        self.window = window
        self.frame = frame
        
        self.fig_size_x = fig_size_x
        self.fig_size_y = fig_size_y
    
    def creatPlot(self, sides, x_label, y_label, title):
        self.figure = Figure(figsize=(self.fig_size_x,self.fig_size_y))
        self.graph = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.show()
        self.figure.set_canvas(self.canvas)
        self.canvas.get_tk_widget().pack(side=sides)
        self.graph.set_title(title,  y=1.08) 
        
        self.toolbar = tkagg.NavigationToolbar2TkAgg(self.canvas, self.window)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=LEFT, expand=1)
        self.graph.grid(True)
        self.graph.set_xlabel(x_label)
        self.graph.set_ylabel (y_label)
        
    def plot(self, x, y, line, labels, markersizes, pickers):
        self.graph.plot(x, y, line, label=labels, markersize=markersizes, picker=pickers)
        self.graph.legend(loc=2)
        
    def annotations(self, xvalues, yvalues):
        ax = self.figure.add_subplot(111)
        for xy in zip(xvalues, yvalues):                        
            ax.annotate('(%.2f, %.1f)' % xy, xy=xy, textcoords='data')
            
    def annotation(self, xvalue, yvalue, text):
        self.ax = self.figure.add_subplot(111)
        self.annotate = self.ax.annotate(text, xy=(xvalue, yvalue),  textcoords='data')

    def remannotation(self):
        self.annotate.remove()
        
    def canvasShow(self):
        self.canvas.show()
    
    def addPickEvent(self, callback):
        self.cid = self.canvas.mpl_connect('pick_event', callback)
    
    def removePickEvent(self):
        self.figure.canvas.mpl_disconnect(self.cid)
        
    def addSecondAss(self, ass, label, start, stop, step):
        self.second_x_ass = self.graph.twiny()
        self.second_x_ass.set_xlabel(label)
        self.graph.tick_params(axis=ass)
        self.second_x_ass.set_xticks(range(start, stop, step))
        
    def addSlider(self, cords, label, start, stop, init, callback):
        self.figure.subplots_adjust(bottom=0.25)
        axcolor = 'lightgoldenrodyellow'
        Axes = self.figure.add_axes(cords, axisbg=axcolor)
        slider=Slider(Axes, label,  start, stop, valinit=init)
        slider.on_changed(callback)
        
    def removePolt(self):
        self.figure.clf()
        self.canvas.get_tk_widget().delete("all")
        self.canvas.get_tk_widget().destroy()
        self.toolbar.destroy()
        
