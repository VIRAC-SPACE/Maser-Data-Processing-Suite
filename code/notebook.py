#! /usr/bin/env python

from tkinter import *
from tkinter.ttk import *
import tkinter as tk
import tkinter.ttk as tkk

from ploting import Plot

def frame(parent, size, sides, **options):
    Width=size[0]
    Height=size[1]
    f=tk.Frame(parent, width=Width,height=Height,**options)
    f.pack(side = sides)
    return (f)

class App(tk.Frame):

    def __init__(self, window):

        tk.Frame.__init__(self)
        self.window = window
        self.noteBook = tkk.Notebook(self.window,  name='master')
	
        self.infoFrame = frame(self.window,(1000,1000), None)
        #self.masterFrame = tkk.Frame(self.noteBook)
        #self.plotFrame = tkk.Frame(self.masterFrame)
        self.masterFrame = frame(self.window,(1000,1000), None)
        self.plotFrame = frame(self.masterFrame,(1000,1000), None)

        self.plot_1 = Plot(6,6, self.plotFrame, self.masterFrame)
        self.plot_1.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "1u Polarization")
        self.plot_1.plot([1,2,3], [1,2,3], 'ko', 'Data Points', 5, 15)
        self.plot_1.addPickEvent(self.onpickU1)
        self.plot_1.addSecondAss("x", "Data points", 0, 4096 + 512, 1024)

        self.plot_2 = Plot(6,6, self.plotFrame, self.masterFrame)
        self.plot_2.creatPlot(LEFT, 'Frequency Mhz', 'Flux density (Jy)', "9u Polarization")
        self.plot_2.plot([1,2,3], [1,2,3], 'ko', 'Data Points', 5, 15)
        self.plot_2.addPickEvent(self.onpickU1)
        self.plot_2.addSecondAss("x", "Data points", 0, 4096 + 512, 1024)

        self.noteBook.add(self.infoFrame, text='Info')
        self.noteBook.add(self.masterFrame, text='Data')
        self.noteBook.pack(fill=BOTH)
	
    def onpickU1():
        print ("pick")
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        p = tuple(zip(xdata[ind], ydata[ind]))
        self.plot_1.plot(p[0][0], p[0][1], 'ro', None, 6, 5)
        print (p[0][0], p[0][1])
        self.plot_1.canvasShow()

def main(): 
    #Create App
    window = tk.Tk() 
    ploting = App(window)
    ploting.mainloop()
    
    sys.exit(0)

if __name__=="__main__":
    #if len(sys.argv) < 3:
        #usage()
        #sys.exit(1)
        
    main()
