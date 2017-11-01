#! /usr/bin/python

import os
import sys
import subprocess
import ast
import numpy as np
from matplotlib import pyplot as plt

from experimentsLogReader import *

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def usage():
    print ('Usage: '+sys.argv[0]+' log file')
    
if __name__=="__main__":
    if len(sys.argv) < 3:
        usage()
        sys.exit(1)

    proc = subprocess.Popen(['python', './experimentsLogReader.py',  sys.argv[1]], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    logfile =  ast.literal_eval(proc.communicate()[0])

    data = np.fromfile(sys.argv[2], dtype="float64", count=-1, sep=" ") .reshape((file_len(sys.argv[2]),5))
    data = np.delete(data, (0), axis=0)
    
    scanNumber = sys.argv[2].split(".")[0].split("_")[1][1:len(sys.argv[2])]
    scan = logfile[scanNumber]
    Systemtemperature1u = float(scan["Systemtemperature"][0])
    Systemtemperature9u = float(scan["Systemtemperature"][1])
    
    plt.plot(data[:, [0]], data[:, [1]] * Systemtemperature1u)
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("1u")
    plt.show()
    
    plt.plot(data[:, [0]], data[:, [2]] * Systemtemperature9u)
    plt.grid(True)
    plt.xlabel('Mhz')
    plt.legend("9u")
    plt.show()
    sys.exit(0)