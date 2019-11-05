import sys
import matplotlib

matplotlib.use('Qt5Agg')
import matplotlib.pyplot  as plt
import datetime
import json
import argparse
import configparser
from operator import itemgetter
import numpy as np

from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QPushButton, QLabel, QLineEdit, QDesktopWidget, QComboBox)
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

from ploting_qt5 import  Plot
from result import  Result

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
from PyQt5.QtWidgets import (QWidget, QGridLayout, QApplication, QPushButton, QMessageBox, QLabel, QLineEdit, QSlider, QDesktopWidget, QLCDNumber)
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QIcon

from PyQt5.QtGui import QColor

from ploting_qt5 import  Plot