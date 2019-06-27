import matplotlib.pyplot as plt
import numpy as np

from help import *
from parsers._configparser import ConfigParser

def getConfigs(key, value):
    configFilePath = "config/config.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)


file = "monitoring/cepa.txt"
compunetCount = len(getConfigs("velocities", "cepa").replace(" ", "").split(","))
data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),compunetCount + 1))
x = correctNumpyReadData(data[:, [0]])
y = correctNumpyReadData(data[:, [1]])
plt.plot(x,y)
plt.show()