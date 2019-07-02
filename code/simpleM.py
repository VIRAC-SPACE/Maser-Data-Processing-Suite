import matplotlib.pyplot as plt
import numpy as np

from help import *
from parsers._configparser import ConfigParser

def getConfigs(key, value):
    configFilePath = "config/config.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

components = [1,2,3,4,5]
file = "monitoring/cepa.txt"
compunetCount = len(getConfigs("velocities", "cepa").replace(" ", "").split(","))
velocity = getConfigs("velocities", "cepa").replace(" ", "").split(",")
data = np.fromfile(file, dtype="float64", count=-1, sep=" ").reshape((file_len(file),compunetCount + 1))
x = correctNumpyReadData(data[:, [0]])

for component in components:
    index = components.index(component)
    plt.plot(x,correctNumpyReadData(data[:, [index +1]]), label=velocity[index])
    plt.legend()
    plt.grid(True)
plt.show()