import numpy as np
from matplotlib import pyplot as plt

from help import *
from parsers._configparser import ConfigParser

def getConfigs(key, value):
    configFilePath = "config/config.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

file1 = getConfigs("paths", "notSmoohtFilePath") + "cepa/" + "cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat" #Befor convolation
file2 = getConfigs("paths", "notSmoohtFilePath") + "cepa/" + "cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat.new" #After deconvulatio
file3 = getConfigs("paths", "outputFilePath") + "cepa/" + "cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat" #After convolation

data1 = np.fromfile(file1, dtype="float64", count=-1, sep=" ").reshape((file_len(file1), 4))
data2 = np.fromfile(file2, dtype="float64", count=-1, sep=" ").reshape((file_len(file2), 4))
data3 = np.fromfile(file3, dtype="float64", count=-1, sep=" ").reshape((file_len(file3), 4))

v1 = correctNumpyReadData(data1[:, [0]])
amp1 = correctNumpyReadData(data1[:, [3]])

v2 = correctNumpyReadData(data2[:, [0]])
amp2 = correctNumpyReadData(data2[:, [3]])

v3 = correctNumpyReadData(data3[:, [0]])
amp3 = correctNumpyReadData(data3[:, [3]])

plt.subplot(1,3,1)
plt.plot(v1, amp1)
plt.title("Befor convolation")

plt.subplot(1,3,2)
plt.plot(v2,amp2)
plt.title("After deconvulatio")

plt.subplot(1,3,3)
plt.plot(v3,amp3)
plt.title("After convolation")

plt.show()


