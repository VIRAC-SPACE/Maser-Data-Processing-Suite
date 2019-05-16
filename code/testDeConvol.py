import numpy as np
from matplotlib import pyplot as plt

from help import *
from parsers._configparser import ConfigParser

def getConfigs(key, value):
    configFilePath = "config/config.cfg"
    config = ConfigParser.getInstance()
    config.CreateConfig(configFilePath)
    return config.getConfig(key, value)

file1 = getConfigs("paths", "notSmoohtFilePath") + "cepa/" + "cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat"
file2 = getConfigs("paths", "notSmoohtFilePath") + "cepa/" + "cepa_23_53_48_22_Dec_2018_IRBENE16_294.dat.new"

data1 = np.fromfile(file1, dtype="float64", count=-1, sep=" ").reshape((file_len(file1), 4))
data2 = np.fromfile(file2, dtype="float64", count=-1, sep=" ").reshape((file_len(file2), 4))

v1 = correctNumpyReadData(data1[:, [0]])
amp1 = correctNumpyReadData(data1[:, [3]])

v2 = correctNumpyReadData(data2[:, [0]])
amp2 = correctNumpyReadData(data2[:, [3]])

print(len(v2), len(amp2))

plt.subplot(1,2,1)
plt.plot(v1, amp1)
plt.subplot(1,2,2)
plt.plot(v2,amp2)
plt.show()


