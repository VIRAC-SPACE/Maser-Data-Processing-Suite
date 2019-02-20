import numpy as np
from astropy.time import Time

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def indexies(array, value):
    indexs = list()
    for i in range(0, len(array)):
        if array[i] == value:
            indexs.append(i)
    return indexs

def correctNumpyReadData(data):
    correctedData = list()
    for d in data:
        correctedData.append(d[0])
    return np.array(correctedData)

def convertDatetimeObjectToMJD(time):
    time=time.isoformat()
    t=Time(time, format='isot')
    return t.mjd