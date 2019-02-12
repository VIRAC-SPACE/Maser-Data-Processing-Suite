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
