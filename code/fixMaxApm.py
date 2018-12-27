#! /usr/bin/python3
# -*- coding: utf-8 -*-

import json
import argparse
import configparser
import pickle
import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve

def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser(description='''Fix Amplitudes in output file. ''',
    epilog="""Fix.""")
    
    # Positional mandatory arguments
    parser.add_argument("outputfilename", help="output file name", type=str, default="")

    # Positional mandatory arguments
    parser.add_argument("-c", "--config", help="Configuration cfg file", type=str, default="config/config.cfg")

    # Print version
    parser.add_argument("-v","--version", action="version", version='%(prog)s - Version 0.1')

    # Parse arguments
    args = parser.parse_args()

    return args


class Result():
    def __init__(self, matrix, specie):
        self.matrix = matrix
        self.specie = specie
        
    def getMatrix(self):
        return self.matrix
    
    def getSpecie(self):
        return self.specie

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1    
    
def getLocalMax(data_file, source_velocities, index_range_for_local_maxima, cuts):
    result = pickle.load(open(data_file, "rb"))
    data = result.getMatrix()
    dataPoints = data.shape[0]
    
    x = data[:, [0]]
    u1 = data[:, [1]] 
    u9 = data[:, [2]]
    xarray = np.zeros(dataPoints)
    y1array = np.zeros(dataPoints)
    y2array = np.zeros(dataPoints)
    
    for i in range (0, dataPoints):
        xarray[i] = x[i]
        y1array[i] = u1[i]
        y2array[i] = u9[i]
        
    xarray =  np.flip(xarray,0)
    y1array =  np.flip(y1array,0)
    y2array =  np.flip(y2array,0)
    
    cutsIndex = list()
    cutsIndex.append(0)

    for cut in cuts:
        cutsIndex.append((np.abs(xarray-float(cut[0]))).argmin()) 
        cutsIndex.append((np.abs(xarray-float(cut[1]))).argmin())
            
        cutsIndex.append(dataPoints)
    
    polyArray_x = list()
    polyArray_u1 = list()
    polyArray_u9 = list()
        
    i = 0
    j = 1
        
    while i != len(cutsIndex):
        polyArray_x.append(xarray[cutsIndex[i] : cutsIndex[j]])
        polyArray_u1.append(y1array[cutsIndex[i] : cutsIndex[j]])
        polyArray_u9.append(y2array[cutsIndex[i] : cutsIndex[j]])
        i = i + 2 
        j = j + 2
    
    poly_x = list()
    poly_u1 = list()
    poly_u9 = list()
        
    for p in polyArray_x:
        for p1 in p:
            poly_x.append(p1)
                
    for p in polyArray_u1:
        for p1 in p:
            poly_u1.append(p1)
        
    for p in polyArray_u9:
        for p1 in p:
            poly_u9.append(p1)
    
    polyx = np.array(poly_x)
    polyu1 = np.array(poly_u1)
    polyu9 = np.array(poly_u9)
    
    z_u1 = np.polyfit(polyx, polyu1, 3)
    p_u1 = np.poly1d(z_u1)
        
    z_u9 = np.polyfit(polyx, polyu9, 3)
    p_u9 = np.poly1d(z_u9)
    
    localMax_Array_u1 = y1array - p_u1(xarray)
    localMax_Array_u9 = y2array - p_u9(xarray)
      
    g1 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
    g2 = Gaussian1DKernel(stddev=3, x_size=19, mode='center', factor=100)
    z1 = convolve(localMax_Array_u1, g1, boundary='extend')
    z2 = convolve(localMax_Array_u9, g2, boundary='extend')
    avg_y = (z1 + z2) / 2
    
    indexies_for_source_velocities = [0] * len(source_velocities)
    for index in range (0, len(source_velocities)):
            indexies_for_source_velocities[index] =  (np.abs(x-float(source_velocities[index]))).argmin()
                      
    max_amplitude_list_u1 = list()
    max_amplitude_list_u9 = list()
    max_amplitude_list_uavg = list()

    for index in indexies_for_source_velocities:
            max_amplitude_list_tmp_u1 = list()
            max_amplitude_list_tmp_u9 = list()
            max_amplitude_list_tmp_uavg = list()
           
            for i in range (index - index_range_for_local_maxima, index + index_range_for_local_maxima):
                max_amplitude_list_tmp_u1.append(z1[i])
                max_amplitude_list_tmp_u9.append(z2[i])
                max_amplitude_list_tmp_uavg.append(avg_y[i])
                  
            max_amplitude_list_u1.append(max_amplitude_list_tmp_u1)
            max_amplitude_list_u9.append(max_amplitude_list_tmp_u9)
            max_amplitude_list_uavg.append(max_amplitude_list_tmp_uavg)
        
    max_apmlitudes_u1 = [np.max(value) for value  in max_amplitude_list_u1]
    max_apmlitudes_u9 = [np.max(value) for value  in max_amplitude_list_u9]
    max_apmlitudes_uavg = [np.max(value) for value  in max_amplitude_list_uavg]
        
    for max in range(0,len(max_apmlitudes_u1)):
        max_apmlitudes_u1[max] = [source_velocities[max], max_apmlitudes_u1[max]]
        max_apmlitudes_u9[max] = [source_velocities[max], max_apmlitudes_u9[max]]
        max_apmlitudes_uavg[max] = [source_velocities[max], max_apmlitudes_uavg[max]]
    
    return  (max_apmlitudes_u1, max_apmlitudes_u9, max_apmlitudes_uavg)

def chanegResultAmplitudes(source, outputfilename, resultFilePath, dataFiles, iteration_number, source_velocities, index_range_for_local_maxima, cuts): 
    result_file =  resultFilePath + source + ".json"
    data_file = dataFiles + outputfilename + ".dat"
    
    with open(result_file) as result:
        result_json = json.load(result)
    
    max_apmlitudes_u1, max_apmlitudes_u9, max_apmlitudes_uavg = getLocalMax(data_file, source_velocities, index_range_for_local_maxima, cuts)
    
    for key in result_json.keys():
        if key.endswith(iteration_number):
            result_json[key]["polarizationU1"] =  max_apmlitudes_u1
            result_json[key]["polarizationU9"] = max_apmlitudes_u9
            result_json[key]["polarizationAVG"] = max_apmlitudes_uavg
            
    resultFile = open (resultFilePath + source + ".json", "w")
    resultFile.write(json.dumps(result_json, indent=2))
    resultFile.close()

def main():
    args = parseArguments()
    outputfilename = str(args.__dict__["outputfilename"])
    
    configFilePath = str(args.__dict__["config"])
    config = configparser.RawConfigParser()
    config.read(configFilePath)
    resultFilePath =  config.get('paths', "resultFilePath")
    dataFiles =  config.get('paths', "dataFilePath")
    source = outputfilename.split("_")[0]
    source_velocities = config.get('velocities', source).replace(" ", "").split(",")
    index_range_for_local_maxima = int(config.get('parameters', "index_range_for_local_maxima"))
    cuts = config.get('cuts', source).split(";")
    cuts = [c.split(",") for c in  cuts]
    iteration_number = outputfilename.split("_")[-1]
    chanegResultAmplitudes(source, outputfilename, resultFilePath, dataFiles, iteration_number, source_velocities, index_range_for_local_maxima, cuts)

if __name__=="__main__":
    main()
