# prelimaries for the benchmark files are:
# A benchmark consists of two lines, a description line and a value line
# The description line holds the identifiers for the benchmark seperated by commatas
# The value line holds the values in the same order as the description line for the benchmark seperated by commatas
# The last entry of the description is "Durations" which may have an arbitrary number of values


import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

import os
import sys
import time, datetime

# if set to false all existing benchmarks are ignored
# if set to true all benchmarks from the data directory are converted to plot diagrams
REBUILD_ALL = False
# input and output directories
DATA_DIR = r"../data"
OUTPUT_DIR = r"../output"

# get path of current execution
ABS_PATH = os.path.dirname(os.path.abspath(__file__))

# get the absolute path to input and output directory
DATA_PATH = os.path.join(ABS_PATH, DATA_DIR)
OUTPUT_PATH = os.path.join(ABS_PATH, OUTPUT_DIR)


"""
Entry point of plotter script
"""
def main():

    global REBUILD_ALL, DATA_PATH

    # get argument from command line
    if "-r" in sys.argv:
        REBUILD_ALL = True
    if REBUILD_ALL:
        print("Rebuild all benchmarks")

    # iterate through speedup plots files and process the files
    speedupPath = os.path.join(DATA_PATH, "speedup")
    for subdir, dirs, files in os.walk(speedupPath):
        for file in files: 
            processSpeedupFile(os.path.join(subdir, file));


"""
This function will take the location of a file with speedup values as an
argument and will convert the file content to a speedup plot diagram
"""
def processSpeedupFile(filePath):
    print("Process Speedup file: ") 
    print(filePath)

    (header, data) = parseBenchmarkFile(filePath)

    # set up seaborn plot
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(30,30))
    ax = fig.add_subplot(111)

    # set up labels
    ax.set_xlabel("Number of threads", fontsize = 20)
    ax.set_ylabel("Speedup", fontsize = 20)
    #ax.set_xscale("log", basex = 2)

    plt.title("Place holder", fontsize = 24)

    plt.scatter(x = 1, y = np.average(data["1"].asArray()))





"""
This funciton will take the location of a file with timing values as an
argument and will convert the file content to a speedup plot diagram
"""
def processTimingFile(filePath):
    pass

"""
This function will take the location of a file as an argument and will 
convert the file content to a plot diagram
"""
def processFile(filePath):

    # a list of all lines in the file
    lines = []    
    
    # open file and iterate over lines
    with open(filePath) as fp:
        for line in fp:
            lines.append(line)

    # a list of all benchmarks from this file
    benchmarks = []

    # iterator
    lineNo = 0
    while lineNo < len(lines) - 1:
        if len(lines[lineNo].rstrip()) > 0:
            benchmarks.append(parseBenchmark(lines[lineNo], lines[lineNo + 1]))
            lineNo = lineNo + 2

        else:
            lineNo = lineNo + 1

    for benchmark in benchmarks:
        plotBenchmark(benchmark)

"""
This function parses a benchmark file from the csv format to a pandas dataframe
@returns: This function returns a tuple of header and data of the parsed csv file
"""
def parseBenchmarkFile(benchmarkFilePath):

    header = pd.read_csv(benchmarkFilePath, nrows = 1)
    data = pd.read_csv(benchmarkFilePath, comment = '#', skiprows = 2)

    return (header, data)


    """
    # parse input to lists of entries
    descriptionList = descriptionLine.rstrip('\n').replace(" ", "").split(',')
    valueList = valueLine.rstrip('\n').replace(" ", "").split(',')

    # new dictionary
    benchmarkDict = {}

    # fill the dictionary with the values
    index = 0
    while (index < len(descriptionList) and descriptionList[index] != "Durations"):
        benchmarkDict[descriptionList[index]] = valueList[index]
        index = index + 1

    # special treatment for durations as there may be arbitrary many
    benchmarkDict[descriptionList[index]] = valueList[index:]

    return benchmarkDict    
    """

"""
This function takes a benchmark dictionary and plots the diagram
"""
def plotBenchmark(benchmarkDict):
    # create filename:
    dateStr = benchmarkDict['TimeStamp']
    dateStruct = time.strptime(dateStr, "%d-%m-%Y%H%M%S")
    dateStr = time.strftime("%Y%m%d_%H%M%S", dateStruct)

    fileName = benchmarkDict['Algorithm'] + '_' + benchmarkDict['GraphName'] + '_' + dateStr

    pass

#df = pd.DataFrame()

#df['x'] = random.sample(range(1, 100), 25)
#df['y'] = random.sample(range(1, 100), 25)

#df.head()

#sns.set_stype('whitegrid')

#sns_plot = sns.lmplot('x', 'y', data=df, fit_reg=False)
#sns_plot.savefig('output.png')    

#df = pd.DataFrame()
    
#df = sns.load_dataset("exercise")

#g = sns.factorplot(x="time", y="pulse", hue="kind", col="diet", data=df,
#                   capsize=.2, palette="YlGnBu_d", size=6, aspect=.75)


if __name__ == "__main__":
    main();
    print("Done")