import matplotlib as mpl
mpl.use("Agg")
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from os import listdir, remove
from os.path import isfile, join, exists

# Constants
REBUILD_DATA = False
COLD_START = 1
DATA_DIR = "../data/speedup/coPaper/"
OUT_DIR = "../output/speedup/coPaper/"
PLOT_DATA = "plotData.csv"
BASELINE = 0.11525239

def main():
    # If REBUILD_DATA is True, delete old data
    # an build a new file
    if REBUILD_DATA:
        f = join(DATA_DIR, PLOT_DATA)
        if exists(f):
            remove(f)

        buildPlotData(DATA_DIR)

    plotSpeedUp(DATA_DIR + PLOT_DATA)

def plotSpeedUp(dataPath):
    data = pd.read_csv(dataPath, index_col=0)

    # Convert durations to speedup value
    # Comment this out for timing plot
    data["Dur"] = BASELINE / data["Dur"].values

    # Select data to plot
    # Syntax: plotData = data[(data.Cloumn==Value)]
    # It's also possible to combine multiple conditions
    # with & (and) and | (or):
    # plotData = data[(data.Algo=="ppf1") | (data.Algo=="ppf2")]
    plotData = data[(data.Algo=="ppf1_KS")]

    # Plot style
    sns.set_style("whitegrid")
    sns.set_context("talk")

    # Plot
    v = sns.factorplot(x="NThread", y="Dur", hue="Algo",
             data=plotData, size=10, aspect=2, kind="violin")

    # Graphical options
    v.set(xticks=np.arange(8, 264, 8))
    v.set(xticklabels=np.arange(8, 264, 8))
    v.ax.yaxis.set_major_formatter(FormatStrFormatter("%d"))
    v.ax.set_ylim(ymin=0)
    v.ax.set_xlabel("Number of threads", fontsize=20)
    v.ax.set_ylabel("Speedup", fontsize=20)
    plt.title("Place holder", fontsize=24)

    # Save plot
    outFileName = OUT_DIR + "coPaper_pf_KS_test.png"
    v.savefig(filename=outFileName)
    print("Saved plot to: " + outFileName)

def buildPlotData(dirPath):
    plotData = pd.DataFrame(columns=["TS", "Algo", "Graph", "NThread", "Dur"])
    dataFrames = []
    fileNames = [f for f in listdir(dirPath) if isfile(join(dirPath,f))]

    for f in fileNames:
        print("Processing file: " + f)
        dataFrames.append(loadTransform(join(dirPath,f)))

    plotData = pd.concat(dataFrames, ignore_index=True)
    plotData.to_csv(dirPath + PLOT_DATA)
    print("Writing plot data to: " + dirPath + PLOT_DATA)


def loadTransform(filePath):
    header = pd.read_csv(filePath, comment='#', nrows=1).dropna(axis=1)
    data = pd.read_csv(filePath, comment='#', skiprows=2).dropna(axis=1)

    newHeader = ["TS", "Algo", "Graph", "NThread", "Dur"]
    newData = pd.DataFrame(columns=newHeader)

    ts = header["TimeStamp"].values[0]
    algo = header["Algorithm"].values[0]
    graph = header["GraphName"].values[0]

    tempDataFrames = []

    for col in range(0, data.shape[1]):
        # Iterate over columns
        nthread = data.columns.values[col]
        durations = data[nthread]

        temp = []
        for row in range(0, len(durations)):
            # Ignore first COLD_START items
            if row < COLD_START:
                continue

            temp.append([ts, algo, graph, int(nthread), float(durations[row])])

        tempDataFrames.append(pd.DataFrame(data=temp, columns=newHeader))

    newData = pd.concat(tempDataFrames, ignore_index=True)
    return newData


if __name__ == "__main__":
    main();
    print("Done")
