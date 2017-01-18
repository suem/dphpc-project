import matplotlib as mpl
mpl.use("Agg")
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.legend_handler import HandlerLine2D
from os import listdir, remove
from os.path import isfile, join, exists


# Constants

# GRAPH = "amazon"
#GRAPH = "coPaperDBLP"
GRAPH = "wikipedia"
# SPEEDUP_PLOT = True
SPEEDUP_PLOT = True

REBUILD_DATA = True
COLD_START = 1
DATA_DIR = "../data/speedup/" + GRAPH + "/"
OUT_DIR = "../output/speedup/" + GRAPH + "/"
PLOT_DATA = "plotData" + GRAPH + ".csv"

def main():
    # If REBUILD_DATA is True, delete old data
    # an build a new file
    if REBUILD_DATA:
        f = join(DATA_DIR, PLOT_DATA)
        if exists(f):
            remove(f)

        buildPlotData(DATA_DIR)

#   plotSpeedUp(DATA_DIR + PLOT_DATA)
    plotBoxplot(DATA_DIR + PLOT_DATA)

def plotBoxplot(dataPath):
    data = pd.read_csv(dataPath, index_col=0)

    # Speedup
    baselineGreedy = data[(data.Algo=='pf_Greedy')].Dur.median()
    data.loc[data['Algo']=='ppfTTAS_Greedy', 'Dur'] = baselineGreedy / data[
            (data.Algo== 'ppfTTAS_Greedy')].Dur.get_values()
    baselineKS = data[(data.Algo=='pf_KS')].Dur.median()
    data.loc[data['Algo']=='ppfTTAS_KS', 'Dur'] = baselineKS / data[
            (data.Algo=='ppfTTAS_KS')].Dur.get_values()

    # Select data to plot
    plotData = data[(data.Algo=='ppfTTAS_Greedy') |
                    (data.Algo=='ppfTTAS_KS')]

    # Plot style
    sns.set_style("whitegrid")
    sns.set_context("talk")
    sns.set_style('whitegrid', {'legend.frameon': True})
    fig = plt.figure(figsize=(20,15))
    plt.title('CopaperDBLP graph', fontsize=32)


    # Plot boxplot
    b = sns.boxplot(x='NThread', y='Dur', hue='Algo', data=plotData,
            notch=True, bootstrap=5000, width=0.8)

    # Plot medians
    ksData = plotData[(plotData.Algo=='ppfTTAS_KS')]
    x_data = np.arange(0, len(ksData.NThread.unique())) + 0.195
    y_data = ksData.groupby('NThread').Dur.median().get_values()
    b.plot(x_data, y_data, c='g', marker='o', ls='')

    grData = plotData[(plotData.Algo=='ppfTTAS_Greedy')]
    x_data = np.arange(0, len(grData.NThread.unique())) - 0.195
    y_data = grData.groupby('NThread').Dur.median().get_values()
    b.plot(x_data, y_data, c='b', marker='o', ls='')

    # Plot baseline
    # b.axhline(y=baselineKS, label='pf_KS', ls='dotted', c='c')
    # b.axhline(y=baselineGreedy, label='pf_Greedy', ls='dotted', c='r')

    # Figure aesthetics
    b.set_ylim(bottom=0)
    # b.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    b.set_xlabel('Number of threads', fontsize=22)
    b.set_ylabel('Speedup', fontsize=22)
    b.legend(title='Algorithm', prop={'size':24}, loc='best')
    b.tick_params(labelsize=22)

    plt.savefig('../output/test/' + GRAPH + '_GREEDYvsKS.png')
    plt.close(fig)

def plotSpeedUp(dataPath):
    data = pd.read_csv(dataPath, index_col=0)


    # Convert durations to speedup value
    # Comment this out for timing plot

    baseline_ks = data[data['Algo'] == 'pf_KS']['Dur'].median()
    print(baseline_ks)
    baseline_greedy = data[data['Algo'] == 'pf_Greedy']['Dur'].median()
    print(baseline_greedy)
    if SPEEDUP_PLOT:
        for index, row in data.iterrows():
            if row['Algo'].endswith('KS'):
                data.set_value(index, 'Dur', baseline_ks / data['Dur'][index])
            if row['Algo'].endswith('Greedy'):
                data.set_value(index, 'Dur', baseline_greedy / data['Dur'][index])

        # for i in range(10,250,10):
        #     line = {"TS": 0, "Algo": "lin", "Graph": GRAPH, "NThread": i, "Dur": i}
        #     print(line)
        #     data = data.append(line, ignore_index=True)


    # Select data to plot
    # Syntax: plotData = data[(data.Cloumn==Value)]
    # It's also possible to combine multiple conditions
    # with & (and) and | (or):
    # plotData = data[(data.Algo=="ppf1") | (data.Algo=="ppf2")]
    #plotData = data[(data.Algo=="ppf3_KS")]


    # plotData = data[(data.Algo=="ppfTAS_Greedy") | (data.Algo=="ppfTTAS_Greedy")]
    # title = "TAS vs TTAS (coPaperDBLP)"
    # outFileName = OUT_DIR + GRAPH + "TASvsTTAS.png"


    # plotData = data[(data.Algo=="pf_Greedy_i7") | (data.Algo=="pf_Greedy") | (data.Algo=="ppfTTAS_Greedy")]
    # title = "Runtimes (coPaperDBLP)"
    # outFileName = OUT_DIR + GRAPH + "runtimes.png"


    # plotData = data[(data.Algo=="ppfTTAS_Greedy")]
    # title = "Speedup (coPaperDBLP)"
    # outFileName = OUT_DIR + GRAPH + "speedup.png"


    # plotData = data[(data.Algo=="ptg_Greedy") | (data.Algo=="ppfTTAS_Greedy")]
    # title = "Speedup (coPaperDBLP)"
    # outFileName = OUT_DIR + GRAPH + "speedup_tg.png"


    plotData = data[(data.Algo=="ppfTTAS_KS") | (data.Algo=="ppfTTAS_Greedy")]
    title = "Speedup (wikipedia)"
    outFileName = OUT_DIR + GRAPH + "speedup_matching.png"

    # plotData = data[(data.Algo=="ppfTTAS_KS") | (data.Algo=="ppfTTAS_Greedy")]
    # title = "Speedup (coPaperDBLP)"
    # outFileName = OUT_DIR + GRAPH + "speedup_matching.png"

    # plotData = data[(data.Algo=="ppfTTAS_KS")]
    # title = "Speedup (amazon)"
    # outFileName = OUT_DIR + GRAPH + "speedupg.png"


    # plotData = data[(data.Algo=="ppfTTAS_KS")]
    # title = "Speedup (amazon)"
    # outFileName = OUT_DIR + GRAPH + "speedupg.png"


    # plotData = data[(data.Algo=="ppfTTAS_KS")]
    # title = "Runtimes (amazon)"
    # outFileName = OUT_DIR + GRAPH + "runtime.png"

    # plotData = data[(data.Algo=="ppfTTAS_KS") | (data.Algo=="ppfTTAS_Greedy")]
    # title = "Runtimes (coPaperDBLP)"
    # outFileName = OUT_DIR + GRAPH + "runtimes_matching.png"

    # plotData = data[(data.Algo=="ppfTTAS_KS") | (data.Algo=="ppfTTAS_Greedy")]
    # title = "Runtimes (wikipedia)"
    # outFileName = OUT_DIR + GRAPH + "runtimes_matching.png"


    # Plot style
    sns.set_style("whitegrid")
    sns.set_context("talk")
    sns.set_style('whitegrid', {'legend.frameon': True})

    # Plot
    # kind = "violin"
    kind = "point"
    # kind = "box"
    v = sns.factorplot(x="NThread", y="Dur", hue="Algo", data=plotData, size=10, aspect=2, kind=kind, legend_out=False)
    plt.legend(loc='upper right', fontsize='20')
    if SPEEDUP_PLOT:
        plt.plot(np.linspace(0,24), np.linspace(10,250), 'k')

    # v.despine(left=True)

    # Graphical options
    #v.set(xticks=np.arange(8, 264, 8))
    #v.set(xticklabels=np.arange(8, 264, 8))
    formatStr = "%.2f"
    label = "Runtime [s]"
    if SPEEDUP_PLOT:
        formatStr = "%d"
        label = "Speedup"
        v.ax.set_ylim(ymin=0, ymax=360)
        # v.set(yticks=np.arange(0, 260, 10))
    else:
        v.ax.set_ylim(ymin=0)


    v.ax.yaxis.set_major_formatter(FormatStrFormatter(formatStr))
    v.ax.set_xlabel("Number of threads", fontsize=24)
    v.ax.set_ylabel(label, fontsize=24)
    plt.title(title, fontsize=32)

    # Save plot
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
