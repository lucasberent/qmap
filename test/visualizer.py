import csv

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plotEC():
    filenameEC = '/home/luca/Desktop/benchmarkEC 17-02-2022.json'
    # fin = open(filenameEC, "rt")
    # data = fin.read()
    # data = data.replace(',', '', 1)
    # fin.close()
    # fin = open(filenameEC, "wt")
    # fin.write(data)
    # fin.close()

    fig, ax = plt.subplots()
    ax.axis('tight')
    ax.axis('off')

    df = pd.read_json(filenameEC)['benchmarks']
    columns = (
        'Qubits', 'Sum Gates', 'Sum SatVars', 'Unique Gens', 'Input States', 'Max Depth', 'Struc An (ms)',
        'Sum constr(ms)',
        'Equivalent', '#Conflicts', 'Solv(ms)')

    nrRows = len(df)
    cellText = []

    for row in range(nrRows):
        cellText.append(
            [df[row]['nrOfQubits'], df[row]['numGates'], df[row]['numSatVarsCreated'], df[row]['numGenerators'],
             df[row]['numInputStates'], df[row]['circDepth'],
             df[row]['preprocTime'], df[row]['satConstructionTime'], df[row]['equivalent'],
             df[row]['z3map']['sat conflicts'] if 'sat conflicts' in df[row]['z3map'] else '-',
             df[row]['solvingTime']])

    table = ax.table(cellText=cellText, loc='center', colLabels=columns)
    table.auto_set_font_size(False)
    table.set_fontsize(14)

    outfile = open('/home/luca/Desktop/ec_out.csv', 'w')
    writer = csv.writer(outfile)
    writer.writerow(columns)

    for row in range(nrRows):
        writer.writerow(
            [df[row]['nrOfQubits'], df[row]['numGates'], df[row]['numSatVarsCreated'], df[row]['numGenerators'],
             df[row]['numInputStates'], df[row]['circDepth'],
             df[row]['preprocTime'], df[row]['satConstructionTime'], df[row]['equivalent'],
             df[row]['z3map']['sat conflicts'] if 'sat conflicts' in df[row]['z3map'] else '-',
             df[row]['solvingTime']])
    plt.show()


def plotQB():
    filenameQB = '/home/luca/Desktop/benchmarkQB 17-02-2022.json'
    # fin = open(filenameQB, "rt")
    # data = fin.read()
    # data = data.replace(',', '', 1)
    # fin.close()
    # fin = open(filenameQB, "wt")
    # fin.write(data)
    # fin.close()

    fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
    df = pd.read_json(filenameQB)['benchmarksQubit']
    data = []
    data2 = []
    x = [z for z in range(1, 128, 1)]

    for row in range(len(df)):
        data.append(df[row]['preprocTime'])
        data2.append(df[row]['satConstructionTime'])

    ax.plot(x, data, label='preprocessing time')
    ax.plot(x, data2, label='sat construction time')
    ax.legend()
    plt.xlabel('Number of Qubits')
    plt.ylabel('Time(ms)')
    plt.title("SAT construction time per #qubits")
    plt.show()


def plotCS():
    filenameCS = '/home/luca/Desktop/benchmarkCS 17-02-2022.json'
    # fin = open(filenameCS, "rt")
    # data = fin.read()
    # data = data.replace(',', '', 1)
    # fin.close()
    # fin = open(filenameCS, "wt")
    # fin.write(data)
    # fin.close()

    fig, ax = plt.subplots(figsize=(5, 2.7), layout='constrained')
    df = pd.read_json(filenameCS)['benchmarksCS']
    data = []
    data2 = []
    x = []

    for i in range(len(df)):
        x.append(df[i]['numGates'])

    for row in range(len(df)):
        data.append(df[row]['preprocTime'])
        data2.append(df[row]['satConstructionTime'])

    ax.plot(x, data, label='preprocessing time')
    ax.plot(x, data2, label='sat construction time')
    plt.xlabel('Number of Gates')
    plt.ylabel('Time(ms)')
    plt.title("SAT construction time per #gates, 10 qubits")
    ax.legend()
    plt.show()


# plotQB()
# plotCS()
plotEC()
