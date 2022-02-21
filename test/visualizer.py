import csv
import math
import statistics
import numpy as np
import scipy.optimize as opt

import matplotlib.pyplot as plt
import pandas as pd


# equivalence checking, eq and not eq instances separate
# nrqubits, nr gates, solver time, conflicts analysis+construction
def plotEC():
    filenameEC = '/home/luca/Desktop/benchmarkEC 19-02-2022.json'

    fig, ax = plt.subplots()
    ax.axis('tight')
    ax.axis('off')

    df = pd.read_json(filenameEC)['benchmarks']
    columns = (
    'Qubits', 'Gates', 'Construction', 'Solving', 'Clauses', 'Equivalent', 'Qubits', 'Gates', 'Construction', 'Solving',
    'Clauses', 'Equivalent')
    nrRows = len(df)
    eqIdx = math.floor(nrRows / 2)
    rows = []

    for row in range(eqIdx):
        rows.append(
            [df[row]['nrOfQubits'],
             df[row]['numGates'],
             df[row]['preprocTime'] + df[row]['satConstructionTime'],
             df[row]['solvingTime'],
             # df[row]['z3map']['sat conflicts'] if 'sat conflicts' in df[row]['z3map'] else '-',
             int(df[row]['z3map']['sat mk clause 2ary'] if 'sat mk clause 2ary' in df[row]['z3map'] else 0 +
                                                                                                         df[row][
                                                                                                             'z3map'][
                                                                                                             'sat mk clause 3ary'] if 'sat mk clause 3ary' in
                                                                                                                                      df[
                                                                                                                                          row][
                                                                                                                                          'z3map'] else 0 +
                                                                                                                                                        df[
                                                                                                                                                            row][
                                                                                                                                                            'z3map'][
                                                                                                                                                            'sat mk clause nary'] if 'sat mk clause nary' in
                                                                                                                                                                                     df[
                                                                                                                                                                                         row][
                                                                                                                                                                                         'z3map'] else 0),
             df[row]['equivalent']
             ])

    ueqIdx = len(df) - eqIdx
    rows2 = []
    for row in range(ueqIdx, len(df), 1):
        rows2.append([
            df[row]['nrOfQubits'],
            df[row]['numGates'],
            df[row]['preprocTime'] + df[row]['satConstructionTime'],
            df[row]['solvingTime'],
            # df[row]['z3map']['sat conflicts'] if 'sat conflicts' in df[row]['z3map'] else '-',
            int(df[row]['z3map']['sat mk clause 2ary'] if 'sat mk clause 2ary' in df[row]['z3map'] else 0 +
                                                                                                        df[row][
                                                                                                            'z3map'][
                                                                                                            'sat mk clause 3ary'] if 'sat mk clause 3ary' in
                                                                                                                                     df[
                                                                                                                                         row][
                                                                                                                                         'z3map'] else 0 +
                                                                                                                                                       df[
                                                                                                                                                           row][
                                                                                                                                                           'z3map'][
                                                                                                                                                           'sat mk clause nary'] if 'sat mk clause nary' in
                                                                                                                                                                                    df[
                                                                                                                                                                                        row][
                                                                                                                                                                                        'z3map'] else 0),
            df[row]['equivalent']
        ])
    outfile = open('/home/luca/Desktop/ec-out.csv', 'w')
    writer = csv.writer(outfile)
    writer.writerow(columns)
    for i in range(len(rows)):
        writer.writerow(rows[i] + rows2[i])


def logFunc(x, a, b):
    return a * np.log2(x) + b


# preprocessing+satconstruction time in nr of qubits
def plotScaling():
    filenameQB = '/home/luca/Desktop/benchmarkQB 19-02-2022.json'
    filenameCS = '/home/luca/Desktop/benchmarkCS 19-02-2022.json'

    fig, ax = plt.subplots(figsize=(10, 5), layout='constrained', nrows=2, ncols=2)
    dfQB = pd.read_json(filenameQB)['benchmarks']
    dfCS = pd.read_json(filenameCS)['benchmarks']
    data = []
    data2 = []
    xDataQB = []
    dataCS = []
    data2CS = []
    xDataCS = []

    for i in range(0, len(dfQB), 10):
        xDataQB.append(dfQB[i]['nrOfQubits'])

    for i in range(0, len(dfCS), 10):
        xDataCS.append(dfCS[i]['numGates'])

    for row in range(0, len(dfQB), 10):
        tmp = []
        tmp2 = []
        for i in range(0, 10, 1):
            tmp.append(dfQB[row + i]['preprocTime'] + dfQB[row + i]['satConstructionTime'])
            tmp2.append(dfQB[row + i]['numGenerators'])
        data.append(statistics.mean(tmp))
        data2.append(statistics.mean(tmp2))

    for row in range(0, len(dfCS), 10):
        tmp = []
        tmp2 = []
        for i in range(0, 10, 1):
            tmp.append(dfCS[row + i]['preprocTime'] + dfCS[row + i]['satConstructionTime'])

            twocls = dfCS[row + i]['z3map']['sat mk clause 2ary'] if 'sat mk clause 2ary' in dfCS[row + i][
                'z3map'] else 0
            threecls = dfCS[row + i]['z3map']['sat mk clause 3ary'] if 'sat mk clause 3ary' in dfCS[row + i][
                'z3map'] else 0
            ncls = dfCS[row + i]['z3map']['sat mk clause nary'] if 'sat mk clause nary' in dfCS[row + i]['z3map'] else 0
            tmp2.append(twocls + threecls + ncls)
        dataCS.append(statistics.mean(tmp))
        data2CS.append(statistics.mean(tmp2))

    ax[0, 0].plot(xDataQB, data, label='SAT construction time')
    ax[0, 0].plot(xDataQB, xDataQB, label='linear')
    ax[0, 0].legend()
    ax[0, 0].set(xlabel='Number of Qubits')
    ax[0, 0].set(ylabel='Time(ms)')

    ax[0, 1].plot(xDataQB, data2, label='Clauses')
    optimizedParameters, pcov = opt.curve_fit(logFunc, xDataQB, data2)
    func = str(math.ceil(optimizedParameters[1])) + 'log(x)+' + str(math.ceil(optimizedParameters[0]))
    ax[0, 1].plot(xDataQB, logFunc(xDataQB, *optimizedParameters), label=func)
    ax[0, 1].legend()
    ax[0, 1].set(xlabel='Number of Qubits')
    ax[0, 1].set(ylabel='Number of Clauses')

    order = np.argsort(xDataCS)
    xDataCS = np.array(xDataCS)[order]
    dataCS = np.array(dataCS)[order]
    data2CS = np.array(data2CS)[order]

    ax[1, 0].plot(xDataCS, dataCS, 'o', label='SAT construction time, 10 Qubits', markevery=11, markersize=1)
    ax[1, 0].plot(xDataCS, np.sqrt(xDataCS), label='sqrt(x)')
    ax[1, 0].legend()
    ax[1, 0].set(xlabel='Number of Gates')
    ax[1, 0].set(ylabel='Time(ms)')

    ax[1, 1].plot(xDataCS, data2CS, 'o', label='#Clauses, 10 Qubits', markevery=5, markersize=1)
    ax[1, 1].plot(xDataCS, xDataCS, label='linear')
    ax[1, 1].legend()
    ax[1, 1].set(xlabel='Number of Gates')
    ax[1, 1].set(ylabel='Number of Clauses')

    plt.show()


def linfunc(x, k, a):
    return x * k + a


# preprocessing+satconstruction time in nr of qubits
def plotGenerators():
    filename = '/home/luca/Desktop/benchmarkCS 19-02-2022.json'
    xData = []
    yData = []
    xDataQB = []
    yDataQB = []
    df = pd.read_json(filename)['benchmarks']
    fig, ax = plt.subplots(figsize=(10, 5), layout='constrained')

    for i in range(0, len(df), 10):
        xData.append(df[i]['numGates'])

    for row in range(0, len(df), 10):
        for i in range(0, 10, 1):
            yData.append(df[row]['numGenerators'])

    order = np.argsort(xData)
    xData = np.array(xData)[order]
    yData = np.array(yData)[order]
    ax.plot(xData, yData, 'o', label='Number of Generators', markersize=2, markevery=5)
    optimizedParameters, pcov = opt.curve_fit(linfunc, xData, yData)
    func = str(math.ceil(optimizedParameters[1])) + 'x+' + str(math.ceil(optimizedParameters[0]))
    ax.plot(xData, linfunc(xData, *optimizedParameters), label=func)
    ax.legend()
    # ax.set_yscale('log')
    plt.xlabel('Number of Gates')
    plt.ylabel('Number of Generators')
    plt.title("Number of generators")
    plt.show()


# plotScaling()
# plotGenerators()
plotEC()
