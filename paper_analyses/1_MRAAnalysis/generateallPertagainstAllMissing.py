#generate data do compare network quality
import sys
sys.path.insert(0, '../../helperScripts/')
from helperFunctions import *
import scmra
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import math
import pandas as pd
from datetime import datetime
import pickle
import seaborn as sns; sns.set_theme()
from matplotlib.pyplot import figure

#GLOBAL SIMULATION PARAMETERS
cellNumArray = [500]
replications = 10

#helper list to select missing totatl proteins
rtot_names = [s for s in wt.columns if "Tot" in s]

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
#noiseArray = [0.5]
noiseArray = [0.2]

header = ['NOISE', 'PERT', 'MISS', 'ETA', 'REPEAT','CELLNUM']
result = MRASimulationClass(header)

logFile = "../bin/logPertMiss0710_allAgainstAll.txt"
f = open(logFile, "w+")

for cellNum in cellNumArray:
    for noise in noiseArray:
        for pertName in allNodes:
            for missName in allNodes:
                count = 0
                for repeat in range(replications):
                    count += 1

                    pertNodes = [pertName]
                    missingNodes = [missName]

                    pertubationIndices = get_index_list(allNodes, pertNodes)
                    rglob, rtot, cell_annot, tx_annot, pDict = combine_control_and_perturbed_cells(pertubationIndices, cellNum = cellNum, noise = noise)
                    if(missingNodes is not None):
                        rtot = rtot.drop(missingNodes, axis=0)
                    #etasInput = [0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.7, 0.9] for simulations woth 500 cells
                    etasInput = [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

                    eta, etaListNew = estimate_eta(rglob, rtot, cell_annot, tx_annot, withFlattening = True, etasInput=etasInput)

                    if(eta <=0): break 

                    #add values to result
                    result.add_parameter([noise, pertName, missName, eta, count, cellNum])

                    #make MRA simulation
                    scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot)
                    p = scmra.ScMraProblem(scd, eta=eta)
                    p.cpx.solve()
                    s = scmra.ScMraResult(p)

                    #add simulation output to result
                    result.add_result(s)
                    del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)

                    #add values to result
                    #GET NUMBER EDGES IN THIS NETWORK
                    numEdges = 0
                    exampleRloc = s.rloc
                    for column_name in exampleRloc:
                        column = exampleRloc[column_name]
                        numEdges += (column != 0).sum()
                        numEdges -=1 #for the diagonal values

                    #update logfile
                    f.write("SIMULATION: noise=" + str(noise) + ", pert=" + str(pertName) + ", missing=" + str(missName) + ", repeat=" + str(count) + ", CellNum=" + str(cellNum) +",eta: ")
                    f.write(str(round(eta,3)) + " edges: " + str(numEdges) + "\n")
                    f.flush()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_OF_PERT_MISS_0710_allAgainstAll"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))
