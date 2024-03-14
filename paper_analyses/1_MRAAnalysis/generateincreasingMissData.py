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
replications = 20

#helper list to select missing totatl proteins
rtot_names = [s for s in wt.columns if "Tot" in s]

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
#noiseArray = [0.5]
noiseArray = [0.2]
pertValues = [0]
#pertValues = [0]

missingNodesValues = [0,1,2,3,4,5,6,7,8,9,10,11]
#missingNodesValues = [11]

header = ['NOISE', 'PERT', 'MISS', 'ETA', 'REPEAT','CELLNUM']
result = MRASimulationClass(header)

logFile = "../bin/logincreasingMiss1609_2.txt"
f = open(logFile, "w+")

for cellNum in cellNumArray:
    for noise in noiseArray:
        for pertValue in pertValues:
            for missingNodesVal in missingNodesValues:
                count = 0
                for repeat in range(replications):
                    count += 1

                    pertNodes = random.sample(allNodes, pertValue)
                    missingNodes = random.sample(allNodes, missingNodesVal)
                    if(len(pertNodes) == 0): pertNodes = None
                    if(len(missingNodes) == 0): missingNodes = None

                    pertubationIndices = get_index_list(allNodes, pertNodes)
                    rglob, rtot, cell_annot, tx_annot, pDict = combine_control_and_perturbed_cells(pertubationIndices, cellNum = cellNum, noise = noise)
                    if(missingNodes is not None):
                        rtot = rtot.drop(missingNodes, axis=0)
                    etasInput = [0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
                    eta, etaListNew = estimate_eta(rglob, rtot, cell_annot, tx_annot, withFlattening = True, etasInput=etasInput)

                    if(eta <=0): break 

                    #add values to result
                    result.add_parameter([noise, pertValue, missingNodesVal, eta, count, cellNum])

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
                    f.write("SIMULATION: noise=" + str(noise) + ", pert=" + str(pertValue) + ", missing=" + str(missingNodesVal) + ", repeat=" + str(count) + ", CellNum=" + str(cellNum) +",eta: ")
                    f.write(str(round(eta,3)) + " edges: " + str(numEdges) + "\n")
                    f.flush()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_OF_increasing_MISS_1609_2"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))
