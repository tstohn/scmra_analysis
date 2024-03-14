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
cellNumArray = [1000]
replications = 20

#helper list to select missing totatl proteins
rtot_names = [s for s in wt.columns if "Tot" in s]

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
#noiseArray = [0.5]
noiseArray = [0.2]

header = ['NOISE', 'PERT', 'MISS', 'ETA', 'REPEAT','CELLNUM','MISSNAMES', 'PERTNAMES', 'EDGENUM']
result = MRASimulationClass(header)

logFile = "../bin/logPertMissPapers2.txt"
f = open(logFile, "w+")

#simulations to run: to keep track of using e.g. same missing total when haviong perturbation fixed (first add 4 missing, then another 4 for 8 missing, ..)
pertMissSimulations = [{"P":[5,8,11], "M":[11]}, {"P":[0], "M":[0,3,6,11]}]

for cellNum in cellNumArray:
    for noise in noiseArray:
        for simulation in pertMissSimulations:
            for repeat in range(replications):

                pertValues = simulation["P"]
                missingNodesValues = simulation["M"]

                oldPertValue = 0
                oldMissValue = 0
                pertNodes = []
                missingNodes = []
                pertNodesSamplePool = allNodes
                missNodesSamplePool = allNodes
                for pertValue in pertValues:
                    pertValueDiff = pertValue - oldPertValue
                    pertNodes += random.sample(pertNodesSamplePool, pertValueDiff)
                    pertNodesSamplePool = [x for x in pertNodesSamplePool if x not in pertNodes]

                    for missingNodesVal in missingNodesValues:
                        missingNodeDiff = missingNodesVal - oldMissValue
                        missingNodes += random.sample(missNodesSamplePool, missingNodeDiff)
                        missNodesSamplePool = [x for x in missNodesSamplePool if x not in missingNodes]

                        #get names of missing and eprt
                        missingString = ""
                        perturbedString = ""
                        for missingNode in missingNodes:
                            missingString += ','
                            missingString += str(missingNode)
                        for perturbedNode in pertNodes:
                            perturbedString += ','
                            perturbedString += str(perturbedNode)

                        #for backwards compatibility of old scripts just convert None and mepty list into each other for different functions
                        if(len(pertNodes) == 0): pertNodes = None
                        if(len(missingNodes) == 0): missingNodes = None
                        pertubationIndices = get_index_list(allNodes, pertNodes)
                        rglob, rtot, cell_annot, tx_annot, pDict = combine_control_and_perturbed_cells(pertubationIndices, cellNum = cellNum, noise = noise)
                        if(missingNodes is not None):
                            rtot = rtot.drop(missingNodes, axis=0)
                        if(pertNodes == None): pertNodes = []
                        if(missingNodes == None): missingNodes = []

                        eta, scMraResult, edgeNumber = get_suitable_solution_fast(rglob, rtot, cell_annot, tx_annot, etaGuesses=30)

                        #add values to result
                        result.add_parameter([noise, pertValue, missingNodesVal, eta, repeat, cellNum, missingString, perturbedString, edgeNumber])

                        #add simulation output to result
                        result.add_result(scMraResult)

                        #update logfile
                        f.write("SIMULATION: noise=" + str(noise) + ", pert=" + str(pertValue) + ", missing=" + str(missingNodesVal) + ", repeat=" + str(repeat) + ", CellNum=" + str(cellNum) +",eta: ")
                        f.write(str(round(eta,3)) + " edges: " + str(edgeNumber) + "\n")
                        f.flush()
                        oldMissValue = missingNodesVal
                    oldPertValue = pertValue

#store pickled result
now = datetime.now()
folder = "./bin/"
name = "MISSINGTOTAL_PERTURBATION_2"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))
