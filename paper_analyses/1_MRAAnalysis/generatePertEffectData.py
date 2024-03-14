'''   THIS SCRIPT IS TO SHOW THE EFFECT OF ADDITIONAL PERTURBATIONS TO NETWORK RECONSTRUCTION
        IT WAS MADE AFTER THE ADDITIONAL OFFSET TERM FOR THE LIENAR ACTIVATION OF PERTURBED DATA SCALED TO CTR
'''

#generate data do compare network quality
import sys
sys.path.insert(0, '../helperScripts/')
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

#suppress all output except error messages
sys.stdout = open(os.devnull, "w")
sys.stderr = open(os.devnull, "w")
sys.stdout = sys.__stdout__
sys.stderr = sys.__stderr__

#GLOBAL SIMULATION PARAMETERS
cellNum =10000
replications = 1

#helper list to select missing totatl proteins
rtot_names = [s for s in wt.columns if "Tot" in s]

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
#noiseArray = [0, 0.5]
noiseArray = [0.02]
pertValues = [11]

#missingNodesValues = [0,5,11]
missingNodesValues = [11]

header = ['NOISE', 'PERT', 'MISS', 'ETA', 'REPEAT']
result = MRASimulationClass(header)

logFile = "../bin/logManyCells.txt"
f = open(logFile, "a")

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
                rglob, rtot, cell_annot, tx_annot, pDict = combine_control_and_perturbed_cells(pertubationIndices, cellNum = cellNum, noise = noise, scaleByTxGroupIndividually = True)
                if(missingNodes is not None):
                    rtot = rtot.drop(missingNodes, axis=0)
                #eta, etaListNew = estimate_eta(rglob, rtot, cell_annot, tx_annot, withFlattening=True, modelPertRloc = True, modelPertStot=True, perturbedNodeDict = pDict,logFile="/DATA/t.stohn/MRA/scMRA-analysis/notebooks/NetworkAnalysis/logFile.txt",
                #                   etasInput= [0.01,0.04,0.08,0.1,0.14,0.18,0.2,0.24,0.28,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.63,0.66,0.7,0.73,0.76,0.8,0.9])
                eta = 0.91

                if(eta <=0): break 
                #update logfile
                f.write("SIMULATION: noise=" + str(noise) + ", pert=" + str(pertValue) + ", missing=" + str(missingNodesVal) + ", repeat=" + str(count) + " eta: " + str(round(eta,3)))
                f.flush()
                #add values to result
                result.add_parameter([noise, pertValue, missingNodesVal, eta, count])

                #make MRA simulation
                scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot)
                p = scmra.ScMraProblem(scd, eta=eta, modelPertRloc = True, modelPertStot=True, perturbedNodeDict = pDict) #0.3765 1pert   #0.5566 11pert #zero 0.02775
                p.cpx.solve()
                s = scmra.ScMraResult(p)

                rmse = calc_rmse(rloc_true, s.rloc)
                f.write(" RMSE: " + str(rmse) + "\n")
                f.flush()

                #add simulation output to result
                result.add_result(s)
                del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)


#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_MANY_CELLS"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))
