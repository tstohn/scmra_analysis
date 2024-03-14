import sys
sys.path.insert(0, '../../helperScripts/')
from IDseqHelperFunctions import *

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
import scmra
import itertools
import umap
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
#distribution of distances in cell space
from sklearn.metrics import pairwise_distances

import random

'''
bootstrapping to check how often edges r detected
or dependant on outliers
'''

def bootstrap_data(data):

    vehicleData = (data[data["treatment"] == "vehicle"])
    vehicleData = vehicleData.set_index(['sample_id',"treatment"])
    vehicleData = vehicleData.pivot(columns='ab_name', values='ab_count_tmm')
    samplingVehicle = vehicleData.sample(frac=1.0, replace=True)
    proteinVars = samplingVehicle.columns
    samplingVehicle = samplingVehicle.reset_index()
    #change sample id bcs of duplications
    samplingVehicle['range'] = (np.arange(len(samplingVehicle))).astype(str)
    samplingVehicle["sample_id"] = samplingVehicle["sample_id"] + samplingVehicle["range"]
    samplingVehicle = samplingVehicle.drop(['range'], axis=1)
    #remelt to origional format
    samplingVehicle = pd.melt(samplingVehicle, id_vars=["sample_id", "treatment"], value_vars=proteinVars)


    treatedData = (data[data["treatment"] == "AG1478"])
    treatedData = treatedData.set_index(['sample_id',"treatment"])
    treatedData = treatedData.pivot(columns='ab_name', values='ab_count_tmm')
    samplingTreated = treatedData.sample(frac=1.0, replace=True)
    proteinVars = samplingTreated.columns
    samplingTreated = samplingTreated.reset_index()
    #change sample id bcs of duplications
    samplingTreated['range'] = np.arange(len(samplingTreated)).astype(str)
    samplingTreated["sample_id"] = samplingTreated["sample_id"] + samplingTreated["range"]
    samplingTreated = samplingTreated.drop(['range'], axis=1)
    #remelt to origional format
    samplingTreated = pd.melt(samplingTreated, id_vars=["sample_id", "treatment"], value_vars=proteinVars)

    frames = [samplingTreated, samplingVehicle]
    result = pd.concat(frames)
    result = result.rename(columns={"value": "ab_count_tmm"})

    return(result)

PERMUTATIONS = 1000
THETA = GLOBALTHETA_EGFRINH

#returns a list
def run_bootstrapping(egfrInhDataProcessed, masterCnrResult):

    cols = ["SIMULATION","GROUP","EDGE","VALUE"]
    results = pd.DataFrame(columns = cols)

    for i in range(PERMUTATIONS):

        Bdata = bootstrap_data(egfrInhDataProcessed)
        rglob, rtot, cell_annot, tx_annot, group_annot = prepare_data_for_EGFR_inhibition(Bdata)

        #SOLVE THE NETWORK
        scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)
        scp = scmra.ScCnrProblem(scd, eta=0.0, theta=THETA, prior_network = CANNONICAL_EDGES) 
        scp.cpx.solve()
        s = scmra.ScCnrResult(scp)

        #get all edges that deviate in this simulation
        deviatingFromToList = []
        for indicator, status in s.vardict.items():
            if(indicator.startswith("IDev") and status==1):
                indSplit = indicator.split("_")
                deviatingFromToList.append([indSplit[1],indSplit[2]])

        #make a matrix storing edge, weight, cluster
        dfEdgeDiff = pd.DataFrame(columns = cols)
        for edge in deviatingFromToList:
            edgeName = edge[0] + "_" + edge[1]
            edgeDiff = (s.rloc["vehicle"].loc[edge[0],edge[1]]) - (s.rloc["AG1478"].loc[edge[0],edge[1]])
            arr = [i, "AG1478", edgeName, edgeDiff]
            dfEdgeDiff = pd.concat([dfEdgeDiff, pd.DataFrame([arr], columns=["SIMULATION","GROUP","EDGE","VALUE"])], axis=0)
        results = pd.concat([results, dfEdgeDiff])

    return(results)

#write two dataframes: one with the true edge differences (between each population and the vehicle)
#and the calcualted edge differences from the label-permuted reconstructions
def write_edge_subpopulation_differences(trueEdges, resultPermSimulations):

    resultPermSimulations.to_csv("../bin/bootstrap_RecoveredDifferences_EGFR.tsv", sep="\t")

    #WRITE THE EDGE DIFFERENCES FORM MASTERCNR TO A TSV
    #get all different edges
    cols = ["GROUP","EDGE","VALUE"]
    colsToKeep = []
    for indicator, status in trueEdges.vardict.items():
        if(indicator.startswith(("IDev")) and status==1):
            indSplit = indicator.split("_")
            colsToKeep.append(indSplit[1] + "_" + indSplit[2])
    truePoint = pd.DataFrame(columns = cols)
    populations = list(trueEdges.rloc.keys())
    populations.remove("vehicle")
    for edge in colsToKeep:
        for pop in populations:
            toFromList = edge.split("_")
            val = trueEdges.rloc["vehicle"].loc[toFromList[0], toFromList[1]] - trueEdges.rloc[pop].loc[toFromList[0], toFromList[1]]
            arr = [pop, edge, val]
            truePoint = pd.concat([truePoint, pd.DataFrame([arr], columns=["GROUP","EDGE","VALUE"])], axis=0)

    truePoint.to_csv("../bin/bootstrap_TrueDifferences_EGFR.tsv", sep="\t")


egfrInhDataProcessed = pd.read_csv("../bin/EGFRInhibitionDataProcessed.tsv",sep='\t')



objectMasterResult = open("../bin/MASTERCNR_EGFRINH_CNRRESULT.pickle", "rb")
masterCnrResult = pickle.load(objectMasterResult)

resultPermSimulations = run_bootstrapping(egfrInhDataProcessed, masterCnrResult)
write_edge_subpopulation_differences(masterCnrResult, resultPermSimulations)
