#ONLY TO RUN PERMUTATIONS AND STORE THE ERROR
#NOT TO ACTAULLY GET EDGE STRENGTH, THEIR DIFFERENCES AND COMPARE TO PERMUTATION DATA

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

##############################
# GENERATE PERMUTATIONS DATA
##############################
    #permute randomly group annotation

    #only generates the pickle file: for edge significance run otherscript

import random

def shuffle_group_label(group_annot):
    keys = list(group_annot.keys())
    keyValueNumber = {}
    valueList= []
    for key in keys:
        keyValueNumber[key] = len(group_annot[key])
        valueList.append(group_annot[key])

    valueList = [item for sublist in valueList for item in sublist]

    group_annot_return = {}
    assert(len(keys) == 2)
    group_annot_return[keys[0]] = random.sample(valueList, keyValueNumber[keys[0]])
    group_annot_return[keys[1]] = list(set(valueList) - set(group_annot_return[keys[0]]))

    return(group_annot_return)

class Result:

    trueError = None
    randomErrorList = None
    theta = None

#RANDOM SHUFFLED LABELS
result = Result()
THETA = GLOBALTHETA_EGFRINH
results = []
errorList = []

data = pd.read_csv("../bin/EGFRInhibitionDataProcessed.tsv",sep='\t')

for i in range(1000):
    rglob, rtot, cell_annot, tx_annot, group_annot = prepare_data_for_EGFR_inhibition(data)
    group_annot = shuffle_group_label_3labels(group_annot)

    #scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)
    scd = scmra.ScData(rglob=rglob, rtot = rtot, group_annot= group_annot)

    p = scmra.ScCnrProblem(scd, eta=0.0, theta=THETA, prior_network = CANNONICAL_EDGES) 

    #fix edges of true network
    objectMasterResult = open("../bin/MASTERCNR_EGFRINH_CNRRESULT.pickle", "rb")
    masterCnrResult = pickle.load(objectMasterResult)
    for indicator, status in masterCnrResult.vardict.items():
        if indicator.startswith(("IDev", "ISDev", "IPDev")):
            #we need to also set contraints to ENFORCE this edge in networks reconstruction
            scmra.cplexutils.set_indicator_status(p.cpx, indicator, status)

    p.cpx.solve()
    s = scmra.ScCnrResult(p)
    del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)
    results.append(s)

    n_res = np.size(s.residuals_complete) + np.size(s.residuals_incomplete)
    error = 0
    error += np.sum(np.array(np.square(s.residuals_complete)))/n_res
    error += np.sum(np.array(np.square(s.residuals_incomplete)))/n_res
    errorList.append(error)

#ONE TRUE EXAMPLE
trueError = 0

trueFile = open("../bin/MASTERCNR_EGFRINH_CNRRESULT.pickle", "rb")
sTrue = pickle.load(trueFile)

n_res = np.size(sTrue.residuals_complete) + np.size(sTrue.residuals_incomplete)
error = 0
error += np.sum(np.array(np.square(sTrue.residuals_complete)))/n_res
error += np.sum(np.array(np.square(sTrue.residuals_incomplete)))/n_res
trueError = error

#safe results for later loading
result.theta = THETA
result.trueError = trueError
result.randomErrorList = errorList

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "PERMUTATION_DATA_EGFR_INHTEST"
pickleFile = folder + name + ".pickle"
pickle.dump(result, open(pickleFile, 'wb'))