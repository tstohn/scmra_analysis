'''
THIS SCRIPT GENERATES MRA DATA FOR TWO MODELS: A PRIOR-KNOWLEDGE and an UN BIASED MODEL

In this simulation we randomly samples true prior-knowledge edges and gradually increase
the number of prior knowledge edges to see if more prior-knowledge gradually improves network reconstruction

'''

#generate data do compare network quality
import sys
sys.path.insert(0, '../..') #path to helperFunctions
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
import scmra

ORTON_EDGES = list([
    # MAPK pathway
    ("Ras", "Sos"), ("Raf1","Ras"),
    ("PI3K", "Ras"), ("Akt", "PI3K"),("Raf1", "Akt"),
    ("bRaf", "Ras"), ("Mek", "Raf1"),
    ("Rap1", "C3G"), ("bRaf", "Rap1"),("Mek", "bRaf"),("Erk", "Mek"),("P90Rsk","Erk"),("Sos","P90Rsk")
])

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noise=0.2
#cellNumArray = [10,50,100,500,1000]
cell = 1000

repeats = 1
#gradually increase the number of prior-knowledge edges from 1 to 13
numberPriorKnowledgeEdges = list(range(1, 14))

#add header to result
header = ['NOISE', 'CELLNUMBER', 'ETA','EDGENUMBER', 'REPEAT', 'MODEL','NUMBER_PRIOR_EDGES']
result = MRASimulationClass(header)
logFile = "../bin/log4b.txt"
f = open(logFile, "a")

#FOR EVERY REPEAT: MAKE A BASE-LINE MODEL (no prior knowledge) 
# AND SEVERAL MODELS OF INCREASING NUMBER OF PRIOR KNOWLEDGE
for repeat in range(repeats):

    count = 0

    #create data for this repeat
    count += 1
    rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)
    scd = scmra.ScData(rglob=rglob, rtot = rtot)

    #UNBIASED MODEL
    etaU, scMraResult, edgeNumberU = get_suitable_solution_fast(rglob, rtot, etaGuesses=30, priorNetwork = None)
    #add values to result
    result.add_parameter([noise, cell, etaU, edgeNumberU, repeat, 'NOVEL', "0"])
    #add simulation output to result
    result.add_result(scMraResult)

    for pkEdgeNum in numberPriorKnowledgeEdges:

        #sample pkedges
        pkEdges = random.sample(ORTON_EDGES, pkEdgeNum)
        #PRIOR KNowLEDGE MODEL
        etaPN = 0.0
        etaPN, scMraResultPN, edgeNumberPN = get_suitable_solution_fast(rglob, rtot, etaGuesses=30, priorNetwork = pkEdges)        
        edgeNumberPN = scMraResultPN.imap.sum().sum()
        result.add_parameter([noise, cell, etaPN, edgeNumberPN, repeat,'PRIOR_KNOWLEDGE', pkEdgeNum])
        #add simulation output to result
        result.add_result(scMraResultPN)
    
        #update logfile
        f.write(str(count) + ".) SIMULATION: cells=" + str(cell) + "    EDGE-Numbers(unbiased, prior-network): (" + str(edgeNumberU)+ "," + str(edgeNumberPN) + ")  REPEAT:" + str(repeat) + "\n")
        f.flush()
        
f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "GRADUAL_PRIOR_KNOWLEDGE_EFFECT_1000cells"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))