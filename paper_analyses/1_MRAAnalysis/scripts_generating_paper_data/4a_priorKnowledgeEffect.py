'''
THIS SCRIPT GENERATEA MRA DATA FOR TWO MODELS: A PRIOR-KNOWLEDGE and an UN BIASED MODEL

In this simulation we use an unbiased and a model relying only on prior knowledge, and we see how those models deal with different
numbers of cells as input data (same noise for the data, several repeats, where every repeat means a new set of data point)
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
noise=0.5
#cellNumArray = [10,50,100,500,1000]
cellNumArray = [10,25,50,75,100,200,300,400,500,600,700,800,900,1000]

repeats = 20

#add header to result
header = ['NOISE', 'CELLNUMBER', 'ETA','EDGENUMBER', 'REPEAT', 'MODEL']
result = MRASimulationClass(header)
logFile = "../bin/log4a.txt"
f = open(logFile, "a")

for cell in cellNumArray:
    count = 0
    #if we want to estimate eta/ however it looks like it does not 'guess' the node number very well
    #better run a bunch of eta guesses
    #etaList = estimate_eta_list(cell, noise)
    for repeat in range(repeats):
        
        #create data
        count += 1
        rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)

        #UNBIASED MODEL
        scd = scmra.ScData(rglob=rglob, rtot = rtot)
        etaU, scMraResult, edgeNumberU = get_suitable_solution_fast(rglob, rtot, etaGuesses=30, priorNetwork = None)
        #add values to result
        result.add_parameter([noise, cell, etaU, edgeNumberU, repeat, 'NOVEL'])
        #make MRA simulation
        result.add_result(scMraResult)
    
        #PRIOR KNowLEDGE MODEL
        etaPN = 0.0
        p = scmra.ScMraProblem(scd, eta = etaPN, prior_network = ORTON_EDGES)
        p.cpx.solve()
        scMraResultPN = scmra.ScMraResult(p)
        edgeNumberPN = scMraResult.imap.sum().sum()
        result.add_parameter([noise, cell, etaPN, edgeNumberPN, repeat,'PRIOR_KNOWLEDGE'])
        #add simulation output to result
        result.add_result(scMraResultPN)
        del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)
    
        #update logfile
        f.write(str(count) + ".) SIMULATION: cells=" + str(cell) + "    EDGE-Numbers(unbiased, prior-network): (" + str(edgeNumberU)+ "," + str(edgeNumberPN) + ")  REPEAT:" + str(repeat) + "\n")
        f.flush()
        
f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "PRIOR_KNOWLEDGE_EFFECT"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))