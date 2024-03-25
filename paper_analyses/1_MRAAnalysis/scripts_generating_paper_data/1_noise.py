'''
THIS SCRIPT GENERATEA MRA DATA FOR DIFFERENT VALUES OF NOISE, CELLNUMBER
with several repeats and different etas

etas are chosen by firstly estimating the eta-edgeNumber function (stretched exponential decay)
and choosing the etas resultign in networks with edges between 1-22 (two times the real edge number), assuming
that for this maximum edge number almost all true edges are found (TPR=100%)

'''



#generate data do compare network quality
import sys
sys.path.insert(0, '../scripts_generating_paper_data/')
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

def estimate_eta_list(cell, noise):
    rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)
    eta, etaList = estimate_eta(rglob, rtot)
    return(etaList)

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0.0, 0.1, 0.2, 0.5]
#cellNumArray = [10,50,100,500,1000]
cellNumArray = [1000]

repeats = 10

#add header to result
header = ['NOISE', 'CELLNUMBER', 'ETA', 'REPEAT']
result = MRASimulationClass(header)
logFile = "../bin/log.txt"
f = open(logFile, "a")

for noise in noiseArray:
    for cell in cellNumArray:
        count = 0
        #if we want to estimate eta/ however it looks like it does not 'guess' the node number very well
        #better run a bunch of eta guesses
        #etaList = estimate_eta_list(cell, noise)
        for repeat in range(repeats):

            #manually run with a couple of etas
            smallEtas = [10**(i) for i in np.arange(-2,-4,-.1)]
            etaList = [0.025, 0.05, 0.075,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9] + smallEtas

            count += 1
            rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)
            scd = scmra.ScData(rglob=rglob, rtot = rtot)
            #update logfile
            f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(count) + " etas: ")
            for etaTmp in etaList:
                if(etaTmp <=0): break
                f.write(str(round(etaTmp,3)) + " ")

                #add values to result
                result.add_parameter([noise, cell, etaTmp, count])

                #make MRA simulation
                p = scmra.ScMraProblem(scd, eta=etaTmp)
                p.cpx.solve()
                s = scmra.ScMraResult(p)

                #add simulation output to result
                result.add_result(s)
                del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)
            f.write("\n")
            f.flush()

f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_OF_NOISE_CELLNUM"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))