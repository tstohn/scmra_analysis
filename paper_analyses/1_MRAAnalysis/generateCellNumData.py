'''
THIS SCRIPT GENERATEA MRA DATA FOR DIFFERENT  CELLNUMBER and eta with constant noise
to estimate needed cell number to ifner network topology

etas are chosen by firstly estimating the eta-edgeNumber function (stretched exponential decay)
and choosing the etas resultign in networks with edges between 1-22 (two times the real edge number), assuming
that for this maximum edge number almost all true edges are found (TPR=100%)

'''



#generate data do compare network quality
import sys
sys.path.insert(0, '../../helperScripts/')
sys.path.insert(0,'/home/t.stohn/Tools/CPLEX_Studio201/cplex/python/3.8/x86-64_linux')

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

def estimate_eta_list(cell, noise, etas, withFlattening = False):
    rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)
    eta, etaList = estimate_eta(rglob, rtot, etasInput = etas, withFlattening = withFlattening)
    return(etaList)

#we use a gradient of etas to generate ROC curves, adressing how cellnumber affects
#as the eta might not reflect the network complexity well, instead of having an array of etas
#we generate etas for different network complexiites (edges in netowork)
def refineEta(etaList, lastEta, edgesExpected, edgesFound, etaChangeIdx, etaRefineArray, edgeRefineArray, inBetween):
    edgesExpectedIndex = edgesExpected - 1
    simplyPushNewEtaBack = 0
    if(edgesExpected < edgesFound and etaChangeIdx <= 0 and not inBetween):
        etaChangeIdx -= 1
        if(edgesExpectedIndex+etaChangeIdx < 0): 
            newEta = (lastEta + (1-lastEta)/2) #we have already the biggest eta in our list (with smallest edges number), so we need to increase last eta
        else:
            newEta = (etaList[edgesExpectedIndex + etaChangeIdx])
        #f.write("   - a")

    elif(edgesExpected > edgesFound and etaChangeIdx >= 0 and not inBetween):
        etaChangeIdx += 1
        if(edgesExpectedIndex+etaChangeIdx >= (len(etaList))): 
            newEta = (lastEta/2)
        else:
            newEta = (etaList[edgesExpectedIndex + etaChangeIdx])
        #f.write("   - b")
    else: #if last eta was in between the previous ones
        #f.write("   - in between")
        inBetween = True
        #assert(edgesExpected[0] != -1 and edgesExpected[1] != -1)
        #assert((edgesExpected[0]-edgesExpected[1])^2 < (edgesFound-edgesExpected[0])^2) #check that we really have the case of the edge beeing in between
        #we have the 2. and 3. last etas a,b in the array {a, x, b}
        # as well as our very last eta x, but we do not know if we came from an edge increasing or decreasing eta
        #therefore we need to check if expectedNumberOfEdges is between a,x or between x,b
        if( (edgesExpected-edgesFound)*(edgesExpected-edgeRefineArray[0]) < 0): #eta is inbetween last eta (edgesFound) and eta[0]
            newEta = ( etaRefineArray[0] + ((lastEta-etaRefineArray[0])/2) )
            #f.write(" 1")
            simplyPushNewEtaBack = 1
        else:
            newEta = ( etaRefineArray[1] + ((lastEta-etaRefineArray[1])/2) )
            #f.write(" 2")
            simplyPushNewEtaBack = 2

    #the eta and edge list must be adjusted depending on where the newEat is inbetween (between the last and 1 or last and 0 index)
    #in the case where we still did not find an in between value yet we keep on reducing/ increasing the eta and just push back the new values
    if(simplyPushNewEtaBack == 0):
        etaRefineArray[1] = etaRefineArray[0]
        etaRefineArray[0] = lastEta

        edgeRefineArray[1] = edgeRefineArray[0]
        edgeRefineArray[0] = edgesFound
    else:
        #special case where the lastEdge is same as the previous one
        if(edgesFound == edgeRefineArray[0]): #in case the edge number is the same, we can not push it into the array list, otherwise we end up not having upper and lower bounds anymore
            #e.g. having edgesnumber [0,7] and new edge 7 and looking for edge 5, we will not find it between [7,7]
            edgeRefineArray[0] = edgesFound
            etaRefineArray[0] = lastEta
        if(simplyPushNewEtaBack == 1):
            edgeRefineArray[1] = edgesFound
            etaRefineArray[1] = lastEta 
        elif(simplyPushNewEtaBack == 2):
            edgeRefineArray[0] = edgesFound
            etaRefineArray[0] = lastEta 

    return (newEta, etaChangeIdx, inBetween)

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0.2] #unused
cellNumArray = [10,50,100,500,1000]

etaArray = [0.0005, 0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

repeats = 20

#add header to result
header = ['NOISE', 'CELLNUMBER', 'ETA', 'REPEAT']
result = MRASimulationClass(header)
logFile = "../bin/logCellnum2009.txt"
global f 
f= open(logFile, "w+")

noise = noiseArray[0]

for cell in cellNumArray:
    count = 1
    etaList = estimate_eta_list(cell, noise, etaArray, withFlattening = True)

    for repeat in range(repeats):
        expectedEdges = 1
        #we iterate over eta, but actually we are interested in the range of expectedEdges that we increment simultaniously
        for eta in etaList: #etaList is the list of etas for the requested range of edges from the helper function

            rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)
            scd = scmra.ScData(rglob=rglob, rtot = rtot)
            #update logfile
            f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(count) + " etas: " + str(round(eta,3)) + " wantedEdges: " + str(expectedEdges))

            rightEdgeNumber = False
            numberOfTries = 0
            etaChangeIdx = 0
            etaRefineArray = np.array([-1.0, -1.0]) # 0 is new, 1 is old
            edgeRefineArray = np.array([-1.0, -1.0]) # 0 is new, 1 is old
            inBetween = False

            newEta = eta
            while(not rightEdgeNumber and numberOfTries < 50):

                if(eta < 0.): eta = 0. #since eta is approximated by an exponential decay function, set it here

                #make MRA simulation
                p = scmra.ScMraProblem(scd, eta=newEta)
                p.cpx.solve()
                s = scmra.ScMraResult(p)

                #check for the number of edges
                numEdges = 0
                for column_name in s.rloc:
                    column = s.rloc[column_name]
                    numEdges += (column != 0).sum()
                    numEdges -=1 #for the diagonal values

                if(numEdges == expectedEdges or eta == 0.0): 
                    rightEdgeNumber = True
                    break

                numberOfTries += 1
                newEta, etaChangeIdx, inBetween = refineEta(etaList, newEta, expectedEdges, numEdges, etaChangeIdx, etaRefineArray, edgeRefineArray, inBetween)

            #add values to result
            result.add_parameter([noise, cell, newEta, count])

            #add simulation output to result
            result.add_result(s)
            del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)
            f.write(", edges=" + str(numEdges) + "\n")
            f.flush()
            expectedEdges += 1
        count += 1

f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_OF_ETA_CELLNUM_2009"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))