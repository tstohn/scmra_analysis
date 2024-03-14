
#generate data do compare network quality
import sys
sys.path.insert(0, '../../helperScripts/')
sys.path.insert(0,'/home/t.stohn/Tools/CPLEX_Studio201/cplex/python/3.8/x86-64_linux')

from helperFunctions import *

#the mutant networks have 12 and 13 edges in the true network, for the analysis we try to solve
# for those networks


#function to simulate a CNR with a defined set of populations among the possible ones=[wt, ras, braf]
def simulate_cnr(cellPopulationList, noise, cell, f, count, result):
    #generate the raw data
    pertIndices = None
    rglob, rtot, group_annot, cell_annot, tx_annot = generate_cnr_data(cell, pertIndices, noise, cellPopulationList)

    theta = 0.0
    scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)

    #eta must be estimated for every simulation bcs. we can have different missing proteins
    #etasInput = [0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.7, 0.9]
    eta, possibleSolution, numEdges = get_suitable_solution_fast(rglob, rtot, cell_annot, tx_annot, group_annot, reconstructionType = "CNR", expectedEdges = 13)

    #update logfile
    f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(count) + 
            ", etas: " + str(round(eta,3)) )
    
    populationString = ""
    for population in cellPopulationList:
        populationString += ','
        populationString += population
    result.add_parameter([noise, cell, eta, count, numEdges, "CNR", populationString])

    f.write(" edges: " + str(numEdges) + "\n")
    f.flush()

    #add simulation output to result
    result.add_result(possibleSolution)
################################################
# CNR SIMULATIONS
################################################

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0.2]
#CELLNUM IS NUMBER PER SUBGROUP: 100 = 2*100 = 200
cellNumArray = [10, 100, 1000]
repeats = 20

#cnrSimulations = [{"N":noiseArray, "CN":[1000]}, {"N":[0.2], "CN":cellNumArray}]
cnrSimulations = [{"N":[0.2], "CN":cellNumArray}]

#add header to result
header = ['NOISE', 'CELLNUM', 'ETA', 'REPEAT', 'EDGES', 'MODELTYPE','POPULATIONS']
result = MRASimulationClass(header)

logFile = "../bin/logCnrNoiseCellNum_FINAL.txt"
f = open(logFile, "w+")

#Simulations to run
# 1.) CNR of ABC, AB, AC, BC
# 2.) MRA of A, B, C with 500 each and 1500 eACH

for sim in cnrSimulations:
    noiseArray = sim["N"]
    cellNumArray = sim["CN"]
    for noise in noiseArray:
        for cellNum in cellNumArray:
            for count in range(repeats):
                #CNR RUNS of all combinations
                #simulate_cnr(["braf", "ras", "wt"], noise, cellNum, f, count, result)
                simulate_cnr(["braf", "wt"], noise, cellNum, f, count, result)
                #simulate_cnr(["ras", "wt"], noise, cellNum, f, count, result)

f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_CNR_Noise_Cellnum_FINAL"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))