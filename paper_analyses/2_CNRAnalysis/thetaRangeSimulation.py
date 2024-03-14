#generate data do compare network quality
import sys
sys.path.insert(0, '../../helperScripts/')
sys.path.insert(0,'/home/t.stohn/Tools/CPLEX_Studio201/cplex/python/3.8/x86-64_linux')
from helperFunctions import *

scriptTitle = "THETA_ESTIMATION_Noise:20_Cellnum:250"

#generate MRA + CNR data for smae population with theta range to
# 1. make theta dependant plot
# 2. estimte good theta for other runs

################################################
# CNR SIMULATIONS
################################################

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0.2] #also tested 0.5 which did not yield great results combined with 2*500 cells
#cellnum is the number of cells per subgroup: in total 500 cells (tested 2*500 before)
cellNumArray = [250]
repeats = 1
thetas = [0.0,0.0001, 0.0005, 0.001, 0.0015, 0.002, 0.0025,0.005, 0.01,0.02, 0.03, 0.06, 0.1, 0.25, 0.3, 0.4, 0.5, 0.75,1]

#add header to result
header = ['NOISE', 'CELLNUM', 'ETA', 'REPEAT', 'EDGES', 'MODELTYPE','POPULATIONS','NUMDIFFEDGES','THETA','RESIDUALS']
result = MRASimulationClass(header)

logFile = "./bin/log_" + scriptTitle + ".txt"
f = open(logFile, "w+")

#Focus on one cell comparison: the difference between Braf - wt was slightly smaller,
#therefore taking this one here to proof CNR outperforms MRA
#choosing 500 cells in each population with a noise 0f 0.2
for noise in noiseArray:
    for cellNum in cellNumArray:
        for count in range(repeats):
  
            #BRAF SIMULATIONS
            rglobWT, rtotWT, rglobMUT, rtotMUT = generate_raw_data_for_MRA_and_CNR(cellNum, noise, "braf")
            #for this simulation we add no treatment/ perturbations
            cell_annot = None
            tx_annot = None
            #combine the mutant & wt data for CNR
            rglobList = [rglobWT, rglobMUT]
            rtotList = [rtotWT, rtotMUT]
            #make group annot: 'wt' ->'cellIds', , 'braf' ->'cellIds'
            group_annot = {}
            group_annot["wt"] = list(rglobWT.columns)
            group_annot["braf"] = list(rglobMUT.columns)
            #combine rglob and rtot
            rglob = pd.concat(rglobList, axis=1)
            rtot = pd.concat(rtotList, axis=1) 
            simulate_mra_from_data("wt",noise, cellNum, f, count, result, rglobWT, rtotWT, cell_annot, tx_annot, "WT:wt_braf")
            simulate_mra_from_data("braf",noise, cellNum, f, count, result, rglobMUT, rtotMUT, cell_annot, tx_annot, "BRAF:wt_braf")

            for thetaTmp in thetas:
                simulate_cnr_from_data(["wt","braf"],noise, cellNum, f, count, result, rglob, rtot, cell_annot, tx_annot, group_annot, thetaTmp)

            #RAS SIMULATIONS
            rglobWT, rtotWT, rglobMUT, rtotMUT = generate_raw_data_for_MRA_and_CNR(cellNum, noise, "ras")
            #for this simulation we add no treatment/ perturbations
            cell_annot = None
            tx_annot = None
            #combine the mutant & wt data for CNR
            rglobList = [rglobWT, rglobMUT]
            rtotList = [rtotWT, rtotMUT]
            #make group annot: 'wt' ->'cellIds', , 'braf' ->'cellIds'
            group_annot = {}
            group_annot["wt"] = list(rglobWT.columns)
            group_annot["ras"] = list(rglobMUT.columns)
            #combine rglob and rtot
            rglob = pd.concat(rglobList, axis=1)
            rtot = pd.concat(rtotList, axis=1) 

            simulate_mra_from_data("wt",noise, cellNum, f, count, result, rglobWT, rtotWT, cell_annot, tx_annot, "WT:wt_ras")
            simulate_mra_from_data("ras",noise, cellNum, f, count, result, rglobMUT, rtotMUT, cell_annot, tx_annot, "RAS:wt_ras")
            
            for thetaTmp in thetas:
                simulate_cnr_from_data(["wt","ras"],noise, cellNum, f, count, result, rglob, rtot, cell_annot, tx_annot, group_annot, thetaTmp)

f.close()

#store pickled result
now = datetime.now()
folder = "./bin/"
name = "SIMULATION_" + scriptTitle
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))