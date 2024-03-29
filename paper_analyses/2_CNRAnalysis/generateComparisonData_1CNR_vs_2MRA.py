#generate data do compare network quality
#while exampleThetaRange iterates over theta ranges, we here iterate over several repeats

scriptTitle = "CnrOverMra_NoiseRange_CellNum250"

import sys
sys.path.insert(0, '..') #path to helperFunctions

from helperFunctions import *

#from exampleThetaRange: 0.01 was roughly 5 to 6 edges with best RMSE of difference
THETA = 0.02

################################################
# CNR SIMULATIONS
################################################

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0.0,0.2,0.5]
#cellnum is the number of cells per subgroup: in total [cellnum * subgroup] cells
cellNumArray = [250]
repeats = 20

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

            simulate_cnr_from_data(["wt","braf"],noise, cellNum, f, count, result, rglob, rtot, cell_annot, tx_annot, group_annot, THETA)
            simulate_mra_from_data("wt",noise, cellNum, f, count, result, rglobWT, rtotWT, cell_annot, tx_annot, "WT:wt_braf")
            simulate_mra_from_data("braf",noise, cellNum, f, count, result, rglobMUT, rtotMUT, cell_annot, tx_annot, "BRAF:wt_braf")


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

            simulate_cnr_from_data(["wt","ras"],noise, cellNum, f, count, result, rglob, rtot, cell_annot, tx_annot, group_annot, THETA)
            simulate_mra_from_data("wt",noise, cellNum, f, count, result, rglobWT, rtotWT, cell_annot, tx_annot, "WT:wt_ras")
            simulate_mra_from_data("ras",noise, cellNum, f, count, result, rglobMUT, rtotMUT, cell_annot, tx_annot, "RAS:wt_ras")

f.close()

#store pickled result
now = datetime.now()
folder = "./bin/"
name = "SIMULATION_" + scriptTitle
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))