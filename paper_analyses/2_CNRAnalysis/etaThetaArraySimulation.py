#generate data do compare network quality
import sys
sys.path.insert(0, '..') #path to helperFunctions
from helperFunctions import *

scriptTitle = "etaThetaArray_Noise:20_Cellnum:250"

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
repeats = 5

#(0.001 took too long for eta=0 => too many possibilites of diff edges)
thetaArray = [0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1] #decreasing number of diff edges
etaArray = [0.005, 0.0075, 0.01, 0.075, 0.2, 0.4, 0.6] #decreasing edge number

#thetaArray = [0.0005, 0.005,0.2, 0.4, 0.6]
#etaArray = [0.025, 0.1,0.3, 0.6]

#add header to result
header = ['NOISE', 'CELLNUM', 'ETA', 'REPEAT', 'EDGES', 'MODELTYPE','POPULATIONS','NUMDIFFEDGES','THETA','RESIDUALS','MODELCOMPLEXITY','RESIDUALS_TRAIN_CALCULATED','RESIDUALS_TEST_CALCULATED']
result = MRASimulationClass(header)

logFile = "./bin/log_" + scriptTitle + ".txt"
f = open(logFile, "w+")

#Focus on one cell comparison: the difference between Braf - wt was slightly smaller,
#therefore taking this one here to proof CNR outperforms MRA
#choosing 500 cells in each population with a noise 0f 0.2
for noise in noiseArray:
    for cellNum in cellNumArray:
        for count in range(repeats):
  
            # 1.) SET UP THE MODEL, SIMULATE TEST/ TRAIN DATA
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
            
            test_rglobWT, test_rtotWT, test_rglobMUT, test_rtotMUT = generate_raw_data_for_MRA_and_CNR(cellNum, noise, "ras")
            #combine the mutant & wt data for CNR
            rglobList_test = [test_rglobWT, test_rglobMUT]
            rtotList_test = [test_rtotWT, test_rtotMUT]
            #make group annot: 'wt' ->'cellIds', , 'braf' ->'cellIds'
            group_annot_test = {}
            group_annot_test["wt"] = list(test_rglobWT.columns)
            group_annot_test["ras"] = list(test_rglobMUT.columns)
            #combine rglob and rtot
            rglob_test = pd.concat(rglobList_test, axis=1)
            rtot_test = pd.concat(rtotList_test, axis=1) 
            
            for thetaTmp in thetaArray:
                for etaTmp in etaArray:
                    
                    # 2a) RUN CNR ON TRAINING DATA
                    countTmp = str(count)
                    indicatorsFrom = None #Training data to discover its own network
                    tmpResult = simulate_cnr_from_data(["wt","ras"],noise, cellNum, f, countTmp, result, rglob, rtot, cell_annot, tx_annot, group_annot, thetaTmp, etaTmp, indicatorsFrom)
                    #as a quality check calculate residuals for origional data
                    get_residuals_on_new_data(tmpResult, rglob, rtot, group_annot, result)

                    # 2b) RUN CNR ON TEST DATA
                    get_residuals_on_new_data(tmpResult, rglob_test, rtot_test, group_annot_test, result)

f.close()

#store pickled result
now = datetime.now()
folder = "./bin/"
name = "SIMULATION_" + scriptTitle
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))

result.write_results(folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.tsv')