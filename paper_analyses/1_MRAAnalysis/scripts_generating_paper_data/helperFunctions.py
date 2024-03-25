from curses import setupterm
from sys import settrace
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import numpy as np
import sympy
import random
import pyprojroot
import os
import pickle
import math
import pickle
import seaborn as sns; sns.set_theme()
from matplotlib.pyplot import figure
from itertools import product
from sklearn.metrics import mean_squared_error 
from numpy import NaN, nan
import os
import matplotlib
from datetime import datetime
from collections import defaultdict
from copy import deepcopy

import scmra

dirname = os.path.dirname(__file__)

class MRASimulationClass:
    
    def __init__(self, header):
        self._headerList = header #list of the data types stored in the list above (e.g. 'noise', 'CellNum',...)
        self._parameterList = [] #list of list with parameter elements (e.g. a list for each simulation, which is a list of noise, cellNUmber, eta, etc.)
        self._resultList = [] #list of result elements

    @property 
    def header(self):
        return(self._headerList)
    @property 
    def parameter(self):
        return(self._parameterList)
    @property 
    def result(self):
        return(self._resultList)

    @header.setter
    def header(self, header):
        self._headerList = header
    def add_parameter(self, param):
        self._parameterList.append(param)
    def add_result(self, newResult):
        self._resultList.append(newResult)

    #returns results from this object
    #keylist stores the key for which to check
    #and valueList contains the correcponding values for the keys (in same order!)
    def getResult(self, keylist, valuelist, checkForSingleResult = True, allowEmptyResults = True):
        headerPositionDict = {}
        results = []
        params = []
        for headerCol in keylist:
            pos = self._headerList
            pos = pos.index(headerCol)
            headerPositionDict[headerCol] = pos

        variablesToCompare = len(keylist)
        simulationLength = len(self._parameterList)
        for simulationIdx in range(simulationLength):
            simulation = self._parameterList[simulationIdx]
            foundSolution = True
            for idx in range(variablesToCompare):
                header = keylist[idx]
                value = valuelist[idx]
                if(simulation[headerPositionDict[header]] != value): 
                    foundSolution = False
                    break
            if(foundSolution):
                results.append(self._resultList[simulationIdx])
                params.append(self._parameterList[simulationIdx])

        if(checkForSingleResult):
            assert(len(results) < 2)    
            if(allowEmptyResults):
                if(len(results) == 0):
                    return None, None
                else:
                    return(results[0], params[0]) 
            else:
                assert(len(results) == 1)    
                return(results[0], params[0])

        #if we do not check for single result return all found ones
        return(results, params)

class MraParameter:
    eta = None
    etaGraph = None
    logFile = None
    modelPertubationValues = False

    def __init__(self, logFileInput = None):
        self.logFile = logFileInput.logFile
        self.folder = logFileInput.folder

class LogFile:

    logFile = None
    folder = None

    def __init__(self, folder = None, file = None):
        self.logFile = folder + "/" + file
        self.folder = folder

        file = open(self.logFile,'a+')
        file.close()

class SimulationResults:

    logFile = None
    resultDF = None
    logFolder = None

    def __init__(self, folder):
        self.logFolder = folder
        self.logFile = folder + "/log.txt"
        column_names = ["SimulationName", "RMSE_rloc", "TP_rloc", "edges"]
        self.resultDF = pd.DataFrame(columns = column_names)

        f = open(self.logFile, "w")
        f.close()
        #create subdir for eta-estimation plots
        path = folder + '/etas/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        #create subdir for scMRA result pickles
        path = folder + '/resultObjects/'
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)

#loading all the data
# Load the ground-truth data
DAT_DIR = pyprojroot.here("simulations/true_network_parameters/")

rloc_true = pd.read_csv(os.path.join(DAT_DIR, "rloc-true.tsv"), 
                        sep='\t', header=0, index_col=0)
stot_true = pd.read_csv(os.path.join(DAT_DIR, "stot-true.tsv"), sep='\t', 
                        header=0, index_col=0)
# Determine the set of true interactions (edges)
imap_filtered = 1*(abs(rloc_true) > 0.01 + np.identity(rloc_true.shape[0]))
imap_true = 1*(abs(rloc_true) > 0 + np.identity(rloc_true.shape[0]))

NODES = tuple(imap_true.columns)

dat_dir = pyprojroot.here("data/simulated_data")
# Unperturbed cells - 1000 cells
wt = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1.tsv"), sep='\t')
wt.index = ["ctr."+str(i) for i in range(wt.shape[0])]

braf = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut.tsv"), sep='\t')
braf.index = ["braf."+str(i) for i in range(braf.shape[0])]

ras = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut.tsv"), sep='\t')
ras.index = ["ras."+str(i) for i in range(ras.shape[0])]

#two dataframes with all the inhibited populations, one for each population of cell types (wt and braf)
PerturbDataList = {}
PerturbDataListbraf = {}
PerturbDataListras = {}

orton_RASi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASi-delta=0.25.tsv"), sep='\t')
orton_RASi.index = ["Rasi."+str(i) for i in range(orton_RASi.shape[0])]
PerturbDataList["Rasi"]=(orton_RASi)

orton_Akti = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Akti-delta=0.25.tsv"), sep='\t')
orton_Akti.index = ["Akti."+str(i) for i in range(orton_Akti.shape[0])]
PerturbDataList["Akti"]=(orton_Akti)

orton_bRafi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-bRafi-delta=0.25.tsv"), sep='\t')
orton_bRafi.index = ["bRafi."+str(i) for i in range(orton_bRafi.shape[0])]
PerturbDataList["bRafi"]=(orton_bRafi)

orton_C3Gi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-C3Gi-delta=0.25.tsv"), sep='\t')
orton_C3Gi.index = ["C3Gi."+str(i) for i in range(orton_C3Gi.shape[0])]
PerturbDataList["C3Gi"]=(orton_C3Gi)

orton_Erki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Erki-delta=0.25.tsv"), sep='\t')
orton_Erki.index = ["Erki."+str(i) for i in range(orton_Erki.shape[0])]
PerturbDataList["Erki"]=(orton_Erki)

orton_Meki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Meki-delta=0.25.tsv"), sep='\t')
#orton_Meki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Meki.tsv"), sep='\t')
orton_Meki.index = ["Meki."+str(i) for i in range(orton_Meki.shape[0])]
PerturbDataList["Meki"]=(orton_Meki)

orton_P90Rski = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-p90Rski-delta=0.25.tsv"), sep='\t')
orton_P90Rski.index = ["P90Rski."+str(i) for i in range(orton_P90Rski.shape[0])]
PerturbDataList["P90Rski"]=(orton_P90Rski)

orton_PI3Ki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-PI3Ki-delta=0.25.tsv"), sep='\t')
orton_PI3Ki.index = ["PI3Ki."+str(i) for i in range(orton_PI3Ki.shape[0])]
PerturbDataList["PI3Ki"]=(orton_PI3Ki)

orton_Raf1i = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Raf1i-delta=0.25.tsv"), sep='\t')
orton_Raf1i.index = ["Raf1i."+str(i) for i in range(orton_Raf1i.shape[0])]
PerturbDataList["Raf1i"]=(orton_Raf1i)

orton_Rap1i = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Rap1i-delta=0.25.tsv"), sep='\t')
orton_Rap1i.index = ["Rap1i."+str(i) for i in range(orton_Rap1i.shape[0])]
PerturbDataList["Rap1i"]=(orton_Rap1i)

orton_Sosi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Sosi-delta=0.25.tsv"), sep='\t')
orton_Sosi.index = ["Sosi."+str(i) for i in range(orton_Sosi.shape[0])]
PerturbDataList["Sosi"]=(orton_Sosi)

#inhibition of braf populations
orton_RASi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Rasi.tsv"), sep='\t')
orton_RASi.index = ["Rasi."+ "braf" + str(i) for i in range(orton_RASi.shape[0])]
PerturbDataListbraf["Rasi"]=(orton_RASi)

orton_Akti = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Akti.tsv"), sep='\t')
orton_Akti.index = ["Akti."+ "braf" + str(i) for i in range(orton_Akti.shape[0])]
PerturbDataListbraf["Akti"]=(orton_Akti)

orton_bRafi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-bRafi.tsv"), sep='\t')
orton_bRafi.index = ["bRafi." + "braf" +str(i) for i in range(orton_bRafi.shape[0])]
PerturbDataListbraf["bRafi"]=(orton_bRafi)

orton_C3Gi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-C3Gi.tsv"), sep='\t')
orton_C3Gi.index = ["C3Gi." + "braf" +str(i) for i in range(orton_C3Gi.shape[0])]
PerturbDataListbraf["C3Gi"]=(orton_C3Gi)

orton_Erki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Erki.tsv"), sep='\t')
orton_Erki.index = ["Erki." + "braf" +str(i) for i in range(orton_Erki.shape[0])]
PerturbDataListbraf["Erki"]=(orton_Erki)

orton_Meki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Meki.tsv"), sep='\t')
#orton_Meki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-Meki.tsv"), sep='\t')
orton_Meki.index = ["Meki." + "braf" +str(i) for i in range(orton_Meki.shape[0])]
PerturbDataListbraf["Meki"]=(orton_Meki)

orton_P90Rski = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-p90Rski.tsv"), sep='\t')
orton_P90Rski.index = ["P90Rski." + "braf" +str(i) for i in range(orton_P90Rski.shape[0])]
PerturbDataListbraf["P90Rski"]=(orton_P90Rski)

orton_PI3Ki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-PI3Ki.tsv"), sep='\t')
orton_PI3Ki.index = ["PI3Ki." + "braf" +str(i) for i in range(orton_PI3Ki.shape[0])]
PerturbDataListbraf["PI3Ki"]=(orton_PI3Ki)

orton_Raf1i = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Raf1i.tsv"), sep='\t')
orton_Raf1i.index = ["Raf1i." + "braf" +str(i) for i in range(orton_Raf1i.shape[0])]
PerturbDataListbraf["Raf1i"]=(orton_Raf1i)

orton_Rap1i = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Rap1i.tsv"), sep='\t')
orton_Rap1i.index = ["Rap1i." + "braf" +str(i) for i in range(orton_Rap1i.shape[0])]
PerturbDataListbraf["Rap1i"]=(orton_Rap1i)

orton_Sosi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-BRAFmut-Sosi.tsv"), sep='\t')
orton_Sosi.index = ["Sosi." + "braf" +str(i) for i in range(orton_Sosi.shape[0])]
PerturbDataListbraf["Sosi"]=(orton_Sosi)

#inhibition of ras populations
orton_RASi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Rasi.tsv"), sep='\t')
orton_RASi.index = ["Rasi."+ "ras" + str(i) for i in range(orton_RASi.shape[0])]
PerturbDataListras["Rasi"]=(orton_RASi)

orton_Akti = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Akti.tsv"), sep='\t')
orton_Akti.index = ["Akti."+ "ras" + str(i) for i in range(orton_Akti.shape[0])]
PerturbDataListras["Akti"]=(orton_Akti)

orton_bRafi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-bRafi.tsv"), sep='\t')
orton_bRafi.index = ["bRafi." + "ras" +str(i) for i in range(orton_bRafi.shape[0])]
PerturbDataListras["bRafi"]=(orton_bRafi)

orton_C3Gi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-C3Gi.tsv"), sep='\t')
orton_C3Gi.index = ["C3Gi." + "ras" +str(i) for i in range(orton_C3Gi.shape[0])]
PerturbDataListras["C3Gi"]=(orton_C3Gi)

orton_Erki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Erki.tsv"), sep='\t')
orton_Erki.index = ["Erki." + "ras" +str(i) for i in range(orton_Erki.shape[0])]
PerturbDataListras["Erki"]=(orton_Erki)

orton_Meki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Meki.tsv"), sep='\t')
orton_Meki.index = ["Meki." + "ras" +str(i) for i in range(orton_Meki.shape[0])]
PerturbDataListras["Meki"]=(orton_Meki)

orton_P90Rski = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-p90Rski.tsv"), sep='\t')
orton_P90Rski.index = ["P90Rski." + "ras" +str(i) for i in range(orton_P90Rski.shape[0])]
PerturbDataListras["P90Rski"]=(orton_P90Rski)

orton_PI3Ki = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-PI3Ki.tsv"), sep='\t')
orton_PI3Ki.index = ["PI3Ki." + "ras" +str(i) for i in range(orton_PI3Ki.shape[0])]
PerturbDataListras["PI3Ki"]=(orton_PI3Ki)

orton_Raf1i = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Raf1i.tsv"), sep='\t')
orton_Raf1i.index = ["Raf1i." + "ras" +str(i) for i in range(orton_Raf1i.shape[0])]
PerturbDataListras["Raf1i"]=(orton_Raf1i)

orton_Rap1i = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Rap1i.tsv"), sep='\t')
orton_Rap1i.index = ["Rap1i." + "ras" +str(i) for i in range(orton_Rap1i.shape[0])]
PerturbDataListras["Rap1i"]=(orton_Rap1i)

orton_Sosi = pd.read_csv(os.path.join(dat_dir, "steady-state-sigma=0.1-RASmut-Sosi.tsv"), sep='\t')
orton_Sosi.index = ["Sosi." + "ras" +str(i) for i in range(orton_Sosi.shape[0])]
PerturbDataListras["Sosi"]=(orton_Sosi)

#further variables of wild type and mutated populations
imap_true = 1*(abs(rloc_true) > 0 + np.identity(rloc_true.shape[0]))
#having list of dataframes and list of all nodes in same order for subsetting by index
PerturbDataList = [v for k,v in (sorted(PerturbDataList.items()))]
PerturbDataListbraf = [v for k,v in (sorted(PerturbDataListbraf.items()))]
PerturbDataListras = [v for k,v in (sorted(PerturbDataListras.items()))]

allNodes = sorted(imap_true.columns)
rloc_true = pd.read_csv(os.path.join(dirname,"simulations/true_network_parameters/rloc-true.tsv"), sep='\t', header=0, index_col=0)
rloc_braf = pd.read_csv(os.path.join(dirname,"simulations/true_network_parameters/rloc-braf-true.tsv"), sep='\t', header=0, index_col=0)
rloc_ras = pd.read_csv(os.path.join(dirname,"simulations/true_network_parameters/rloc-ras-true.tsv"), sep='\t', header=0, index_col=0)

def reorder_rloc(rlocOld):
    colOrder = rloc_true.columns
    rowOrder = rloc_true.index
    rlocNew = rlocOld.reindex(columns=colOrder)
    rlocNew = rlocNew.reindex(rowOrder)

    return(rlocNew)
rloc_braf = reorder_rloc(rloc_braf)
rloc_ras = reorder_rloc(rloc_ras)

imap_true = 1*(abs(rloc_true) > 0 + np.identity(rloc_true.shape[0]))

#rloc dict
rlocDict = {"wt": rloc_true, "braf" : rloc_braf, "ras" : rloc_ras}

#cellTypeDict: mapping of string for cell type (like wt or braf) to the actuall population
cellTypeDict = {"wt": wt, "braf" : braf, "ras" : ras}

def get_to_dict(matrix):
    edgeDict = {}
    rows = list(matrix.index)
    for col in matrix.columns:
        colVec = matrix[col]
        edges = []
        for i in range(0, len(rows)):
            if(colVec[i] == 1):
                edges.append(rows[i])
        edgeDict[col] = edges    
    return(edgeDict)

def get_from_dict(matrix):
    edgeDict = {}
    cols = list(matrix.columns)
    for row in matrix.index:
        colVec = matrix[row]
        edges = []
        for i in range(0, len(cols)):
            if(colVec[i] == 1):
                edges.append(cols[i])
        edgeDict[row] = edges    
    return(edgeDict)

#returns None oi the subset is None
def get_index_list(list_all, list_subset):
    if(list_subset is None): return(None)
    return_list = []
    for el in list_subset:
        x = list_all.index(el)
        return_list.append(x)
    return(return_list)

toDict = get_to_dict(imap_true)
fromDict = get_from_dict(imap_true)

#get a vairable with a certain prefix from singleCellProblem (not from the result)
def get_variable(scp, solidx=0, prefix=""):
    allvars = scp.cpx.variables.get_names()
    vars_lst = [var for var in allvars if var.startswith(prefix)]
    vars_vals = scp.cpx.solution.pool.get_values(solidx, vars_lst)
    vars_dict = dict(zip(vars_lst, vars_vals))

    return(vars_dict)

def sample_raw_mra_data(data = wt, indexListOfPertData=None, cellNum=NaN, perturbedData = PerturbDataList, scaleByTxGroupIndividually = False):
    dat = None
    pertubeData = None
    cell_annot = None
    tx_annot = None
    maxCellNum = min([len(data.index),len(perturbedData[0].index)])
    if((indexListOfPertData is not None) and (indexListOfPertData)):
        assert(max(indexListOfPertData) < len(perturbedData))

        limit = int(maxCellNum/(len(indexListOfPertData)+1))
        if(not math.isnan(cellNum)):
            limit = int(cellNum/(len(indexListOfPertData)+1))

        pertubeData = [perturbedData[i].sample(n = limit) for i in indexListOfPertData]

        pertubeData.append(data.sample(n = limit))
        dat = pd.concat(pertubeData, axis=0)
        
        # # # Mapping of treatments to cells and treatment to target
        cell_annot = {}
        treatmentList = []
        for s in list(dat.index):
            key = s.split(".")[0]
            treatmentList.append(key)
            if key in cell_annot:
                # append the new cell to the existing array
                cell_annot[key].append(s)
            else:
                # create a new array
                cell_annot[key] = [s]

        #add keys to treament annotation list (remove i from back only)
        tx_annot = {s:s[:-1] for s in treatmentList}
        if(data.equals(wt)):
            tx_annot["ctr"] = None
        elif(data.equals(braf)):
            tx_annot["braf"] = None
        elif(data.equals(ras)):
            tx_annot["ras"] = None

    elif(cellNum is not NaN):
        dat = data.sample(n = cellNum)
    else:
        dat = data.sample(n = maxCellNum)
    return(dat, cell_annot, tx_annot)

#only for the wildtype population (not for braf mutated population)
def combine_control_and_perturbed_cells(indexListOfPertData, noise = 0.0, cellNum=NaN, scaleByTxGroupIndividually = False, cellPop = "wt"):
    # Combine data
    if(cellPop == "wt"):
        dat, cell_annot, tx_annot = sample_raw_mra_data(wt, indexListOfPertData, cellNum, PerturbDataList, scaleByTxGroupIndividually)
    elif(cellPop == "braf"):
        dat, cell_annot, tx_annot = sample_raw_mra_data(braf, indexListOfPertData, cellNum, PerturbDataList, scaleByTxGroupIndividually)
    elif(cellPop == "ras"):
        dat, cell_annot, tx_annot = sample_raw_mra_data(ras, indexListOfPertData, cellNum, PerturbDataList, scaleByTxGroupIndividually)

    # Get rglob and rtot, use mean of CTR expression as reference 
    rglob, rtot, pDict = prepare_orton_data(dat,'all', noise, cell_annot, tx_annot, 'ctr', scaleByTxGroupIndividually)
    return(rglob, rtot, cell_annot, tx_annot, pDict)

def performance_metrics(imap, ref_imap):
    """Caluculate performance metrics
    
    True positive rate (tpr): True positive / Condition positive
    False positive rate (fpr): False positive / Condition negative
    Precision : True positive / Predicted positive
    Specificity: True negative / Condition negative
    
    
    """
    assert imap.shape == ref_imap.shape
    assert imap.shape[0] == imap.shape[1]
    n_diag = imap.shape[0]
    tp = ((imap == 1) & (ref_imap == 1)).sum().sum() # True positive
    fp = ((imap == 1) & (ref_imap == 0)).sum().sum() # False positive
    tn = ((imap == 0) & (ref_imap == 0)).sum().sum() - n_diag # True negative
    fn = ((imap == 0) & (ref_imap == 1)).sum().sum() # False negative

    return {"tp": tp, "fp": fp, "tn": tn, "fn" : fn}

def add_noise(mat_input, noise=.1):
    """Add noise to the elements of a matrix.
​
    Matrix should be either numpy array pandas DataFrame.
    Noise is added by multiplying values with a random variable with mu = 1.
    Sigma van be given as input. Default value is 0.1.
    Returns a copy of the matrix. Input is unchanged.
    """
    mat = mat_input.copy()
    dims = mat.shape
    for i in range(dims[0]):
        for j in range(dims[1]):
            if isinstance(mat, np.ndarray):
                mat[i][j] = mat[i][j] * random.normalvariate(1, noise)
            elif isinstance(mat, pd.DataFrame):
                mat.iloc[i, j] = mat.iloc[i, j] * random.normalvariate(1, noise)
            else:
                raise TypeError("argument mat should be numpy array or pandas \
                DataFrame.")
    return(mat)

def prepare_orton_data(dat, n_cells='all', noise=0, cell_annot=None, tx_annot = None, populationsType = 'ctr', scaleByTxGroupIndividually = False, noiseOnRawData=False):

    dat_scaled = None
    perturbedNodeDict = {}
    if n_cells == 'all':
        dat_sample = dat.copy(deep=True)
    else:
        dat_sample = dat.copy(deep=True).sample(n_cells)

    #ADD NOISE ALREADY HERE, NOT AFTER SCALING (after scaling introduces different variabilities for perturbations data
    # as we scale perturbed group by control - this means the noise is based on difference to ct group not reflecting true variability of the raw data)
    if(noiseOnRawData):
        dat_sample = add_noise(dat_sample, noise)

    if cell_annot is None:
        dat_scaled = ((dat_sample - dat_sample.median())/dat_sample.median()).transpose()
    else:

        #scale by difference to ctr group
        if(not scaleByTxGroupIndividually):
            ctr_cells = cell_annot[populationsType]
            dat_ctr = dat_sample.loc[ctr_cells]
            dat_scaled = ((dat_sample - dat_ctr.median())/dat_ctr.median()).transpose()
        else:

            ctr_cells = cell_annot[populationsType]
            dat_ctr = dat_sample.loc[ctr_cells]
            dat_ctr_scaled = ((dat_ctr - dat_ctr.median())/dat_ctr.median())

            #scaling each subgroup seperately
            for tx, node in tx_annot.items():
                tx_cells = cell_annot[tx]
                dat_tx = dat_sample.loc[tx_cells]

                #scale by subgroup
                #dat_scaled_tmp = ((dat_tx - dat_tx.median())/dat_tx.median())

                #crazy scaling theme (minus ctr and div by scaled)
                dat_scaled_tmp = ((dat_tx - dat_tx.median())/dat_tx.median())

                dat_tmp_scaled_to_ctr = ((dat_tx - dat_ctr.median())/dat_ctr.median())

                if(tx != populationsType):
                    diff = dat_tmp_scaled_to_ctr.median()
                    dat_scaled_tmp = dat_scaled_tmp + diff

                #alternative to the uncommented code
                dat_scaled_tmp = dat_scaled_tmp.transpose() 

                #add new scaled data to final dataframe
                if(dat_scaled is None):
                    dat_scaled = dat_scaled_tmp
                else:
                    dat_scaled = pd.concat([dat_scaled, dat_scaled_tmp], axis=1)

    if(not noiseOnRawData):
        dat_noise = add_noise(dat_scaled, noise)
    else:
        dat_noise = dat_scaled
    
    rglob = dat_noise.loc[[n for n in list(dat_noise.index) if ("Active" in n)]]
    rglob.index = [i.replace("Active", "") for i in rglob.index]
    
    rtot = dat_noise.loc[[n for n in list(dat_noise.index) if ("Tot" in n)]]
    rtot.index = [i.replace("Tot", "") for i in rtot.index]

    #reorder columns to have wt and perturbed mixed in columns
    columnReorder = rglob.columns.values.tolist()
    random.shuffle(columnReorder)

    rglob = rglob[columnReorder]
    rtot = rtot[columnReorder]

    return rglob, rtot, perturbedNodeDict

#CNR data generation with perturbaitons and different cell populations (for now only wt and braf mutant)
#CELLNUM is the number of cells in EACH POPULATION seperately, so cellnum=1000 will e.g. make 2000 cells if we model WT and BRAF
def generate_cnr_data(cellNum, perturbedNodeIndices, noise, populationList = ["wt", "braf", "ras"]):
    modelWt = False
    modelBraf = False
    modelRas = False

    for el in populationList:
        if(el not in ["wt", "braf", "ras"]):
            print("ERROR IN POPULATION LIST\n")
            exit()

    #prepare orton data for every inhibition/ cell population
    cellAnnotList = []
    txAnnotList = []
    if("wt" in populationList):
        modelWt = True
        datWt, cell_annotWt, tx_annotWt = sample_raw_mra_data(wt, perturbedNodeIndices, cellNum, PerturbDataList)
        cellAnnotList.append(cell_annotWt)
        txAnnotList.append(tx_annotWt)
    if("braf" in populationList):
        modelBraf = True
        datBraf, cell_annotBraf, tx_annotBraf = sample_raw_mra_data(braf, perturbedNodeIndices, cellNum, PerturbDataListbraf)
        cellAnnotList.append(cell_annotBraf)
        txAnnotList.append(tx_annotBraf)
    if("ras" in populationList):
        modelRas = True
        datRas, cell_annotRas, tx_annotRas = sample_raw_mra_data(ras, perturbedNodeIndices, cellNum, PerturbDataListras)
        cellAnnotList.append(cell_annotRas)
        txAnnotList.append(tx_annotRas)

    #combine the cell annotations to one
    #cell annot form: 'ctr' -> 'cellIds', 'rasi' -> 'cellIds'
    cell_annot = defaultdict(list)
    for d in cellAnnotList: # you can list as many input dicts as you want here
        if(not d): 
            cell_annot = None
            continue
        for key, value in d.items():
            for cellID in value:
                cell_annot[key].append(cellID)

    #tx_annot = {'ctr': None, 'rasi' : 'Ras'}
    tx_annot = defaultdict(list)
    for d in txAnnotList: # you can list as many input dicts as you want here
        if(not d): 
            tx_annot = None
            continue
        for key, value in d.items():
            tx_annot[key] = value

    #actually prepare the orton data (e.g. adding noise) seperately for all ecll populations(wt, braf)
    #this means each populations is scaled by its own control group (e.g. for inhibitions)
    #therefore be careful to also use the cell annotations/ data of the specific population
    rglobList = []
    rtotList = []
    if(modelWt):
        rglobWt, rtotWt, pDictWt = prepare_orton_data(datWt,'all', noise, cell_annotWt, tx_annotWt)
        rglobList.append(rglobWt)
        rtotList.append(rtotWt)
    if(modelBraf):
        rglobBraf, rtotBraf, pDictBraf = prepare_orton_data(datBraf,'all', noise, cell_annotBraf, tx_annotBraf, "braf")
        rglobList.append(rglobBraf)    
        rtotList.append(rtotBraf)    
    if(modelRas):
        rglobRas, rtotRas, pDictRas = prepare_orton_data(datRas,'all', noise, cell_annotRas, tx_annotRas, "ras")
        rglobList.append(rglobRas)
        rtotList.append(rtotRas)

    #make group annot: 'wt' ->'cellIds', , 'braf' ->'cellIds'
    group_annot = {}
    if(modelWt):
        group_annot["wt"] = list(rglobWt.columns)
    if(modelBraf):
        group_annot["braf"] = list(rglobBraf.columns)
    if(modelRas):
        group_annot["ras"] = list(rglobRas.columns)

    #combine rglob and rtot
    rglob = pd.concat(rglobList, axis=1)
    rtot = pd.concat(rtotList, axis=1) 

    return(rglob, rtot, group_annot, cell_annot, tx_annot)


def estimate_eta_for_network_complexity(etaTonodesDict, folder, name, logFile):
    lists = (etaTonodesDict.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists)

    from scipy.optimize import curve_fit

    def func(x, a, b, c):
        return a * np.exp(-b* x) + c

    etaList = []

    guess = (1.0,0.1,0.0)

    popt, pcov = curve_fit(func, np.array(y), np.array(x), p0 = guess)
    a,b,c= popt

    #plot points and function
    fig, ax = plt.subplots()
    ax.scatter(y, x)
    x = np.linspace(0,110,100)
    y = func(x, a,b,c)
    plt.plot(x,y, 'r')
    plt.title("Fitting of eta (edge penalty)")
    plt.xlabel("number of edges in network")
    plt.ylabel("eta")
    plt.show()
    plt.close()

    if(folder is not None):
        now = datetime.now()
        etaGraphObjName = folder + "/etas/" + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
        pickle.dump(ax, open(etaGraphObjName, 'wb'))

    #get actual eta guess
    etaGuess = func(13, a,b,c)
    if(etaGuess < 0):
        etaGuess = 0
    elif(etaGuess > 1):
        etaGuess = 1

    #get 22 values for network with a complexity of 1 to 22 (we have 11 nodes, so we go up to 2XnodeNUmber - as we ideally have almost all TP edges by then, resulting in 100% TPR)
    for i in range(1,23):
        result = func(i, a,b,c)
        etaList.append(result)

    return(etaGuess, etaList)

def make_mra(rglob, rtot, cell_annot, tx_annot, eta):
    if(cell_annot is None):
        scd = scmra.ScData(rglob=rglob, rtot = rtot)
    else:
        scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot)

    p = scmra.ScMraProblem(scd, eta = eta)
    p.cpx.solve()
    result = scmra.ScMraResult(p)
    return(result)

def make_cnr(rglob, rtot, cell_annot, tx_annot, group_annot, eta, theta, priorNetwork=None):
    scd = scmra.ScData(rglob=rglob, rtot=rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot=group_annot)
    p = scmra.ScCnrProblem(scd, eta=eta, theta=theta, prior_network=priorNetwork)
    p.cpx.solve()
    result = scmra.ScCnrResult(p)
    return(result)

#perfectly estimate the eta, by slowly narrowing the possible eta values down
#abort after certain number of tries and return so far best value
def get_suitable_solution_fast(rglob, rtot, cell_annot=None, tx_annot=None, group_annot = None, 
noBidirectionality = False, folder = None, logFile=None, name = None, startingEta = 0.1, 
expectedEdges = 13, etaGuesses=10, reconstructionType = "MRA", 
theta=0.0, allowedEdgesFromEstimate = 2, fitExponentialDecayIfNoSolutionFound = False):
    etaTonodesDict = {}

    returnedProblem = None # if we already found a perfect solution, do not calculate new but return
    returnedNumEdges = 0
    class Boundary:
        lower: float = 0.0
        middle: float = 0.0
        upper: float = 0.0

        def __init__(self, lower, middle, upper):
            self.lower = lower
            self.middle = middle
            self.upper = upper

    etaBoundary = Boundary(0., 0.1, 1.)

    #estimate first edges for initial eta
    if(reconstructionType == "MRA"):
        result = make_mra(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, eta = startingEta)
    elif(reconstructionType == "CNR"):
        result = make_cnr(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot=group_annot, eta = etaBoundary.middle, theta=theta)

    edges = result.imap.sum().sum()
    edgeBoundary = Boundary(0, edges, 110)

    correctEtaFound = False
    currentEtaTest = 1
    #according to expectedEdge increase, or decrease eta x-times

    solutionFound = False
    #while we did not find a solution
    while(not solutionFound):
        while(not correctEtaFound and currentEtaTest < etaGuesses):
            lastEdgeBoudary = deepcopy(edgeBoundary)
            lastEtaBoundary = deepcopy(etaBoundary)
            currentEtaTest += 1
            #we found just the perfect edge number
            if(expectedEdges == edgeBoundary.middle):
                correctEtaFound = True
                returnedProblem = result
                solutionFound = True
                a = etaBoundary.middle
                returnedNumEdges = expectedEdges
            elif(expectedEdges < edgeBoundary.middle):
                etaBoundary.lower = lastEtaBoundary.middle
                etaBoundary.middle = lastEtaBoundary.middle + ((lastEtaBoundary.upper-lastEtaBoundary.middle)/2)
                etaBoundary.upper = lastEtaBoundary.upper
                #set edge limits, new middle has to be calculated
                edgeBoundary.lower = lastEdgeBoudary.lower
                edgeBoundary.upper = lastEdgeBoudary.middle
            else:
                etaBoundary.lower = lastEtaBoundary.lower
                etaBoundary.middle = lastEtaBoundary.lower + ((lastEtaBoundary.middle-lastEtaBoundary.lower)/2)
                etaBoundary.upper = lastEtaBoundary.middle

                edgeBoundary.lower = lastEdgeBoudary.middle
                edgeBoundary.upper = lastEdgeBoudary.upper
            #otherwise continue x-times
            if(reconstructionType == "MRA"):
                result = make_mra(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, eta = etaBoundary.middle)
            elif(reconstructionType == "CNR"):
                result = make_cnr(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot=group_annot, eta = etaBoundary.middle, theta=theta)

            edges = result.imap.sum().sum()
            edgeBoundary.middle = edges

            etaTonodesDict[etaBoundary.middle] = edges

        #if we do not have sth after x tries we use the eta-edge data we generated so far
        #to fit a curve through data and check if the result is within error margin...
        #we can continue forever, until we find a reasonably good eta, then we could fit a function
        #through the previously found etas...
        if(not correctEtaFound and fitExponentialDecayIfNoSolutionFound):
            #try to get an eta
            try:
                a, b = estimate_eta_for_network_complexity(etaTonodesDict, folder, name, logFile)
            
                #if we got an eta, calculate the result
                if(reconstructionType == "MRA"):
                    result = make_mra(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, eta = startingEta)
                    exampleRloc = result.rloc
                elif(reconstructionType == "CNR"):
                    result = make_cnr(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot=group_annot, eta = etaBoundary.middle, theta=theta)
                    exampleRloc = list(result.rloc.values())[0]
                #get the number of edges from this result
                numEdges = 0
                for column_name in exampleRloc:
                    column = exampleRloc[column_name]
                    numEdges += (column != 0).sum()
                    numEdges -=1 #for the diagonal values

                if( (numEdges - expectedEdges)**2 < allowedEdgesFromEstimate**2):
                    solutionFound = True
                    returnedProblem = result
                    returnedNumEdges = numEdges
                else:
                    solutionFound = False 
                    currentEtaTest = 1 
            except:
                solutionFound = False
                currentEtaTest = 1
        else:
            #we accept whatever we found until now, we can increase the number of tries to get a better eta
            #if we do not mind longer running time
            correctEtaFound = True
            returnedProblem = result
            solutionFound = True
            a = etaBoundary.middle
            returnedNumEdges = edgeBoundary.middle
    return(a, returnedProblem, returnedNumEdges)

#the estimation metod for eta
def estimate_eta(rglob, rtot, cell_annot=None, tx_annot=None, noBidirectionality = False,  
                folder = None, logFile=None, name = None, perturbedNodeDict = None, etasInput = None,
                 modelPertRloc = False, modelPertStot=False):
    etaTonodesDict = {}
    etas = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    if(etasInput is not None):
        etas = etasInput

    for etatest in etas: 
        print("SOLVE FOR ETA: " + str(etatest)) 
        if(logFile is not None):
            f = open(logFile, "a")
            f.write("       'in eta estimation' : testEta=" + str(etatest) + " => ")
            f.close()
        if(cell_annot is None):
            scd = scmra.ScData(rglob=rglob, rtot = rtot)
        else:
            scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot)

        p = scmra.ScMraProblem(scd, eta = etatest, noBidirectionality = noBidirectionality,
                               modelPertRloc = modelPertRloc, modelPertStot = modelPertStot)
        p.cpx.solve()
        result = scmra.ScMraResult(p)
        edges = result.imap.sum().sum()
        etaTonodesDict[etatest] = edges

        if(logFile is not None):
            f = open(logFile, "a")
            f.write(str(edges) + " \n")
            f.close()

    etaTonodesDict[0.0] = 110
    etaTonodesDict[1.0] = 0

    a, b = estimate_eta_for_network_complexity(etaTonodesDict, folder, name, logFile)

    for i in range(len(b)):
        lastPosEta = 1.
        if(b[i] > 0.): 
            lastPosEta = b[i]
        elif(b[i] < 0.): 
            lastPosEta = lastPosEta/2
            b[i] = lastPosEta

    return(a, b)

# eta estimation for CNR analysis
def estimate_eta_cnr(rglob, rtot, cell_annot=None, tx_annot=None, group_annot = None,  
                folder = None, name = None, etasInput = None):
    etaTonodesDict = {}
    etas = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    if(etasInput is not None):
        etas = etasInput

    for etatest in etas: 
        print("SOLVE FOR ETA: " + str(etatest)) 
            
        scd = scmra.ScData(rglob=rglob, rtot=rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot=group_annot)

        theta = 0.0
        p = scmra.ScCnrProblem(scd, eta=etatest, theta=theta)
        p.cpx.solve()
        result = scmra.ScCnrResult(p)
        edges = result.imap.sum().sum()
        etaTonodesDict[etatest] = edges

    etaTonodesDict[0.0] = 110
    etaTonodesDict[1.0] = 0

    a, b = estimate_eta_for_network_complexity(etaTonodesDict, folder, name, None)

    return(a, b)



def calc_rmse(true_values, estimated_values):
    assert(np.array_equal(true_values.columns, estimated_values.columns))
    rmse = mean_squared_error(estimated_values.values.flatten(), true_values.values.flatten(), squared=False)
    return(rmse)

#the pickle file must be of MRASimulationClass form
#it contains a list of the results as well as an additional variable storing the parameters for those simulations
def write_metrics_to_csv_from_pickle_MRA(pickleFile, outputFile):
    #load pickle
    result = pickle.load(open(pickleFile, 'rb'))

    # Create the pandas DataFrame
    headers = result.header
    df = pd.DataFrame(columns=headers)

    #add all the existing parameters
    resultNumber = len(result.parameter)

    print(resultNumber)

    for i in range(resultNumber):
        params = result.parameter[i]
        df_length = len(df)
        df.loc[df_length] = params
        
    #extend columns by new variables (TPR, ...)
    #VARIABLES TO FILL
    tpr = []
    fpr = []
    edges = []
    rmseList = []
    precision = []
    recall = []

    for i in range(resultNumber):
        scMraResult = result.result[i]

        #GET BASIC METRICS TPR, FPR
        outputMetrics = performance_metrics(scMraResult.imap, imap_true)
        tpr.append(outputMetrics['tp']/(outputMetrics['tp'] + outputMetrics['fn']))
        fpr.append(outputMetrics['fp']/(outputMetrics['fp'] + outputMetrics['tn']))
        precision.append(outputMetrics['tp']/(outputMetrics['tp'] + outputMetrics['fp']))
        recall.append(outputMetrics['tp']/(outputMetrics['tp'] + outputMetrics['fn']))

        #GET NUMBER EDGES IN THIS NETWORK
        numEdges = 0
        for column_name in scMraResult.rloc:
            column = scMraResult.rloc[column_name]
            numEdges += (column != 0).sum()
            numEdges -=1 #for the diagonal values
        edges.append(numEdges)

        #GET RMSE (makes mostly sense for network quantification)
        scMraResult.rloc = reorder_rloc(scMraResult.rloc)
        rmse = calc_rmse(rloc_true, scMraResult.rloc)
        rmseList.append(rmse)

    #FINALLY ADD NEW COLUMNS
    df["TPR"] = tpr
    df["FPR"] = fpr
    df["EDGES"] = edges
    df["RMSE"] = rmseList
    df["PRECISION"] = precision
    df["RECALL"] = recall

    df.to_csv(outputFile, sep="\t", index=False)

#the pickle file must be of MRASimulationClass form
#it contains a list of the results as well as an additional variable storing the parameters for those simulations
#it writes an column for braf, wt, ras. ONLY these three options can have an rloc, if a simulation did not contain one of them
#the corresponding column is empty
def write_metrics_to_csv_from_pickle_CNR(pickleFile, outputFile, metricCanContainMRA = True):
    #load pickle
    result = pickle.load(open(pickleFile, 'rb'))

    # Create the pandas DataFrame
    headers = result.header
    print(headers)
    df = pd.DataFrame(columns=headers)

    #add all the existing parameters
    resultNumber = len(result.parameter)
    for i in range(resultNumber):
        params = result.parameter[i]
        df_length = len(df)
        df.loc[df_length] = params
        
    #extend columns by new variables (TPR, ...)
    #VARIABLES TO FILL
    tpr = []
    fpr = []

    rmse_dict = {'wt':[], 'braf':[], 'ras':[]}
    edges = []

    previousRepeatNumber = -1
    for i in range(resultNumber):
        scCnrResult = result.result[i]

        #check if we are really analyzing a dict for CNR data
        if(not metricCanContainMRA and not type(scCnrResult.rloc) is dict):
            ("The result of rlocs is not a dictionary mapping cell population to rlocs, are you sure this is CNR data!\n")
            exit()
            
        outputMetrics = performance_metrics(scCnrResult.imap, imap_true)
        tpr.append(outputMetrics['tp']/(outputMetrics['tp'] + outputMetrics['fn']))
        fpr.append(outputMetrics['fp']/(outputMetrics['fp'] + outputMetrics['tn']))

        #GET NUMBER EDGES IN THIS NETWORK (mut and wt matrices have same number of edges, using simply the first)
        numEdges = 0
        if(type(scCnrResult.rloc) is dict):
            wt_rloc = list(scCnrResult.rloc.values())[0] #the first value in the rloc dict
        else:
            #assume this line is MRA data
            wt_rloc = scCnrResult.rloc
        for column_name in wt_rloc:
            column = wt_rloc[column_name]
            numEdges += (column != 0).sum()
            numEdges -=1 #for the diagonal values
        edges.append(numEdges)

        #GET RMSE (makes mostly sense for network quantification)
        if(df.iloc[i,5] == "CNR"):
            foundKeys = []
            for key, value in scCnrResult.rloc.items():
                value = reorder_rloc(value)
                rmse = calc_rmse(rlocDict[key], value)
                rmse_dict[key].append(rmse)
                foundKeys.append(key)
            if("wt" not in foundKeys): rmse_dict["wt"].append("")
            if("braf" not in foundKeys): rmse_dict["braf"].append("")
            if("ras" not in foundKeys): rmse_dict["ras"].append("")

        else:
            assert(df.iloc[i,5] == "MRA")
            #assume we have MRA data
            scCnrResult.rloc = reorder_rloc(scCnrResult.rloc)
            print("MRA DETECTED____")
            print(df.iloc[i,3])
            print(previousRepeatNumber)
            print("____________")

            #get simulation type
            mraSimName = df.iloc[i,6].split(":")

            if(mraSimName[0] == "BRAF"):
                print("PREV VALUES")
                print(rmse_dict["braf"][i-1])
                print("RMSE calculation from BRAF: ",df.iloc[i,6])


                rmse = calc_rmse(rloc_braf, scCnrResult.rloc)

                rmse_dict["braf"].append(rmse)
                rmse_dict["wt"].append("NULL")
                rmse_dict["ras"].append("NULL")
            
            elif(mraSimName[0] == "RAS"):

                print("RMSE calculation from RAS: ", df.iloc[i,6])
                rmse = calc_rmse(rloc_ras, scCnrResult.rloc)

                rmse_dict["ras"].append(rmse)
                rmse_dict["wt"].append("NULL")
                rmse_dict["braf"].append("NULL")
            else:
                assert(mraSimName[0] == "WT")
                print("RMSE calculation from WT: ", df.iloc[i,6])

                rmse = calc_rmse(rloc_true, scCnrResult.rloc)
                rmse_dict["wt"].append(rmse)
                rmse_dict["braf"].append("")
                rmse_dict["ras"].append("")
            print(previousRepeatNumber)

    #FINALLY ADD NEW COLUMNS
    df["TPR"] = tpr
    df["FPR"] = fpr
    df["EDGES"] = edges
    df["RMSE_WT"] = rmse_dict["wt"]
    df["RMSE_BRAF"] = rmse_dict["braf"]
    df["RMSE_RAS"] = rmse_dict["ras"]

    df.to_csv(outputFile, sep="\t", index=False)

#simulation of orton model, with theta=0 and eta is estimated for 13 edges
def simulate_cnr(cellPopulationList, noise, cell):
    #generate the raw data
    pertIndices = None
    rglob, rtot, group_annot, cell_annot, tx_annot = generate_cnr_data(cell, pertIndices, noise, cellPopulationList)

    theta = 0.0
    scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)

    #eta must be estimated for every simulation bcs. we can have different missing proteins
    #etasInput = [0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.7, 0.9]
    eta, possibleSolution, numEdges = get_suitable_solution_fast(rglob, rtot, cell_annot, tx_annot, group_annot, reconstructionType = "CNR")
    
    return(possibleSolution, numEdges, eta)

#function to generate raw data for CNR & MRA comparison
#function generating rglob, rtot etc with noise for BOTH MRA and CNR
def generate_raw_data_for_MRA_and_CNR(cellNum, noise, mutant):
    #generate the raw data
    #sample data for both, then separately calculate rglob, rtot,...
    pertNodes = None
    missingNodes = None
    pertubationIndices = get_index_list(allNodes, pertNodes)
    dat, cell_annot, tx_annot = sample_raw_mra_data(wt, pertubationIndices, cellNum, PerturbDataList)
    if(mutant == "braf"):
        datB, cell_annotM, tx_annotM = sample_raw_mra_data(braf, pertubationIndices, cellNum, PerturbDataList)
    elif(mutant == "ras"):
        datB, cell_annotM, tx_annotM = sample_raw_mra_data(ras, pertubationIndices, cellNum, PerturbDataList)

    sampledData = {"wt": [dat, cell_annot, tx_annot], mutant : [datB, cell_annotM, tx_annotM]}

    #generate rglob, rloc
    wtData = sampledData['wt']
    datWt = wtData[0]
    cell_annot = wtData[1]
    tx_annot = wtData[2]

    mutData = sampledData[mutant]
    datMut = mutData[0]
    cell_annotM = mutData[1]
    tx_annotM = mutData[2]

    #we implement this test without any perturbations
    #we r only interested in how MRA compares to CNR
    #therefore we make our lifes easier and do not handle these two objects
    assert(cell_annot is None)
    assert(tx_annot is None)
    assert(cell_annotM is None)
    assert(tx_annotM is None)

    #actually prepare the orton data (e.g. adding noise) seperately for all ecll populations(wt, braf)
    #this means each populations is scaled by its own control group (e.g. for inhibitions)
    rglobWT, rtotWT, pDictWt = prepare_orton_data(datWt,'all', noise)
    rglobMUT, rtotMUT, pDictBraf = prepare_orton_data(datMut,'all', noise, None, None, mutant)

    return(rglobWT, rtotWT, rglobMUT, rtotMUT)


#function to simulate a simple MRA with a certain population
def simulate_mra_from_data(cellPopulation, noise, cell, f, count, result, rglob, rtot, cell_annot, tx_annot, PopulationID):

    #etasInput = [0.005, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.5, 0.7, 0.9]
    #eta, etaListNew = estimate_eta(rglob, rtot, cell_annot, tx_annot, withFlattening = True, etasInput=etasInput)
    edgesToReconstruct = 13
    if(cellPopulation == "ras"):
        edgesToReconstruct = 12
    eta, scMraResult, edgeNumber = get_suitable_solution_fast(rglob, rtot, cell_annot, tx_annot, etaGuesses=30, expectedEdges = edgesToReconstruct)
    
    #update logfile
    f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(count) + 
            ", etas: " + str(round(eta,3)) )

    #add simulation output to result
    result.add_result(scMraResult)

    #add values to result
    numDiffEdges = 0
    residuals = 0
    n_res = np.size(scMraResult.residuals_complete) + \
            np.size(scMraResult.residuals_incomplete)
    residuals += np.sum(np.array(np.square(scMraResult.residuals_complete)))/n_res
    residuals += np.sum(np.array(np.square(scMraResult.residuals_incomplete)))/n_res
    result.add_parameter([noise, cell, eta, count, edgeNumber, "MRA", PopulationID, numDiffEdges, "NO_THETA",residuals])

    #update logfile
    f.write(str(round(eta,3)) + " edges: " + str(edgeNumber) + "\n")
    f.flush()

#function to simulate a CNR with a defined set of populations among the possible ones=[wt, ras, braf]
def simulate_cnr_from_data(cellPopulationList, noise, cell, f, count, result, rglob, rtot, cell_annot, tx_annot, group_annot, theta):
    
    #theta = THETA
    #scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)

    eta, scCnrResult, edgeNumber = get_suitable_solution_fast(rglob, rtot, cell_annot, tx_annot, group_annot, 
                etaGuesses=30, reconstructionType = "CNR", theta = theta)

    #get number of different edges
    vars_lst = [var for var in scCnrResult.vardict if var.startswith('IDev')]
    vars_vals = [scCnrResult.vardict[v] for v in vars_lst]
    lengthIdcs = (np.count_nonzero(vars_vals))
    numDiffEdges = lengthIdcs

    #update logfile
    f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(count) + 
            ", etas: " + str(round(eta,3)) + " diff Edges: " + str(lengthIdcs) )
    
    s = scCnrResult

    populationString = ""
    for population in cellPopulationList:
        populationString += population
        populationString += ','
    populationString = populationString[:-1]

    residuals = 0
    n_res = np.size(scCnrResult.residuals_complete) + \
            np.size(scCnrResult.residuals_incomplete)
    residuals += np.sum(np.array(np.square(scCnrResult.residuals_complete)))/n_res
    residuals += np.sum(np.array(np.square(scCnrResult.residuals_incomplete)))/n_res

    result.add_parameter([noise, cell, eta, count, edgeNumber, "CNR", populationString, numDiffEdges, theta, residuals])

    f.write(" edges: " + str(edgeNumber) + "\n")
    f.flush()

    #add simulation output to result
    result.add_result(s)