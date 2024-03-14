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
import umap
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
#distribution of distances in cell space
from sklearn.metrics import pairwise_distances
import getopt 

import random

'''
Run the generation of permutations to check significance of edges:
using theta of 0.1: an over view is in analysis.ipynb. Chosen to simply find around 10 diff edges
'''

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

def shuffle_group_label_labels(group_annot):
    keys = list(group_annot.keys())
    keyValueNumber = {}
    valueList= []
    for key in keys:
        keyValueNumber[key] = len(group_annot[key])
        valueList.append(group_annot[key])

    valueList = [item for sublist in valueList for item in sublist]

    group_annot_return = {}
    group_annot_return[keys[0]] = random.sample(valueList, keyValueNumber[keys[0]])
    group_annot_return[keys[1]] = random.sample(list(set(valueList) - set(group_annot_return[keys[0]])), keyValueNumber[keys[1]])
    
    if(len(keys) >= 3):
        listSoFar = list(set(valueList) - set(group_annot_return[keys[0]]))
        listSoFar = list(set(listSoFar) - set(group_annot_return[keys[1]]))
        group_annot_return[keys[2]] = random.sample(listSoFar, keyValueNumber[keys[2]])

    if(len(keys) == 4):
        listSoFar = list(set(listSoFar) - set(group_annot_return[keys[2]]))
        group_annot_return[keys[3]] = random.sample(listSoFar, keyValueNumber[keys[3]])

    return(group_annot_return)

def get_pairwise_distances(data, method = 'cosine', seperate_z_scale = False):
    def upper_tri_indexing(A):
        m = A.shape[0]
        r,c = np.triu_indices(m,1)
        return A[r,c]


    # Z-NORMALIZE DATA
    scaledForGroup = None
    if(seperate_z_scale):
        for group in np.unique(data['treatment']):
            data_values = data[data['treatment'] == group].drop(['sample_id', 'treatment'], axis=1).values
            scaledForGroupAddon = StandardScaler().fit_transform(data_values)
            if(scaledForGroup is None):
                scaledForGroup = scaledForGroupAddon
            else:
                scaledForGroup = np.concatenate((scaledForGroup, scaledForGroupAddon), axis=0)
        data = scaledForGroup
    else:
        data_values = data.drop(['sample_id', 'treatment'], axis=1).values
        data = StandardScaler().fit_transform(data_values)

    #EGF_DATA ROW=cell columns=FEATURES (cells in rows)
    #pairwiseDistances diagonal is zero, upper triangular matrix: row=from, col=to (dsitance of cell in row 0 to cell 9 is pairwise[0,9])
    pairwiseDistances = pairwise_distances(data, metric=method)
    
    #list: only to visualize distances (other than that not an interesting data structure)
    disList = upper_tri_indexing(pairwiseDistances)

    import matplotlib.pyplot as plt
    plt.hist(disList, bins=100)
    plt.show()

    return(pairwiseDistances)



def map_string_to_ints(stringList):
    mydict={}
    i = 0
    for item in stringList:
        if(i>0 and item in mydict):
            continue
        else:    
            i = i+1
            mydict[item] = i
    k=[]
    for item in stringList:
        k.append(mydict[item])
    return(k)

def plot_graph(data, subset = None, removeIsolates = True, threshold = 0.65, method = 'cosine', protein = 'JNK-p', treatmentColor = False,
                rescaleToGroup = False, labelBYRank = False, seperate_z_scale = False, subset_data = None, edgeMap = False):

    localData = data.copy()

    adjMat = get_pairwise_distances(localData, method = method, seperate_z_scale = seperate_z_scale)
    adjMat[adjMat > threshold] = 0

    G = nx.from_numpy_array(adjMat)
    G.edges(data=True)

    #remove isolate nodes, they r not useful for graph matching
    if(removeIsolates):
        G.remove_nodes_from(list(nx.isolates(G)))
        #components smaller than...
        theshold_cluster_size = 3
        for component in list(nx.connected_components(G)):
            if len(component)<theshold_cluster_size:
                for node in component:
                    G.remove_node(node)

    #get colour scheme
    plotData = localData.copy()

    plotData = plotData.reset_index()
    index = plotData.index

    color_map = []
    if(treatmentColor):
        for node in G:
            color_map.append(plotData['treatment'][node])
        color_map = map_string_to_ints(color_map)
        minVal = 1
        maxVal = len(np.unique(color_map))
    else:
        ranksPerGroup = {}
        if(rescaleToGroup):
            for group in subset:
                ranksPerGroup[group] = ss.rankdata(plotData[plotData['treatment'] == group][protein])
                plotData.loc[plotData['cluster_id'] == group, 'rank'] = ranksPerGroup[group]
            for node in G:
                color_map.append(plotData['rank'][node]) #bcs ranks are 1-based indexed
        else:  
            for node in G:
                color_map.append(plotData[protein][node])
                #MIN and MAX values are values in range of WHOLE dataset and not only subset to see differences
        minVal = np.min(plotData[protein])
        maxVal = np.max(plotData[protein])

    if(edgeMap):
        edge_map = []
        for node in G:
            edge_map.append(len(G.edges(node)))
        minVal = np.min(edge_map)
        maxVal = np.max(edge_map)
        meanEdge = np.mean(edge_map)
        removeNodes = []
        for node in G:
            if(len(G.edges(node)) < meanEdge/4):
                removeNodes.append(node)
        #remove node and calculte edge map new
        for node in removeNodes:
            G.remove_node(node)
        edge_map = []
        for node in G:
            edge_map.append(len(G.edges(node)))
            
        color_map = edge_map

    #treamtnt col overwriting edge color, so we can reduce edges but color by treatment :D 
    if(treatmentColor):
        color_map = []
        for node in G:
            color_map.append(plotData['treatment'][node])
        color_map = map_string_to_ints(color_map)
        minVal = 1
        maxVal = len(np.unique(color_map))

    labeldict = {}
    if(labelBYRank):
        for group in subset:
            ranksPerGroup = ss.rankdata(plotData[plotData['cluster_id'] == group][protein])
            plotData.loc[plotData['cluster_id'] == group, 'rank'] = ranksPerGroup

        keys = []
        vals = []
        for node in G:
            val = plotData['rank'][node]
            strVal = ""
            #if(val < 10 or val > 150):
            if(val):
                strVal = str(val)
            vals.append(strVal)
            keys.append(node)
        labeldict = dict(zip(keys,vals))


    #LRP6_P
    #NFKB_P65_P
    #RIBOSOMAL_S6_P
    fig,ax = plt.subplots(figsize=(15, 10), dpi=80)
    pos = nx.spring_layout(G, seed=225)  # Seed for reproducible layout
    nx.draw(G, pos, with_labels=True, labels=labeldict, node_size = 350, alpha = 0.7, cmap='viridis', node_color=color_map, vmin=minVal, vmax=maxVal)
    #nx.draw(G, pos, with_labels=False, node_size = 150, alpha = 0.5, cmap=plt.cm.Reds, node_color = edge_map)

    plt.show(fig)

    return G

def plot_3d(xs, ys, zs):  
    fig,ax = plt.subplots(figsize=(15, 10), dpi=80)
    ax = fig.add_subplot(projection='3d')
    n = 100
    ax.scatter(xs, ys, zs, s=80)

    ax.set_xlabel('DIM 1')
    ax.set_ylabel('DIM 2')
    ax.set_zlabel('DIM 3')

    plt.show()

def plot_net(dfData, threshold=0.6):
    import igraph as ig
    import leidenalg as la

    adjMat = get_pairwise_distances(dfData, method = "cosine")

    adjMat[adjMat > threshold] = 0
    adjMat[adjMat > 0] = 1

    adjMat = np.asmatrix(adjMat)
    G = ig.Graph.Adjacency(adjMat)
    #make an edge list with weights for Leiden clustering
    partition = la.find_partition(G, la.ModularityVertexPartition,n_iterations=20)#, weights = weightVec)

    #label the cluster name
    partitionVec = []
    clusterDict = {}
    pNum = (0)
    for p in partition:
        if(len(p) < 5): continue
        curList = []
        for node in p:
            curList.append(str(node))
            clusterDict[node] = str(pNum)
        partitionVec.append((curList))
        pNum += 1
    tmp = dfData.copy()
    tmp.reset_index(inplace=True)
    tmp = tmp.rename(columns = {'index':'ab_name'})
    x = tmp["ab_name"].map(clusterDict)

    ig.plot(partition, vertex_label = x, seed = 1)



PERMUTATIONS = 1000
THETA = GLOBALTHETA_EGFRINH

#returns a list
def run_permutation_analysis_with_fixed_differences(egfrInhDataProcessed, masterCnrResult):
    #clusterList = ["vehicle", "vehiclepRB", "AG1478"]
    clusterList = ["vehicle", "AG1478"]

    rglob, rtot, cell_annot, tx_annot, group_annot = prepare_data_for_EGFR_inhibition(egfrInhDataProcessed)
    edgeStrengthList = [] #list of DFs, where every df stores the edge differences for different clusters

    for i in range(PERMUTATIONS):

        #shuffle the labels...
        group_annot_shuffled = shuffle_group_label_labels(group_annot)
        scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot_shuffled)
        #solve the new problem
        scp = scmra.ScCnrProblem(scd, eta=0.0, theta=THETA, prior_network = CANNONICAL_EDGES) 

        #matrix for every simulation: storing in row-col (cluster-edge) the edge strength in this cluster
        #this is done for every edge for now; later we subset the ones with status == 1 (diff is set)
        deviatingFromToList = []
        for indicator, status in masterCnrResult.vardict.items():
            if indicator.startswith(("IDev", "ISDev", "IPDev")):
                #we need to also set contraints to ENFORCE this edge in networks reconstruction
                scmra.cplexutils.set_indicator_status(scp.cpx, indicator, status)
                if(indicator.startswith("IDev")):
                    indSplit = indicator.split("_")
                    deviatingFromToList.append([indSplit[1],indSplit[2]])
        rows = clusterList
        cols = [x[0] + "_" + x[1] for x in deviatingFromToList]
        clusterEdgeMatrix = pd.DataFrame(columns = cols, index = rows)

        scp.cpx.solve()
        s = scmra.ScCnrResult(scp)

        #make a matrix storing edge, weight, cluster
        for clusterID in clusterList:
            for edge in deviatingFromToList:
                clusterEdgeMatrix.loc[clusterID, edge[0] + "_" + edge[1]] = s.rloc[clusterID].loc[edge[0],edge[1]]
        edgeStrengthList.append(clusterEdgeMatrix)

    return(edgeStrengthList)

#write two dataframes: one with the true edge differences (between each population and the vehicle)
#and the calcualted edge differences from the label-permuted reconstructions
def write_edge_subpopulation_differences(trueEdges, resultPermSimulations, runOnBootrappedEdges):
    #WRITE ASLSO THE SEPERATE ANALYSES INTO A TSV FILE
    #get all different edges
    colsToKeep = []
    for indicator, status in trueEdges.vardict.items():
        if(indicator.startswith(("IDev")) and status==1):
            indSplit = indicator.split("_")
            colsToKeep.append(indSplit[1] + "_" + indSplit[2])
    
    #subset these column now
    resultEdgePerm = []
    for i in resultPermSimulations:
        df = i[colsToKeep]
        resultEdgePerm.append(df)

    #WRITE THE EDGE DIFFERENCES BETWEEN ITGB1low/PROLIFERATING AND VEHICLE INTO A TSV
    cols = ["GROUP","EDGE","VALUE"]
    dfEdgeDiff = pd.DataFrame(columns = cols)
    for i in resultEdgePerm:
        colName = i.columns
        rowNames = i.index
        for col in colName:
            for row in rowNames:
                if(not row=="vehicle"):
                    arr = [row, col, i.loc["vehicle", col] - i.loc[row, col]]
                    dfEdgeDiff = pd.concat([dfEdgeDiff, pd.DataFrame([arr], columns=["GROUP","EDGE","VALUE"])], axis=0)

    if(runOnBootrappedEdges):
        dfEdgeDiff.to_csv("../bin/edgeSignificanceEGFRINH_RecoveredDifferences_BootrapEdges.tsv", sep="\t")
    else:
        dfEdgeDiff.to_csv("../bin/edgeSignificanceEGFRINH_RecoveredDifferences.tsv", sep="\t")

    #WRITE THE EDGE DIFFERENCES FORM MASTERCNR TO A TSV
    cols = ["GROUP","EDGE","VALUE"]
    truePoint = pd.DataFrame(columns = cols)
    populations = list(trueEdges.rloc.keys())
    populations.remove("vehicle")
    for edge in colsToKeep:
        for pop in populations:
            toFromList = edge.split("_")
            val = trueEdges.rloc["vehicle"].loc[toFromList[0], toFromList[1]] - trueEdges.rloc[pop].loc[toFromList[0], toFromList[1]]
            arr = [pop, edge, val]
            truePoint = pd.concat([truePoint, pd.DataFrame([arr], columns=["GROUP","EDGE","VALUE"])], axis=0)

    if(runOnBootrappedEdges):
        truePoint.to_csv("../bin/edgeSignificanceEGFRINH_TrueDifferences_BootrapEdges.tsv", sep="\t")
    else:
        truePoint.to_csv("../bin/edgeSignificanceEGFRINH_TrueDifferences.tsv", sep="\t")


def get_edges_from_bootstrap_data():
    edgesBootstrap = pd.read_csv("../bin/bootstrap_RecoveredDifferences_EGFR.tsv",sep='\t')
    edges = edgesBootstrap['EDGE'].value_counts().head(12)
    edges = edges.index.tolist()

    print(edges)

    return(edges)

def make_bootstrap_master_cnr_result(egfrInhDataProcessed, edges):
    rglob, rtot, cell_annot, tx_annot, group_annot = prepare_data_for_EGFR_inhibition(egfrInhDataProcessed)

    scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)
    #solve the new problem
    scp = scmra.ScCnrProblem(scd, eta=0.0, theta=THETA, prior_network = CANNONICAL_EDGES) 

    exampleResult = open("../bin/MASTERCNR_EGFRINH_CNRRESULT.pickle", "rb")
    exampleResult = pickle.load(exampleResult)
    for indicator, status in exampleResult.vardict.items():
        if indicator.startswith(("IDev")):
            indSplit = indicator.split("_")
            currentEdge = indSplit[1] + "_" + indSplit[2]
            if(currentEdge in edges):
                scmra.cplexutils.set_indicator_status(scp.cpx, indicator, 1)
            else:
                scmra.cplexutils.set_indicator_status(scp.cpx, indicator, 0)

    scp.cpx.solve()
    s = scmra.ScCnrResult(scp)

    return(s)

def main(argv):
    runOnBootrappedEdges = False
    egfrInhDataProcessed = pd.read_csv("../bin/EGFRInhibitionDataProcessed.tsv",sep='\t')

    try:
        opts, args = getopt.getopt(argv,"hb")

    except getopt.GetoptError:
        print('scMRAResult_to_tsv.py -b <if we run with bootrstrapped edges>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('scMRAResult_to_tsv.py -b <if we run with bootrstrapped edges>')
            sys.exit()
        elif opt == '-b':
            runOnBootrappedEdges = True

    if(runOnBootrappedEdges):
        edges = get_edges_from_bootstrap_data()
        print(edges)
        masterCnrResult = make_bootstrap_master_cnr_result(egfrInhDataProcessed, edges)
    else:
        objectMasterResult = open("../bin/MASTERCNR_EGFRINH_CNRRESULT.pickle", "rb")
        masterCnrResult = pickle.load(objectMasterResult)

    resultPermSimulations = run_permutation_analysis_with_fixed_differences(egfrInhDataProcessed, masterCnrResult)
    write_edge_subpopulation_differences(masterCnrResult, resultPermSimulations, runOnBootrappedEdges)


if __name__ == "__main__":
   main(sys.argv[1:])