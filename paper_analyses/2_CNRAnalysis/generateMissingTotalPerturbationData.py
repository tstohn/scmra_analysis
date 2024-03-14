
#generate data do compare network quality
import sys
sys.path.insert(0, '../../helperScripts/')
sys.path.insert(0,'/home/t.stohn/Tools/CPLEX_Studio201/cplex/python/3.8/x86-64_linux')

from helperFunctions import *

def prepare_orton_data(dat, n_cells='all', noise=0, cell_annot=None):
    if n_cells == 'all':
        dat_sample = dat.copy(deep=True)
    else:
        dat_sample = dat.copy(deep=True).sample(n_cells)
  
    if cell_annot is None:
        dat_scaled = ((dat_sample - dat_sample.median())/dat_sample.median()).transpose()
    else:
        ctr_cells = [c for c, tx in cell_annot.items() if tx == "CTR"]
        dat_ctr = dat_sample.loc[ctr_cells]
        dat_scaled = ((dat_sample - dat_ctr.median())/dat_ctr.median()).transpose()
        
    dat_noise = add_noise(dat_scaled, noise)
    
    rglob = dat_noise.loc[[n for n in list(dat_noise.index) if ("Active" in n)]]

    rglob.index = [i.replace("Active", "") for i in rglob.index]
    rtot = dat_noise.loc[[n for n in list(dat_noise.index) if ("Tot" in n)]]
    rtot.index = [i.replace("Tot", "") for i in rtot.index]
    return rglob, rtot

################################################
# CNR SIMULATIONS
################################################

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0., 0.2]
cellNumArray = [1000]
repeats = 10
missingTotal = [0,5,11]
perturbations = [0,5,11]

#add header to result
header = ['NOISE', 'CELLNUM', 'ETA', 'REPEAT', 'MISSINGTOTAL', 'EDGES', 'MISSINGSTRING','PERTURBATIONS','PERTURBATIONSTRING']
result = MRASimulationClass(header)
logFile = "../bin/logCnrMissingTotalPerturbationsShallowInh.txt"
f = open(logFile, "w+")

for noise in noiseArray:
    for cell in cellNumArray:
        for missingNodesVal in missingTotal:
            for pert in perturbations:
                for repeat in range(repeats):

                    #generate the raw data
                    pertNodes = random.sample(allNodes, pert)
                    pertubationIndices = get_index_list(allNodes, pertNodes)
                    rglob, rtot, group_annot, cell_annot, tx_annot = generate_cnr_data(cell, pertubationIndices, noise)

                    #eta must be estimated for every simulation bcs. we can have different missing proteins
                    #for runtime reasons estimate it here
                    if(repeat == 0):
                        etasInput = [0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7]
                        eta, etaListNew = estimate_eta_cnr(rglob, rtot, cell_annot, tx_annot, group_annot, withFlattening = True, etasInput=etasInput)

                    #delete missing nodes
                    missingNodes = random.sample(allNodes, missingNodesVal)
                    if(missingNodes is not None):
                        rtot = rtot.drop(missingNodes, axis=0)

                    theta = 0.0
                    scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)

                    #update logfile
                    f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(repeat) + 
                            ", etas: " + str(round(eta,3)) + ", " + " missing: " + str(missingNodesVal) + ", perturbed: " + str(pert))
                    
                    #make MRA simulation
                    p = scmra.ScCnrProblem(scd, eta=eta, theta=theta) 
                    p.cpx.solve()
                    s = scmra.ScCnrResult(p)

                    #add values to result
                    #GET NUMBER EDGES IN THIS NETWORK
                    numEdges = 0
                    exampleRloc = s.rloc["wt"]
                    for column_name in exampleRloc:
                        column = exampleRloc[column_name]
                        numEdges += (column != 0).sum()
                        numEdges -=1 #for the diagonal values
                    joined_missingTotal_string = ",".join(missingNodes)
                    joined_pert_string = ",".join(pertNodes)

                    result.add_parameter([noise, cell, eta, repeat, missingNodesVal, numEdges, joined_missingTotal_string, pert, joined_pert_string])

                    f.write(" edges: " + str(numEdges) + "\n")
                    f.flush()

                    #add simulation output to result
                    result.add_result(s)
                    del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)
f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_CNR_MISSINGTOTAL_WITH_PERTURBATIONS_ShallowInh"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))