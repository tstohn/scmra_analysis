
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

def estimate_eta_list(cell, noise):
    rglob, rtot, pDict = prepare_orton_data(wt, cell, noise)
    eta, etaList = estimate_eta(rglob, rtot)
    return(etaList)

#more sophisitcated evaluation how food rloc/stot plus new scaling actually works (we have some values from previous tests, just do some more now...)
noiseArray = [0., 0.1, 0.2, 0.5]
cellNumArray = [1000]
etaArray = [0.001, 0.01, 0.025, 0.05, 0.075,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
repeats = 10

#add header to result
header = ['NOISE', 'CELLNUM', 'ETA', 'REPEAT', 'EDGES']
result = MRASimulationClass(header)
logFile = "../bin/logCnrDataGeneration.txt"
f = open(logFile, "w+")

for noise in noiseArray:
    for cell in cellNumArray:
        count = 0
        for eta in etaArray:
            for repeat in range(repeats):

                count += 1

                #generate the raw data
                rglob_wt, rtot_wt     = prepare_orton_data(wt, cell, noise)
                rglob_braf, rtot_braf = prepare_orton_data(braf, cell, noise)
                rglob = pd.concat([rglob_wt, rglob_braf], axis=1)
                rtot = pd.concat([rtot_wt, rtot_braf], axis=1) 
                group_annot = {
                    "wt": list(rglob_wt.columns),
                    "braf": list(rglob_braf.columns) 
                    }
                cell_annot = None
                # cell_annot = {
                #     'ctr': list(rglob_wt.columns) + list(rglob_braf.columns), 
                #     'rasi': list(rglob_wt_rasi.columns) + list(rglob_braf_rasi.columns), 
                #     'p90i': list(rglob_wt_p90i.columns) + list(rglob_braf_p90i.columns)
                # }
                tx_annot = None
                # tx_annot = {'ctr': None, 'rasi': 'Ras', 'p90i': 'P90Rsk'}

                theta = 0.0
                scd = scmra.ScData(rglob=rglob, rtot = rtot, cell_annot=cell_annot, tx_annot=tx_annot, group_annot= group_annot)
                
                #update logfile
                f.write("SIMULATION: noise=" + str(noise) + ", cells=" + str(cell) + ", repeat=" + str(count) + " etas: " + str(round(eta,3)))

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
                result.add_parameter([noise, cell, eta, count, numEdges])

                #add simulation output to result
                result.add_result(s)
                del(p) #delete bcs. cplex complaints for some reason (probably still accessing an old instance)
                f.write("\n")
                f.flush()
f.close()

#store pickled result
now = datetime.now()
folder = "../bin/"
name = "SIMULATION_CNR"
pickleFile = folder + name + "_" + now.strftime("%m-%d-%Y|%H:%M:%S") + '.pickle'
pickle.dump(result, open(pickleFile, 'wb'))