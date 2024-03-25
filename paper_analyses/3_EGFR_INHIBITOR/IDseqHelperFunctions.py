from helperFunctions import *

#global analysis variables: 8 different edges
#GLOBALTHETA_EGFRINH = 0.08
#GLOBALTHETA_EPIDIFF = 0.15 #FOR 3 SUBGROUPS

#for treated vs. untreated
GLOBALTHETA_EPIDIFF = 0.027

#for 12 different edges with 3 groups: 0.053
#with four groups
GLOBALTHETA_EGFRINH_4groups = 0.094
GLOBALTHETA_EGFRINH = 0.0075

#definitions of proteins for data preparation
#problematic  EGFR-p, GDK3b=? GSK3b, stat1_p_p91 :? stat1_p
ab_mapping_phospho = {
"RSK-p" : "RSK",
'cFos-p' : "CFOS",
'cJun-p' : "CJUN",
'Fak-p' : "FAK",
#EDFRY1045 is not 100% known, but guessed from old panel
#"EGFR-p" : "EGFRY1045",  "EGFR-p_AF1095" : "EGFRY1173" ,
"EGFR-p" : "EGFRY1045",  "EGFR-p_AF1095" : "EGFRY1173" ,
"ERK1-2-p" : "ERK12",
"RSK-p90-p" : "p90RSK",
"Akt-p_4060" : "AKT1", "Akt123-p" : "AKT123",
"mTOR-p": "MTOR",
"Ribosomal-S6-p" : "RPS6",
"ITGB1" : "ITGB1",
"JNK-p" : "JNK",
"Rb-p" : "RB",
"Cdc2" : "CDC2",
"cMyc" : "CMYC",
'IKB-a-p' : "IKBA",
'H2A-p' : "H2A",
'H3-p' : "H3",
"CDK4": "CDK4",
"BMPRII": "BMPR",
'CREB-p' : "CREB",


"MAPK-p38-p2" : "P38D",
"MKK3-6-p" : "MKK36",
"NFKB-p65-p" :  "P65",
"GDK3b-p" : "GSK3B",

"MAPK-APK2-p" : "MAPKAP2",
"MAPK-p38-p" : "P38",
"BMP2_4" : "BMP",

"Jak1-p" : "JAK1",

"STAT1-p" : "STAT1",
"STAT1-p_p91" : "STAT1",
"STAT3-p" : "STAT3",
"STAT5-p" : "STAT5",

"Frizzled-3" : "FRZ",
"LRP6-p" : "LRP6",

"Src-p" : "SRC",
"smad1-5-9-p" : "SMAD159",
"smad1-5-p" : "SMAD15",

"Cyclin-B1" : "CYCLINB1",
"Cyclin-E" : "CYCLINE",

"RelA-p" : "RELA",

######### NEW PROTEINS ######
"PLC-g1-p" : "PLCG1",
"PLC-g2-p" : "PLCG2",
'Ephrin-B3' : "EPHRINB3",
'IGFBP-3' : "IGFBP3"
}

#should I duplocate total AKT for both phospho, but same total???...
#mapping of: {TOTAL -> PHOSPHO(new name as assigned above)}
ab_mapping_tot = {
#'GAPDH', 
# #'Notch-1'   'Notch-1' 'Notch-1-cleaved' 'Notch-2' 'Notch-3'
#'smad2-3' 'smad3',
# 'p63'
# 'ITGA6'
#   'Jagged-1'
#  'KLF4' 
# 'Kallikrein-6'  
# 'RNAPOLII'
#    'TGM1' 

"ERK1-2" : "ERK12",
"Akt2": ["AKT1", "AKT123"], 
"JNK" : "JNK",
"Fak" : "FAK",
"Jak1" : "JAK1",
"Src" : "SRC",
"RelA" : "RELA",
"smad1" : ["SMAD159", "SMAD15"],
'IKB-a': "IKBA",

############# NEW PROTEINS #########

}
ab_mapping_tot_inv = {
"ERK12":"ERK1-2",
"AKT1" : "Akt2",
"AKT123" : "Akt2",

"JNK" : "JNK",
"FAK" : "Fak",
"JAK1" : "Jak1",
"SRC" : "Src",
"RELA" : "RelA",
'IKBA': "IKB-a",

"SMAD159" : "smad1",
"SMAD15" : "smad1"
}

ab_mapping = {**ab_mapping_tot, **ab_mapping_phospho}
ab_use = list(ab_mapping_tot.keys()) + list(ab_mapping_phospho.keys())
ab_use_phospho = list(ab_mapping_phospho.keys())
ab_use_tot = list(ab_mapping_tot.keys())

def flatten(A):
    rt = []
    for i in A:
        if isinstance(i,list) or isinstance(i, tuple): rt.extend(flatten(i))
        else: rt.append(i)
    return rt

#MAPPING OF PROTEINS
#for orton data we had Active/ Tot at AB-name as suffix, but it is removed in the final step (so we can directly name our protein same for tot/ active)

NODE_NAMES = set(flatten(ab_mapping.values()))

OUTPUT_NODES = ["CMYC", "CFOS", "CYCLINE", "CYCLINB1", "H2A", "H3"]
OUTPUT_EDGES = []
for out_node in OUTPUT_NODES:
    OUTPUT_EDGES += [(node, out_node) for node in NODE_NAMES if node is not out_node]
OUTPUT_EDGES

INPUT_NODES = ["ITGB1", "BMP"]
INPUT_EDGES = []
for in_node in INPUT_NODES:
    INPUT_EDGES += [(in_node, node) for node in NODE_NAMES if node is not in_node]

ALL_EDGES = list(itertools.permutations(NODE_NAMES, 2))

set(flatten(OUTPUT_EDGES) + flatten(INPUT_EDGES)) - set(NODE_NAMES)
set(flatten(ALL_EDGES)) - set(NODE_NAMES)

#REMOVED FZD edges
CANNONICAL_EDGES = [
    # MAPK pathway
    ("ERK12", "EGFRY1173"), ("ERK12","EGFRY1045"),
    ("AKT1", "EGFRY1173"), ("AKT1", "EGFRY1045"),
    ("AKT123", "EGFRY1173"), ("AKT123", "EGFRY1045"), 
    ("JNK", "EGFRY1173"), ("JNK", "EGFRY1045"), 
    ("SRC", "EGFRY1173"), ("SRC", "EGFRY1045"), 

    #("ERK12", "EGFR"),
    #("AKT1", "EGFR"),
    #("AKT123", "EGFR"),
    #("JNK", "EGFR"),
    #("SRC", "EGFR"),

    ("RSK", "ERK12"),("CMYC", "ERK12"), ("CFOS", "ERK12"), #("RSK1S380", "ERK12"),
    ("STAT3", "ERK12"), # ERK
    ("RPS6", "RSK"), # RSK("RPS6", "RSK1S380"), 
    # AKT pathway
    ("MTOR", "AKT1"), ("MTOR", "AKT123"), ("GSK3B", "AKT1"), ("GSK3B", "AKT123"),
    ("IKBA", "AKT1"), ("IKBA", "AKT123"), # AKT
    ("RPS6", "MTOR"), # MTOR
    # Other
    ("CJUN", "JNK"), # JNK
    ("SRC", "ITGB1"), # ITGB1
    ("FAK", "SRC"), # SRC
    ("STAT3", "JAK1"), ("STAT5", "JAK1"), ("STAT1", "JAK1"), # JAK
#    ("FZD", "LRP6"), # LRP6
#    ("GSK3B", "FZD"), # FZD
    # MKK signaling
    ("MKK36", "GSK3B"), #GSK3B
    ("P38", "MKK36"),("P38D", "MKK36"), # MKK36
    ("MAPKAP2", "P38"), ("MAPKAP2", "P38D"), 
    ("CREB", "P38"), ("CREB", "P38D"), # P38
    ("P65", "IKBA"), ("RELA", "IKBA"), # IKBA
    ("BMPR", "BMP"), # BMP24
    ("SMAD159", "BMPR"), ("SMAD15", "BMPR"),
    
    ("CYCLINE", "CDC2"), ("CYCLINE", "GSK3B") # These phosphorylate CCNE1 --> Protein stability
]

# add edges of pps
psp = pd.read_csv("/DATA/t.stohn/MRA/scMRA-analysis/data/resources/substrate-kinase-psp.tsv", sep="\t")
psp = psp[psp['SUB_NODE_NAME'].isin(NODE_NAMES)]
psp = psp[psp['KIN_NODE_NAME'].isin(NODE_NAMES)]
PSP_INTERACTIONS = list(zip(psp.SUB_NODE_NAME, psp.KIN_NODE_NAME))

#some tests
CANNONICAL_EDGES = list(set(PSP_INTERACTIONS).union(set(CANNONICAL_EDGES)))
assert set(flatten(CANNONICAL_EDGES)).issubset(set(NODE_NAMES))
set(flatten(CANNONICAL_EDGES)) ^ set(NODE_NAMES)

#make CNR
def prepare_data_for_epidermal_differentiation(dat):
    #data table
    #sample_id  treatment   ab_count_tmm    ab_name
    dat.sample_id = [s.replace("_", ".") for s in dat.sample_id]
    dat = dat.set_index('sample_id')

    group_annot = dict()
    for group in dat.treatment.unique():
        group_annot[str(group)] = np.unique(list(dat[dat.treatment == group].index))

    #for completeness how cells and treatment would be annotated
    #HOWEVER: we assign NO TREATMENTS and handle the three populations as independant groups
    cell_annot = None
    tx_annot = None
    '''
    cell_annot = dict()
    for tx in dat.treatment.unique():
        cell_annot[tx] = np.unique(list(dat[dat.treatment == tx].index))
    tx_annot = 
    {
        "differentiated(ITGB1low)" : None, 
        "proliferating" : None, 
        "vehicle" : None
    }
    '''

    #scale every group by its own median
    #subset data by group
    dat_scaled = None
    for group in dat.treatment.unique():
        datTmp = dat.loc[group_annot[group]]

        datTmp = datTmp.pivot(columns='ab_name', values='ab_count_tmm')
        datTmp.reset_index()

        ctr_median = datTmp.median()
        datScaledTmp = (datTmp - ctr_median)/ctr_median
        if(dat_scaled is None):
            dat_scaled = datScaledTmp
        else:
            dat_scaled = pd.concat([dat_scaled, datScaledTmp])

    #set rglob_cnr for ONLY THE PHOSPHO CELLS
    rglob_cnr = dat_scaled[ab_use_phospho].transpose()
    # Make clean names for the node names
    rglob_cnr.index = [ab_mapping[ab_name] for ab_name in rglob_cnr.index]
    rglob_cnr = rglob_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

    #set rtot for ONLY TOTAL PROTEINS
    #rtot_cnr = dat_scaled[ab_use_tot].transpose()
    #rtot_cnr.index = [ab_mapping[ab_name] for ab_name in rtot_cnr.index]

    nn, an = zip(*[(nn, an) for nn, an in ab_mapping_tot_inv.items()])
    rtot_cnr = dat_scaled[list(an)].transpose()
    rtot_cnr.index = nn
    rtot_cnr = rtot_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

    return(rglob_cnr, rtot_cnr, cell_annot, tx_annot, group_annot)

#due to precision errors substracting rlocs causes UNTRUE differences: therefore store differences only if
#cplexs indicator variable for edge differences is set
def return_true_difference_rloc(rlocA, rlocB, result):
    numDifferences = 0
    differenceRloc = None

    #calculate a first difference rloc
    differenceRloc = (rlocA - rlocB)

    #get all TRUELY different edges
    colsToKeep = []
    for indicator, status in result.vardict.items():
        if(indicator.startswith(("IDev")) and status==1):
            indSplit = indicator.split("_")
            colsToKeep.append([indSplit[1], indSplit[2]])

    #make values same for the non different ones: going over current differences
    #and resetting FALSE ones to zero
    for row in differenceRloc.index:
        for col in differenceRloc.columns:
            if([row, col] not in colsToKeep):
                differenceRloc.loc[row, col] = 0

    numDifferences = np.count_nonzero(differenceRloc)

    return (differenceRloc, numDifferences)

#return the variables for two MRAs, one on the vehicle and one on the EGFR inhibition data
def prepare_data_for_separate_MRAs(dat):
    dat.sample_id = [s.replace("_", ".") for s in dat.sample_id]
    dat = dat.set_index('sample_id')

    resultDict = dict()
    for group in dat.treatment.unique():
        groupRDicts = dict()
        datTmp = dat[dat.treatment == group]

        datTmp = datTmp.pivot(columns='ab_name', values='ab_count_tmm')
        datTmp.reset_index()

        ctr_median = datTmp.median()
        dat_scaled = (datTmp - ctr_median)/ctr_median

        #set rglob_cnr for ONLY THE PHOSPHO CELLS
        rglob_cnr = dat_scaled[ab_use_phospho].transpose()
        # Make clean names for the node names
        rglob_cnr.index = [ab_mapping[ab_name] for ab_name in rglob_cnr.index]
        rglob_cnr = rglob_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

        nn, an = zip(*[(nn, an) for nn, an in ab_mapping_tot_inv.items()])
        rtot_cnr = dat_scaled[list(an)].transpose()
        rtot_cnr.index = nn
        rtot_cnr = rtot_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

        groupRDicts["tot"] =rtot_cnr
        groupRDicts["glob"] = rglob_cnr

        resultDict[group] = groupRDicts

    return(resultDict)

#to ru an MRA looking at the whole data, BE AWARE: we also calcualte deviation from mean based on the average mean over all groups
def prepare_data_for_one_MRA(dat):
    dat.sample_id = [s.replace("_", ".") for s in dat.sample_id]
    dat = dat.set_index('sample_id')

    dat = dat.pivot(columns='ab_name', values='ab_count_tmm')
    dat.reset_index()

    ctr_median = dat.median()
    dat_scaled = (dat - ctr_median)/ctr_median

    #set rglob_cnr for ONLY THE PHOSPHO CELLS
    rglob_cnr = dat_scaled[ab_use_phospho].transpose()
    # Make clean names for the node names
    rglob_cnr.index = [ab_mapping[ab_name] for ab_name in rglob_cnr.index]
    rglob_cnr = rglob_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

    nn, an = zip(*[(nn, an) for nn, an in ab_mapping_tot_inv.items()])
    rtot_cnr = dat_scaled[list(an)].transpose()
    rtot_cnr.index = nn
    rtot_cnr = rtot_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

    groupRDicts = {}
    groupRDicts["tot"] =rtot_cnr
    groupRDicts["glob"] = rglob_cnr

    return(groupRDicts)

def prepare_data_for_EGFR_inhibition(dat):
    #data table
    #sample_id  treatment   ab_count_tmm    ab_name
    dat.sample_id = [s.replace("_", ".") for s in dat.sample_id]
    dat = dat.set_index('sample_id')

    group_annot = dict()
    for group in dat.treatment.unique():
        group_annot[str(group)] = np.unique(list(dat[dat.treatment == group].index))

    print(group_annot)
    #for completeness how cells and treatment would be annotated
    #HOWEVER: we assign NO TREATMENTS and handle the three populations as independant groups
    cell_annot = None
    tx_annot = None
    '''
    cell_annot = dict()
    for tx in dat.treatment.unique():
        cell_annot[tx] = np.unique(list(dat[dat.treatment == tx].index))
    tx_annot = 
    {
        "AG1478" : None, #or "EGFR" if with perturbations
        "vehiclepRB" : None, 
        "vehicle" : None
    }
    '''

    #scale every group by its own median
    #subset data by group
    dat_scaled = None
    for group in dat.treatment.unique():
        datTmp = dat.loc[group_annot[group]]

        datTmp = datTmp.pivot(columns='ab_name', values='ab_count_tmm')
        datTmp.reset_index()

        ctr_median = datTmp.median()
        datScaledTmp = (datTmp - ctr_median)/ctr_median
        #print("just scaled: ")
        #print(datScaledTmp)
        if(dat_scaled is None):
            dat_scaled = datScaledTmp
        else:
            dat_scaled = pd.concat([dat_scaled, datScaledTmp])

    #set rglob_cnr for ONLY THE PHOSPHO CELLS
    rglob_cnr = dat_scaled[ab_use_phospho].transpose()
    # Make clean names for the node names
    rglob_cnr.index = [ab_mapping[ab_name] for ab_name in rglob_cnr.index]
    rglob_cnr = rglob_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

    #set rtot for ONLY TOTAL PROTEINS
    #rtot_cnr = dat_scaled[ab_use_tot].transpose()
    #rtot_cnr.index = [ab_mapping[ab_name] for ab_name in rtot_cnr.index]

    nn, an = zip(*[(nn, an) for nn, an in ab_mapping_tot_inv.items()])
    rtot_cnr = dat_scaled[list(an)].transpose()
    rtot_cnr.index = nn
    rtot_cnr = rtot_cnr.sum(level=0) #sum rows of same index (e.g. two ABs measure same protein)

    return(rglob_cnr, rtot_cnr, cell_annot, tx_annot, group_annot)

def shuffle_group_label_3labels(group_annot):
    keys = list(group_annot.keys())
    keyValueNumber = {}
    valueList= []
    for key in keys:
        keyValueNumber[key] = len(group_annot[key])
        valueList.append(group_annot[key])

    valueList = [item for sublist in valueList for item in sublist]

    group_annot_return = {}
    assert(len(keys) == 3)
    group_annot_return[keys[0]] = random.sample(valueList, keyValueNumber[keys[0]])
    group_annot_return[keys[1]] = random.sample(list(set(valueList) - set(group_annot_return[keys[0]])), keyValueNumber[keys[1]])
    listSoFar = list(set(valueList) - set(group_annot_return[keys[0]]))
    listSoFar = list(set(listSoFar) - set(group_annot_return[keys[1]]))
    group_annot_return[keys[2]] = random.sample(listSoFar, keyValueNumber[keys[2]])

    return(group_annot_return)
