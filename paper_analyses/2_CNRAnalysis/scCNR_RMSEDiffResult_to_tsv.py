'''
Calculates the RMSE of the true and reconstructed differences between wt and mut rlocs
used in the paper mainly for the MRA & CNR comparison
'''

import sys, getopt, math
sys.path.insert(0, '../../helperScripts/')
from helperFunctions import *

OUTFILE = "CnrMraComparison.tsv"

def get_edge_number(rloc):
   numEdges = 0
   for column_name in rloc:
      column = rloc[column_name]
      numEdges += (column != 0).sum()
      numEdges -=1 #for the diagonal values
   return(numEdges)

def calculate_rmse_of_network_differences(rlocList, populationStringList, diffDict):

   rloc_wt = None
   rloc_mut = None

   for popIdx in range(len(populationStringList)):
      pop = populationStringList[popIdx]
      if(pop == "wt"):
         rloc_wt = rlocList[popIdx]
      else:
         mutString = pop
         rloc_mut = rlocList[popIdx]

   assert(rloc_wt is not None)
   assert(rloc_mut is not None)
   assert(mutString is not "wt")

   #true difference rloc matrix
   trueDiff = diffDict[mutString]

   #reconstructed difference rloc matrix
   rloc_wt = reorder_rloc(rloc_wt)
   rloc_mut = reorder_rloc(rloc_mut)
   reconstructedDiff = ((rloc_wt - rloc_mut)**2)**0.5

   #from these final two rlocs calculate RMSE
   trueDiff = reorder_rloc(trueDiff)
   reconstructedDiff = reorder_rloc(reconstructedDiff)

   rmse = calc_rmse(trueDiff, reconstructedDiff)
   return(rmse)

def get_population_list_from_string(string):
   popList = []
   x = string.split(",")
   for pop in x:
      if pop == "": 
         continue
      else:
         popList.append(pop)
   return(popList)

def calculate_diff_dict():
   dict = {}
   dict['braf'] = ((rloc_true - rloc_braf)**2)**0.5
   dict['ras'] = ((rloc_true - rloc_ras)**2)**0.5

   return(dict)

#not all against all: for MRA we have used same data as for CNR
#and therefore can directly compare the differences of the specific networks
#networks need to be mapped by REPEAT: so we have a ctr and braf population MRA
#with the same REPEAT number
def write_rloc_differences_to_file(pickleFile, outputFile,):
   result = pickle.load(open(pickleFile, 'rb'))

   print(result.header)
   print(result.parameter)

   diffDict = calculate_diff_dict()

   # Create the pandas DataFrame
   df = pd.DataFrame(columns=["POPULATIONS", "DIFFERENCE_RMSE"] + result.header)
   
   # for every combination of
   # NOISE CELLNUM REPEAT of CNR
   headers = result.header
   resultFrame = pd.DataFrame(columns=headers)
   #add all the existing parameters
   resultNumber = len(result.parameter)
   for i in range(resultNumber):
      params = result.parameter[i]
      df_length = len(resultFrame)
      resultFrame.loc[df_length] = params
   cnrData = resultFrame[resultFrame["MODELTYPE"] == "CNR"]

   #get noise cellnum and repeats for CNR
   noiseLevels = cnrData["NOISE"].unique()
   cellnumLevels = cnrData["CELLNUM"].unique()
   repeatLevels = cnrData["REPEAT"].unique()
   thetaLevels = cnrData["THETA"].unique()

   print("Variables to analyze: ")
   print(noiseLevels)
   print(cellnumLevels)
   print(repeatLevels)

   #get populations that r compared for CNR (taking only the ones with 2 populations for now)
   populationLevelsTmp = cnrData["POPULATIONS"].unique()
   populationLevels = []
   for level in populationLevelsTmp:
      x = level.split(",")
      populationNumber = 0
      for pop in x:
         if pop == "": continue
         populationNumber += 1

      if(populationNumber == 2):
         populationLevels.append(level)

   #calculate Difference RMSE for all combinations of CNRs
   #for now only look at combinations with 1 wt and 1 mutated (removing 3 populations & e.g. ras,braf ones)
   keyList = ["MODELTYPE", "NOISE", "CELLNUM", "POPULATIONS","REPEAT","THETA","EDGES"]
   for noise in noiseLevels:
      for cellnum in cellnumLevels:
         for pop in populationLevels:
            for repeat in repeatLevels:
               for theta in thetaLevels:
                  print(populationLevels)
                  print(repeatLevels)
                  valueList = ["CNR", noise, cellnum, pop, repeat, theta, 13]
                  print(valueList)
                  popList = get_population_list_from_string(pop) # list of populations like ["wt", "braf"]
                  if(len(popList) != 2): 
                     print("NOT EXACTLY 2 POPULATIONS FOR MRA -> SKIP MRA RLOC DIFF CALCUALTION\n")
                     continue
                  if(not "wt" in popList and not "WT" in popList): 
                     print("NO wildtype <wt> in population list (needed for rloc Diff) -> SKIP MRA RLOC DIFF CALCUALTION \n")
                     continue

                  resultObject, params = result.getResult(keyList, valueList) #the result from the simulation
                  if(resultObject == None): continue

                  vars_lst = [var for var in resultObject.vardict if var.startswith('IDev')]
                  print("NUMBER OF DIFFERENT EDGES: ")
                  vars_vals = [resultObject.vardict[v] for v in vars_lst]
                  lengthIdcs = (np.count_nonzero(vars_vals))
                  print(lengthIdcs)

                  assert(len(popList) == 2)
                  rloc_list = []
                  #we assign population by wt, braf, ras. Those MUST BE LOWER CASE NAMES.
                  #just in case convert to lower case beforehand
                  lowerPopList = []
                  for popTmp in popList:
                     lowerPopList.append(popTmp.lower())
                  popList = lowerPopList
                  for popStringIdx in range(len(popList)):
                     popString = popList[popStringIdx]
                     rloc_tmp = resultObject.rloc[popString]
                     rloc_tmp = reorder_rloc(rloc_tmp)
                     rloc_list.append(rloc_tmp)

                  #get a list of network RMSE differences for two rlocs
                  diffRMSE = calculate_rmse_of_network_differences(rloc_list, popList, diffDict)

                  #write new row in result data frame
                  # dataframe struture is: columns=["NOISE", "CELLNUM", "REPEAT", "POPULATIONS", "DIFFERENCE_RMSE"]
                  list1 = [str(popList[0]) + "_" + str(popList[1]), diffRMSE]
                  list2 = params
                  rowList = [*list1, *list2] 
                  print(rowList)
                  print(df.columns)
                  df.loc[len(df)] = rowList

   #---------------------------
   #MRA rloc Differences
   #---------------------------
   #select all the MRAs with same values (same populations, same noise)
   #cellnumber: in CNR cellnumber x means that x cells are picked for each population, so we take THE SAME cellnum for each independant MRA

   #get MRA data
   headers = result.header
   resultFrame = pd.DataFrame(columns=headers)
   #add all the existing parameters
   resultNumber = len(result.parameter)
   for i in range(resultNumber):
      params = result.parameter[i]
      df_length = len(resultFrame)
      resultFrame.loc[df_length] = params
   mraData = resultFrame[resultFrame["MODELTYPE"] == "MRA"]

   if(not mraData.empty):

      #get all possible populations to look into for MRA
      keyList = ["MODELTYPE", "NOISE", "CELLNUM", "POPULATIONS","REPEAT"]
      #iterate over all repeats: 
      #   get the TWO MRA simulations for the repeat and compare the rloc differences
      #   to ground truth

      #get the string of the two populations to test
      resultStrings = []
      print("For MRA we have following population sets to compare: ")
      print(populationLevels)
      '''for pop in populationLevels:
         if("," in pop):
            resultStrings.append(pop)
      for i in range(len(resultStrings)):
         
         popListString = resultStrings[i]'''
      for popList in [ ["WT:wt_braf", "BRAF:wt_braf"], ["WT:wt_ras","RAS:wt_ras"]]:
         for repeatId in repeatLevels:
            for noise in noiseLevels:
               for cellnum in cellnumLevels:
                  #popList = get_population_list_from_string(popListString)
                  print(popList)
                  print(repeatId)
                  print(noise)
                  #if(len(popList) != 2): 
                  #   print("NOT EXACTLY 2 POPULATIONS FOR MRA -> SKIP MRA RLOC DIFF CALCUALTION\n")
                  #   continue
                  #if(not "wt" in popList and not "WT" in popList): 
                  #   print("NO wildtype <wt> in population list (needed for rloc Diff) -> SKIP MRA RLOC DIFF CALCUALTION \n")
                  #   continue

                  #calculate the RMSE of rloc differences for ALL COMBINATIONS of the two populations (for now wt, mut)
                  #by that we end up with many more samples for MRA in the end than for CNR
                  #now we have to compare all to all of the two populations
                  #pop1 = popList[0]
                  #pop2 = popList[1]
                  

                  #get result returns an array of result objects
                  #however the last parameter enforces a check of getResult that ONLY ONE result is returned
                  #this still has to then be subsetted for further rloc diff calculation

                  #GET FIRST RESULT
                  valueList = ["MRA", noise, cellnum, popList[0], repeatId]
                  result1, params1 = result.getResult(keyList, valueList, False, True) #the result from the simulation

                  assert(len(result1) == 1)           
                  result1 = result1[0]
                  params1 = params1[0]
                  popListSingle = popList[0].split(":")
                  pop1 = popListSingle[0]
                  assert(pop1 == "WT")

                  #GET SECOND RESULT
                  valueList = ["MRA", noise, cellnum, popList[1], repeatId]
                  result2, params2 = result.getResult(keyList, valueList, False, True) #the result from the simulation
                  assert(len(result2) == 1)           
                  result2 = result2[0]
                  params2 = params2[0]
                  popListSingle = popList[1].split(":")
                  pop2 = popListSingle[0]
                  assert(pop2 == "BRAF" or pop2 == "RAS")

                  print(pop1)
                  print(pop2)

                  populationList = []

                  populationList.append(pop1.lower())
                  populationList.append(pop2.lower())
                  #calculate RMSE for ALL COMBINATIONS

                  rloc_list = []
                  rloc_tmp1 = reorder_rloc(result1.rloc)
                  rloc_tmp2 = reorder_rloc(result2.rloc)
                  rloc_list.append(rloc_tmp1)
                  rloc_list.append(rloc_tmp2)
                  #get a list of network RMSE differences for two rlocs
                  diffRMSE = calculate_rmse_of_network_differences(rloc_list, populationList, diffDict)

                  #write new row in result data frame
                  # dataframe struture is: columns=["NOISE", "CELLNUM", "REPEAT", "POPULATIONS", "DIFFERENCE_RMSE"]
                  list1 = [str(pop1) + "_" + str(pop2), diffRMSE]
                  list2 = params1
                  rowList = [*list1, *list2] 
                  df.loc[len(df)] =  rowList

   #finally print the dataframe
   df.to_csv(outputFile, index=False, sep="\t")

def main(argv):
   inputfile = ''
   outputfile = ''

   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])

   except getopt.GetoptError:
      print('scMRAResult_to_tsv.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('scMRAResult_to_tsv.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg

   write_rloc_differences_to_file(inputfile,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])