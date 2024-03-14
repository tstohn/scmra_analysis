'''
This script uses the data generated in generateCellNumData:
We look at all the reconstructed data and store for every eta all the edge strength, and weather its a FP or TP edge,
to later visualize the persistence of TP/FP edges with increasing eta.
(When increasing the edge penalty do the FP edges drop out first, and are those edge strength stronger thatn the TP edges???)
'''

INPUT = "../bin/SIMULATION_OF_ETA_CELLNUM_2009_09-25-2022|12:40:45.pickle"
OUTPUT = "../bin/edgeStrengthForIncreasingEta_1000.tsv"
CELLNUM = 1000

import sys, getopt
sys.path.insert(0, '../../helperScripts/')
from helperFunctions import *

def write_edge_strength_for_all_eta_to_file(pickleFile, outputFile):
   #load pickle
   result = pickle.load(open(pickleFile, 'rb'))

   # Create the pandas DataFrame
   headers = result.header
   df = pd.DataFrame(columns=headers)

   #add all the existing parameters
   resultNumber = len(result.parameter)
      
   #extend columns by new variables (TPR, ...)
   #VARIABLES TO FILL
   edgestrengths = []
   istps = []
   edges = []

   #all the param values need to be multiplied, bcs. we have several edges now for the same simulation
   # (e.g. same eta, CellNumber, Noise with X-different edges)
   newDublicatedParams = []
   for i in range(resultNumber):
      #only store the cullNUmber we are interested in: otherwise data overflows
      index = result.header.index("CELLNUMBER")
      if(not (result.parameter[i][index] == CELLNUM)):
         continue

      print(i)

      scMraResult = result.result[i]

      #GET BASIC METRICS TPR, FPR
      rloc =scMraResult.rloc

      colNames = rloc.columns.values
      rowNames = rloc.index.values
      for col in colNames:
         for row in rowNames:
            if(col == row):
               continue
            edges.append(row + "_" + col)
            edgestrengths.append(rloc.loc[row, col])
            if(not rloc_true.loc[row, col] == 0):
               istps.append(1)
            else:
               istps.append(0)
            #also attach the new params
            newDublicatedParams.append(result.parameter[i])
 
   newResultNumber = len(newDublicatedParams)
   for i in range(newResultNumber):
      params = newDublicatedParams[i]
      df_length = len(df)
      df.loc[df_length] = params
   #FINALLY ADD NEW COLUMNS
   df["EDGENAME"] = edges
   df["ISTP"] = istps
   df["EDGESTRENGTH"] = edgestrengths

   df.to_csv(outputFile, sep="\t", index=False)

def main(argv):
   inputfile = ''
   outputfile = ''

   opts = None
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('no input according to: scMRAResult_to_tsv.py -i <inputfile> -o <outputfile> => default to pre-generated files!')
   
   if(opts):
      for opt, arg in opts:
         if opt == '-h':
            print('scMRAResult_to_tsv.py -i <inputfile> -o <outputfile>')
            sys.exit()
         elif opt in ("-i", "--ifile"):
            inputfile = arg
         elif opt in ("-o", "--ofile"):
            outputfile = arg
   else:
      inputfile = INPUT
      outputfile = OUTPUT

   write_edge_strength_for_all_eta_to_file(inputfile,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])