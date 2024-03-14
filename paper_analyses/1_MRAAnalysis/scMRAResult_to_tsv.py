import sys, getopt
sys.path.insert(0, '../../helperScripts/')
from helperFunctions import *

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

   write_metrics_to_csv_from_pickle_MRA(inputfile,outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])