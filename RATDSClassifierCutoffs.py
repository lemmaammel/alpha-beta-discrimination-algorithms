# compiles information about the classifier cutoff values across positions

#Useful generic python imports
from __future__ import print_function
from string import Template
import os, sys, time, math, csv, argparse, glob
import numpy as np
from array import array

#Imports particularly important for our purposes
import ROOT as r
import rat

loopvars = ["0_0_2","0_0_3", "0_0_4", "1_0_2","1_0_3","1_0_4","2_0_2","2_0_3","2_0_4","3_0_2","3_0_3","4_0_2"]

#alphaFilename = raw_input("Please enter the alpha filename:")
#betaFilename = raw_input("Please enter the beta filename:")

parser = argparse.ArgumentParser()
parser.add_argument('--alphafile', '-a', type = str, default = '', help = 'Alpha file to use')
parser.add_argument('--betafile', '-b', type = str, default = '', help = 'Beta file to use')
parser.add_argument('--outfile', '-o', type = str, default = 'output.csv', help = 'Name of file to output')

args = parser.parse_args()


r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/rat/example/root/NhitHistogram.cpp+")

fields = ["location", "classifierCutoffGeneral"]
rows = []
#filename = "summaryClassifierCutoffsPartial.csv"
filename = args.outfile

#We can phase out the loopvars array as predefined if we're smart about naming patterns
#We can use the glob package to our advantage, for example
#files = glob.glob('myfile_*.root') would return a list of all files in the current directory
#that begin with the string 'myfile_' and end with '.root'
#so it could be 'myfile_1.root', 'myfile_alpha.root', etc.
#but we'll have to come up with some conventions, which could be flexible and specified
#as another argument on the command line

for loopvar in loopvars:

        #name1 =  "alphaFilename*".format(loopvar)
        #name2 = "betaFilename*".format(loopvar)
        name1 = args.alphafile 
        name2 = args.betafile
        
        values = r.AlphaRejectionInfo("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", 1)
        
        currentRow = [loopvar, values[2]]
        rows.append(currentRow)

with open(filename, "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)
