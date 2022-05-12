# compiles information about the classifier cutoff values across positions

#Useful generic python imports
from __future__ import print_function
from string import Template
import os, sys, time, math
import numpy as np
from array import array

#Imports particularly important for our purposes
import ROOT as r
import rat
import csv
import math

loopvars = ["0_0_2","0_0_3", "0_0_4", "1_0_2","1_0_3","1_0_4","2_0_2","2_0_3","2_0_4","3_0_2","3_0_3","4_0_2"]

alphaFilename = raw_input("Please enter the alpha filename:")
betaFilename = raw_input("Please enter the beta filename:")


r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/rat/example/root/NhitHistogram.cpp+")

fields = ["location", "classifierCutoffGeneral"]
rows = []
filename = "summaryClassifierCutoffsPartial.csv"

for loopvar in loopvars:

        name1 =  "alphaFilename*".format(loopvar)
        name2 = "betaFilename*".format(loopvar)
        
        values = r.AlphaRejectionInfo("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", 1)
        
        currentRow = [loopvar, values[2]]
        rows.append(currentRow)

with open(filename, "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)

