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

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/rat/example/root/NhitHistogram.cpp+")

fields = ["location", "classifierCutoffGeneral"]
rows = []
filename = "summaryClassifierCutoffsPartial.csv"

for loopvar in loopvars:

	name1 =  "/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_e-*".format(loopvar)
	name2 = "/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_alpha*".format(loopvar)
	
	values = r.AlphaRejectionInfo("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", 1)
	
	currentRow = [loopvar, values[2]]
	rows.append(currentRow)

with open(filename, "w") as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(fields)
	csvwriter.writerows(rows)

