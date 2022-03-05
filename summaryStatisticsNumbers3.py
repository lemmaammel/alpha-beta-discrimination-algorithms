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

loopvars = ["/stuff_for_emma/april2021_","/stuff_for_emma/sept2021_","/snoplusdata/jan2022skimmedbipo/jan2022_"]

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/realEventData2NewUnits.cpp+")

fields = ["name ","rejection_general ","acceptance_general ","rejection_youden ","acceptance_youden ","cutoff_general ", "cutoff_youden ", "mean_rho_alpha ", "mean_z_alpha ", "mean_rho_beta", "mean_z_beta"]
rows = []
filename = "summaryStatistics3.csv"

for loopvar in loopvars:
	
	name1 =  "/data/snoplus/home/masmiley{}po214".format(loopvar)
	name2 = "/data/snoplus/home/masmiley{}bi214".format(loopvar)
	print(name1)
	values = r.rejectionInfo("{}.root".format(name1), "{}.root".format(name2), 0.0, 0.0,1.0)
	alphaValues = r.averageValues("{}.root".format(name1))
	betaValues = r.averageValues("{}.root".format(name2))		
	currentRow = [loopvar, values[6], values[7], values[4],values[5], values[2], values[0], alphaValues[3], alphaValues[2], betaValues[3], betaValues[2]]
	rows.append(currentRow)

with open(filename, "w") as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(fields)
	csvwriter.writerows(rows)

