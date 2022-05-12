# still working on commenting and refactoring!!!

#Useful generic python imports
from __future__ import print_function
from string import Template
import os, sys, time, math
import numpy as np
from array import array
import matplotlib as m
m.use('Agg')

#Imports particularly important for our purposes
import ROOT as r
import rat
import csv
import math
import matplotlib.pyplot as p
import significantFigures as s
import findRatio as f
import os
import csv

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/realDataValues.cpp+")

loopvars = ["0_0_2","0_0_3","0_0_4", "1_0_2","1_0_3","1_0_4", "2_0_2","2_0_3","2_0_4","3_0_2","3_0_3","4_0_2"]

# set ratio Alpha/Beta
#ratio = f.findRatio(0.78,0.78)
#ratio2 = str(s.significantFigures(ratio,3)).replace(".","-")

YoudenCutValue =[]
GeneralCutValue=[]
AlphaRejection=[]
BetaAcceptance=[]

coordinates1 = []
coordinates2 = []

graphs = [YoudenCutValue, GeneralCutValue, AlphaRejection, BetaAcceptance]

titles = ["Youden Statistic Value","General Statistic Value", r"$\alpha$ Rejection", r"$\beta$ Acceptance"]

graphs2 = ["YoudenCutValue", "GeneralCutValue", "AlphaRejection", "BetaAcceptance"]

colorbar = ["Youden Statistic Value", "General Statistic Value", r"$\alpha$ Rejection", r"$\beta$ Acceptance"]

titles2 = []

cutoffs = []

for title in titles:
        newTitle = title + r" (#$\alpha$/#$\beta$ = " + str(1) + ")"
        titles2.append(newTitle)

for loopvar in loopvars:

        cutoff = 0

        with open("summaryClassifierCutoffsPartial.csv") as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=",")
                for row in csv_reader:
                        if row[0] == loopvar:
                                cutoff = row[1]
                                cutoffs.append(cutoff)

        values = r.cutPerformanceValues("/data/snoplus/home/masmiley/skimmedbiposreprocess/partialFill_data_Po214_reprocess.ntuple.root", "/data/snoplus/home/masmiley/skimmedbiposreprocess/partialFill_data_Bi214_reprocess.ntuple.root", float(loopvar[0])*1000, float(loopvar[4])*1000, 9, float(cutoff))
        
        YoudenCutValue.append(values[0])
        GeneralCutValue.append(values[1])
        AlphaRejection.append(values[2])
        BetaAcceptance.append(values[3])
        
        coordinates1.append(float(loopvar[0]))
        coordinates2.append(float(loopvar[4]))

for graph in graphs:
        graph.extend([-100, -100, -100])

coordinates1.extend([3, 4, 4])
coordinates2.extend([4, 4, 3])

for i in range(0,4):
        p.hist2d(coordinates1, coordinates2, bins=(5,3), range=((-0.5,4.5),(1.5,4.5)), weights = graphs[i], cmap=p.cm.viridis, cmin=-10)

        for k in range(12):
                array = graphs[i]
                p.text(coordinates1[k], coordinates2[k], s.significantFigures(array[k], 4), ha="center",va="center",color="w")

        p.xlabel(r"$\rho$ coordinate (m)")
        p.ylabel(r"$z$ coordinate (m)")
        p.title(titles2[i])
        p.colorbar(label=colorbar[i])
        if i == 0 or i == 1: p.clim(0.004, 0.016)
        if i == 2: p.clim(0.4, 0.5)
        if i == 3: p.clim(6,38)
        if i==4 or i==5: p.clim(0.2, 1.1)
        if i==6 or i==7: p.clim(0.2, 1.1)
        p.show()
        p.savefig("realDataPDFs/ComparisonRealData{}.pdf".format(graphs2[i]))
        p.clf()
