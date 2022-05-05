# Generates histograms displaying the distribution of alpha rejection and beta acceptance values (among others) around the detector given montecarlo simulated data
# asks user for filename in format "{filename}{x_coordinate}_{y_coordinate}_{z_coordinate}.root"

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

loopvars = ["0_0_2","0_0_3", "0_0_4","1_0_2","1_0_3","1_0_4","2_0_2","2_0_3","2_0_4","3_0_2","3_0_3","4_0_2"]

filename = raw_input("Please enter the filename in format: {filename}{x_coordinate}_{y_coordinate}_{z_coordinate}.root")

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/rat/example/root/NhitHistogram.cpp+")

# set ratio Alpha/Beta
#ratio = f.findRatio(0.78,0.78)
#ratio2 = str(s.significantFigures(ratio,3)).replace(".","-")

ClassifierYoudenArray = []
ValueYoudenArray = []
ClassifierGeneralArray = []
ValueGeneralArray = []

AlphaRejectionYoudenArray = []
BetaAcceptanceYoudenArray = []
AlphaRejectionGeneralArray = []
BetaAcceptanceGeneralArray = []

coordinates1 = []
coordinates2 = []

graphs = [ClassifierYoudenArray, ClassifierGeneralArray, ValueYoudenArray, ValueGeneralArray, AlphaRejectionYoudenArray, BetaAcceptanceYoudenArray, AlphaRejectionGeneralArray, BetaAcceptanceGeneralArray, ClassifierAlphaArray, NhitAlphaArray, ClassifierBetaArray, NhitBetaArray]
titles = ["Youden Classifier Cut Value","General Classifier Cut Value", "Youden Cut Value", "General Cut Value", r"Youden $\alpha$ Rejection", r"Youden $\beta$ Acceptance", r"General $\alpha$ Rejection", r"General $\beta$ Acceptance", r"Classification $\alpha$ Summary",r"$N_{\mathrm{hit}}$ $\alpha$ Summary", r"Classification $\beta$ Summary", r"$N_{\mathrm{hit}}$ $\beta$ Summary"]
graphs2 = ["YoudenClassifierCut", "GeneralClassifierCut", "YoudenCutValue", "GeneralCutValue", "YoudenAlphaRejection", "YoudenBetaAcceptance", "GeneralAlphaRejection", "GeneralBetaAcceptance", "ClassifierAlphaArray", "NhitAlphaArray", "ClassifierBetaArray", "NhitBetaArray"]
colorbar = ["Classifier Value", "Classifier Value", "Youden Statistic Value", "General Statistic Value", r"$\alpha$ Rejection", r"$\beta$ Acceptance", r"$\alpha$ Rejection", r"$\beta$ Acceptance", "Classifier Value", r"$N_{\mathrm{hit}}$ Value", "Classifier Value", r"$N_{\mathrm{hit}}$ Value"]


for loopvar in loopvars:

	name1 =  "filename*".format(loopvar)
	name2 = "filename*".format(loopvar)
   
	values = r.AlphaRejectionInfo("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", 9)
	
	ClassifierYoudenArray.append(values[0])
	ValueYoudenArray.append(values[1])
	ClassifierGeneralArray.append(values[2])
	ValueGeneralArray.append(values[3])
	AlphaRejectionYoudenArray.append(values[4])
	BetaAcceptanceYoudenArray.append(values[5])
	AlphaRejectionGeneralArray.append(values[6])
	BetaAcceptanceGeneralArray.append(values[7])
	
	htot = r.NhitHistogram("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", "alpha")
	htot2 = r.NhitHistogram("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", "beta")
	
	ClassifierAlphaArray.append(htot.GetMean(2))
	NhitAlphaArray.append(htot.GetMean(1))
	ClassifierBetaArray.append(htot2.GetMean(2))
	NhitBetaArray.append(htot2.GetMean(1))
	
	coordinates1.append(float(loopvar[0]))
	coordinates2.append(float(loopvar[4]))

for graph in graphs:
	graph.extend([-100, -100, -100])

coordinates1.extend([3, 4, 4])
coordinates2.extend([4, 4, 4])

for i in range(0,12):
	p.hist2d(coordinates1, coordinates2, bins=(5,3), range=((-.5,4.5),(1.5,4.5)), weights = graphs[i], cmap=p.cm.viridis, cmin=-10)

	for k in range(12):
		array = graphs[i]
		p.text(coordinates1[k], coordinates2[k], s.significantFigures(array[k], 4), ha="center",va="center",color="w")

	p.xlabel(r"$\rho$ coordinate (m)")
	p.ylabel(r"$z$ coordinate (m)")
	p.title(titles[i])
	p.colorbar(label=colorbar[i])
	if i == 0 or i == 1: p.clim(0.009, 0.012)
	if i == 2: p.clim(0.85, 1.01)
	if i == 3: p.clim(85,100)
	if i==4 or i==6: p.clim(0.1, 0.4)
	if i==5 or i==7: p.clim(0.95, 1.0)
	if(i==8 or i==9): p.clim(-0.01, 0.01)
	else: p.clim(240,250)
	p.show()
	p.savefig("partialFillPDFs/SummaryPartialFill{}.pdf".format(graphs2[i]))
	p.clf()
