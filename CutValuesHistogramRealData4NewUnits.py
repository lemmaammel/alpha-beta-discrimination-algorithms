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

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/realEventDataNewUnits.cpp+")

rhoCoordinates = [0.056, 0.056, 0.056, 0.056, 0.165, 0.165, 0.165, 0.165, 0.280, 0.280, 0.280, 0.280, 0.392, 0.392, 0.392, 0.504, 0.504]
zCoordinates = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 1, 2]

# set ratio Alpha/Beta
ratio = 9
ratio2 = str(s.significantFigures(ratio,3)).replace(".","-")

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

graphs = [ClassifierYoudenArray, ClassifierGeneralArray, ValueYoudenArray, ValueGeneralArray, AlphaRejectionYoudenArray, BetaAcceptanceYoudenArray, AlphaRejectionGeneralArray, BetaAcceptanceGeneralArray]

titles = ["Youden Classifier Cut Value","General Classifier Cut Value", "Youden Cut Value", "General Cut Value", r"Youden $\alpha$ Rejection", r"Youden $\beta$ Acceptance", r"General $\alpha$ Rejection", r"General $\beta$ Acceptance"]

graphs2 = ["YoudenClassifierCut", "GeneralClassifierCut", "YoudenCutValue", "GeneralCutValue", "YoudenAlphaRejection", "YoudenBetaAcceptance", "GeneralAlphaRejection", "GeneralBetaAcceptance"]

colorbar = ["Classifier Value", "Classifier Value", "Youden Statistic Value", "General Statistic Value", r"$\alpha$ Rejection %", r"$\beta$ Acceptance %", r"$\alpha$ Rejection %", r"$\beta$ Acceptance %"]

titles2 = []

for i in range(0,17):

	values = r.rejectionInfo("/data/snoplus/home/masmiley/snoplusdata/skimmedbiposreprocess/partialFill_data_Po214_reprocess.ntuple.root", "/data/snoplus/home/masmiley/snoplusdata/skimmedbiposreprocess/partialFill_data_Bi214_reprocess.ntuple.root", rhoCoordinates[i], zCoordinates[i]*1000, ratio)
	
	ClassifierYoudenArray.append(values[0])
	ValueYoudenArray.append(values[1])
	ClassifierGeneralArray.append(values[2])
	ValueGeneralArray.append(values[3])
	AlphaRejectionYoudenArray.append(values[4])
	BetaAcceptanceYoudenArray.append(values[5])
	AlphaRejectionGeneralArray.append(values[6])
	BetaAcceptanceGeneralArray.append(values[7])

for graph in graphs:
	graph.extend([-100, -100, -100])

rhoCoordinates.extend(0.392, 0.504, 0.504)
zCoordinates.extend([4, 4, 3])

for i in range(0,8):
	p.hist2d(rhoCoordinates, zCoordinates, bins=(5,4), range=((0.0,0.56),(0.5,4.5)), weights = graphs[i], cmap=p.cm.viridis, cmin=-10)

	for k in range(17):
		array = graphs[i]
		p.text(rhoCoordinates[k], zCoordinates[k], s.significantFigures(array[k], 4), ha="center",va="center",color="w")

	p.xlabel(r"$\rho^2/36$ coordinate")
	p.ylabel(r"$z$ coordinate (m)")
	p.title(titles[i])
	p.xticks([0.112, 0.224, 0.336, 0.448, 0.560])
	p.yticks([1.0, 2.0, 3.0, 4.0])
	p.colorbar(label=colorbar[i])
	if i == 0 or i == 1: p.clim(-0.008, 0.004)
	if i == 2: p.clim(0.65, 0.75)
	if i == 3: p.clim(20,50)
	if i==4 or i==5: p.clim(0.7, 1.0)
	if i==6 or i==7: p.clim(0.7, 1.0)
	p.show()
	p.savefig("SummaryPartialFillRealDataSkimmedBiPosReprocessed_{}.pdf".format(graphs2[i]))
	p.clf()
