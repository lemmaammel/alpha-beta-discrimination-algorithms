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

loopvars = ["0_0_2","0_0_3", "0_0_4","1_0_2","1_0_3","1_0_4","2_0_2","2_0_3","2_0_4"]

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

graphs = [ClassifierYoudenArray, ClassifierGeneralArray, ValueYoudenArray, ValueGeneralArray, AlphaRejectionYoudenArray, BetaAcceptanceYoudenArray, AlphaRejectionGeneralArray, BetaAcceptanceGeneralArray]
titles = ["Youden Classifier Cut Value","General Classifier Cut Value", "Youden Cut Value", "General Cut Value", r"Youden $\alpha$ Rejection", r"Youden $\beta$ Acceptance", r"General $\alpha$ Rejection", r"General $\beta$ Acceptance"]
graphs2 = ["YoudenClassifierCut", "GeneralClassifierCut", "YoudenCutValue", "GeneralCutValue", "YoudenAlphaRejection", "YoudenBetaAcceptance", "GeneralAlphaRejection", "GeneralBetaAcceptance"]
colorbar = ["Classifier Value", "Classifier Value", "Youden Statistic Value", "General Statistic Value", r"$\alpha$ Rejection", r"$\beta$ Acceptance", r"$\alpha$ Rejection", r"$\beta$ Acceptance"]
titles2 = []

fields = ["x ","y ","z ","Youden Classifier Cut ","General Classifier Cut ","Youden Cut Value ","General Cut Value ", "Youden Alpha Rejection ", "Youden Beta Acceptance ", "General Alpha Rejection ", "General Beta Acceptance "]
rows = [1,0,3]
filename = "summaryPartialFill.csv"


for title in titles:
	newTitle = title
	titles2.append(newTitle)


for loopvar in loopvars:

	name1 =  "/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_e-*".format(loopvar)
	name2 = "/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_alpha*".format(loopvar)
   
	values = r.AlphaRejectionInfo("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", 9)
	
	ClassifierYoudenArray.append(values[0])
	ValueYoudenArray.append(values[1])
	ClassifierGeneralArray.append(values[2])
	ValueGeneralArray.append(values[3])
	AlphaRejectionYoudenArray.append(values[4])
	BetaAcceptanceYoudenArray.append(values[5])
	AlphaRejectionGeneralArray.append(values[6])
	BetaAcceptanceGeneralArray.append(values[7])

for i in range(0,8):
	number = 0.0
	for k in range(9):
		array = graphs[i]
		number+=array[k]

	number = number/9
	rows.append(number)

	p.hist2d([1], [3], bins=(1,1), range=((-.5,2.5),(1.5,4.5)), weights = [number], cmap=p.cm.viridis, cmin=-10)
		
	number = number/9
	p.text(1, 3, s.significantFigures(number, 4), ha="center",va="center",color="w")

	p.xlabel(r"$\rho$ coordinate (m)")
	p.ylabel(r"$z$ coordinate (m)")
	p.title(titles2[i])
	p.colorbar(label=colorbar[i])
	if i == 0 or i == 1: p.clim(0.009, 0.012)
	if i == 2: p.clim(0.85, 1.01)
	if i == 3: p.clim(85,100)
	if i==4 or i==6: p.clim(0.1, 0.4)
	if i==5 or i==7: p.clim(0.95, 1.0)
	p.show()
	p.savefig("partialFillPDFs/SummaryPartialFill{}_Total.pdf".format(graphs2[i]))
	p.clf()

with open(filename, "w") as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerow(fields)
	csvwriter.writerows([rows])

