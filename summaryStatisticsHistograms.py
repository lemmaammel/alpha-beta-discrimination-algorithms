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

copies = 3
loopvars = ["0_0_2","0_0_3", "0_0_4", "1_0_2","1_0_3","1_0_4","2_0_2","2_0_3","2_0_4","3_0_2","3_0_3","4_0_2"]

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/rat/example/root/NhitHistogram.cpp+")


ClassifierAlphaArray = []
NhitAlphaArray = []
ClassifierBetaArray = []
NhitBetaArray = []
coordinates1 = []
coordinates2 = []

graphs = [ClassifierAlphaArray, NhitAlphaArray, ClassifierBetaArray, NhitBetaArray]
titles = [r"Classification $\alpha$ Summary",r"$N_{\mathrm{hit}}$ $\alpha$ Summary", r"Classification $\beta$ Summary", r"$N_{\mathrm{hit}}$ $\beta$ Summary"]
graphs2 = ["ClassifierAlphaArray", "NhitAlphaArray", "ClassifierBetaArray", "NhitBetaArray"]
colorbar = ["Classifier Value", r"$N_{\mathrm{hit}}$ Value", "Classifier Value", "$N_{\mathrm{hit}}$ Value"]

for loopvar in loopvars:

	name1 =  "/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_e-*".format(loopvar)
	name2 = "/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_alpha*".format(loopvar)
   
	htot = r.NhitHistogramAlpha("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood")
	htot2 = r.NhitHistogramBeta("{}.root".format(name1), "{}.root".format(name2), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood")
	
	ClassifierAlphaArray.append(htot.GetMean(2))
	NhitAlphaArray.append(htot.GetMean(1))
	ClassifierBetaArray.append(htot2.GetMean(2))
	NhitBetaArray.append(htot2.GetMean(1))
	coordinates1.append(float(loopvar[0]))
	coordinates2.append(float(loopvar[4]))

for graph in graphs:
	graph.append(-100)
	graph.append(-100)
	graph.append(-100)
	print(len(graph))

coordinates1.append(3)
coordinates1.append(4)
coordinates1.append(4)
coordinates2.append(4)
coordinates2.append(4)
coordinates2.append(3)

for i in range(0,4):
	p.hist2d(coordinates1, coordinates2, bins=(5,3), range=((-0.5,4.5),(1.5,4.5)), weights = graphs[i], cmap=p.cm.viridis, cmin=-10)
	

	for k in range(12):
		array = graphs[i]
		p.text(coordinates1[k], coordinates2[k], s.significantFigures(array[k], 4), ha="center",va="center",color="w")

	p.xlabel(r"$\rho$ coordinate (m)")
	p.ylabel(r"$z$ coordinate (m)")
	p.title(titles[i])
	p.colorbar(label=colorbar[i])
	if(i==0 or i==2): p.clim(-0.01, 0.01)
	else: p.clim(240,250)
	p.show()
	p.savefig("Summary{}.pdf".format(graphs2[i]))
	p.clf()
