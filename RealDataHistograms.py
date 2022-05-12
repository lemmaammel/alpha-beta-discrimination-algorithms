# Generates histograms displaying the distribution of alpha rejection and beta acceptance values (among others) around the detector given real simulated data (NEW UNITS)

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
r.gROOT.LoadMacro("/data/snoplus/home/ammel/RealDataValues.cpp+")

alphaFilename = raw_input("Please enter the alpha filename:")
betaFilename = raw_input("Please enter the beta filename:")
inputType = raw_input("Please enter if you would like to type your own list for the rho and z coordinates or if you would like to use a square template ('list' OR 'square'):")

rhoCoordinates = []
zCoordinates = []
distance = 1

if inputType == "square"
        squareLength = int(raw_input("Please enter the side length of your square:"))
        distance = int(raw_input("Please enter the distance between the coordinates:))
        for i in range(-math.floor(squareLength/distance, math.floor(squareLength/distance)
                for j in range(-math.floor(squareLength/distance, math.floor(squareLength/distance)
                        rhoCoordinates.extend(i)
                        zCoordinates.extend(j)
                       
if inputType == "list"
        rhoCoordinates = list(map(int, raw_input("Please enter the rho coordinates in the form '3 2 3 4 8':").split()))
        zCoordinates = list(map(int, raw_input("Please enter the z coordinates in the form '4 6 3 8 8':").split()))
        distance = int(raw_input("Please enter the distance between the coordinates:))
                                 
xTicks = []
yTicks = []

for i in range(0, math.floor(min(rhoCoordinates)-max(rhoCoordinates)/distance))
        xTicks.append(min(rhoCoordinates) + (distance*i))
                                 
for i in range(0, math.floor(min(zCoordinates)-max(zCoordinates)/distance))
        yTicks.append(min(zCoordinates) + (distance*i))
                

# rhoCoordinates = [0.056, 0.056, 0.056, 0.056, 0.165, 0.165, 0.165, 0.165, 0.280, 0.280, 0.280, 0.280, 0.392, 0.392, 0.392, 0.504, 0.504]
# zCoordinates = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 1, 2]

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

        values = r.rejectionInfo(alphaFilename, betaFilename, rhoCoordinates[i], zCoordinates[i], ratio, distance)
        
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
        p.xticks(xTicks)
        p.yticks(yTicks)
        p.colorbar(label=colorbar[i])
        p.show()
        p.savefig("SummaryPartialFillRealDataSkimmedBiPosReprocessed_{}.pdf".format(graphs2[i]))
        p.clf()
