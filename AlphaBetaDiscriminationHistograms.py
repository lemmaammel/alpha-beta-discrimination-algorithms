# Generates histograms displaying the distribution of alpha rejection and beta acceptance values (among others) around the detector given real simulated data (NEW UNITS)

# Python imports
from __future__ import print_function
from string import Template
import os, sys, time, math
import numpy as np
from array import array
import matplotlib as m
import argparse
m.use('Agg')
import ROOT as r
import rat
import csv
import math
import matplotlib.pyplot as p
import significantFigures as s
import os

# Function that calculates relevant statistics from given alpha and beta event histograms
def rejectionInfo(alpha_hist, beta_hist, ratio):
                                           
    alpha_histogram = alpha_hist
    beta_histogram = beta_hist  

    allAlphas = alpha_histogram.Integral()
    allBetas = beta_histogram.Integral()

    if allAlphas==0:
        allAlphas = 1e-15
   
    if allBetas == 0:
        allBetas = 1e-15                                     
                                           
    meanNhit = (alpha_histogram.GetMean(1)+beta_histogram.GetMean(1))/(allAlphas+allBetas)

    # cut selection histograms
    youden_histogram = r.TH1D("Youden's J Statistic", "Youden's J Statistic", alpha_histogram.GetNbinsY(),  alpha_histogram.GetMinimum(), alpha_histogram.GetMaximum())
    general_histogram = r.TH1D("General Cut Statistic", "General Cut Statistic", alpha_histogram.GetNbinsY(), alpha_histogram.GetMinimum(), alpha_histogram.GetMaximum())
                                           
    alpha_rejection = 0
    alpha_acceptance = 0
    beta_rejection = 0
    beta_acceptance = 0
    x = general_histogram.GetXaxis().GetXmin()

    youden_statistic = 0
    general_statistic = 0

    for k in range(0, youden_histogram.GetNbinsX()):
        alpha_rejection = ratio*alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), k, youden_histogram.GetNbinsY())
        alpha_acceptance = ratio*allAlphas - alpha_rejection
        beta_acceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, k)
        beta_rejection = allBetas - beta_acceptance

        if not ((beta_acceptance == 0 and alpha_rejection == 0) or (beta_acceptance+beta_rejection == 0 or alpha_acceptance+alpha_rejection == 0)):
            youden_statistic = beta_acceptance/(beta_acceptance+beta_rejection) + alpha_rejection/(alpha_acceptance+alpha_rejection)
            general_statistic = beta_acceptance/math.sqrt(beta_acceptance+alpha_rejection)

        else:
            youden_statistic = 0
            general__statistic = 0

        #fill for cut selection stats
        youden_histogram.Fill(x, youden_statistic)
        general_histogram.Fill(x, general_statistic)

        x += abs(youden_histogram.GetXaxis().GetXmin() - youden_histogram.GetXaxis().GetXmax())/youden_histogram.GetNbinsX()
    

    youdenClassifierBin = youden_histogram.GetMaximumBin()
    youdenClassifierMax = youden_histogram.GetXaxis().GetBinCenter(youdenClassifierBin)
    generalClassifierBin = general_histogram.GetMaximumBin()
    generalClassifierMax = general_histogram.GetXaxis().GetBinCenter(generalClassifierBin)
    youdenNhitMax = youden_histogram.GetMaximum()
    generalNhitMax = general_histogram.GetMaximum()

    youdenAlphaRejection = alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), youdenClassifierBin, youden_histogram.GetNbinsY()) / allAlphas
    youdenBetaAcceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, youdenClassifierBin) / allBetas
    generalAlphaRejection = alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), generalClassifierBin, youden_histogram.GetNbinsY()) / allAlphas
    generalBetaAcceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, generalClassifierBin) /allBetas

    return [youdenClassifierMax, youdenNhitMax, generalClassifierMax, generalNhitMax, youdenAlphaRejection, youdenBetaAcceptance, generalAlphaRejection, generalBetaAcceptance, meanNhit]
                                                                                  
# Collects user input from the command line for files and customization                                
parser = argparse.ArgumentParser()
parser.add_argument('--filetype', '-f', type = str, default = '', help = 'File type of alpha and beta files ("ntuple" OR "ratds")')
parser.add_argument('--alphafile', '-a', type = str, default = '', help = 'Alpha file to use (no ".root" extension)')
parser.add_argument('--betafile', '-b', type = str, default = '', help = 'Beta file to use (no ".root" extension)')
parser.add_argument('--shape', '-s', type = str, default = 'list', help = 'Shape for plot ("list" OR "square")')
parser.add_argument('--rhoCoordinates', '-r', type = int, nargs='+', help = 'List of rho coordinates')
parser.add_argument('--zCoordinates', '-z', type = int, nargs='+', help = 'List of z coordinates')
parser.add_argument('--distance', '-d', type = int, default = 1, help = 'Distance between coordinates')
parser.add_argument('--sideLength', '-l', type = int, help = 'Side length of square')
parser.add_argument('--ratio', '-ra', type = int, default = 1, help = 'Ratio of alpha to beta events')
parser.add_argument('--name', '-n', type = str, default = '', help = 'Name of files produced')

args = parser.parse_args()
alphaFile = "{}*.root".format(args.alphafile)
betaFile = "{}*.root".format(args.betafile)
ratio = args.ratio

r.gROOT.SetBatch(1) 

# Depending on the filetype the user inputted, opens the corresponding c++ file
if args.filetype == "ntuple":
    r.gROOT.LoadMacro("./NtupleValues.cpp+")
elif args.filetype == "ratds":
    r.gROOT.LoadMacro("./RATDSValues.cpp+")
else:
    raise Exception("Filetype must be either 'ntuple' OR 'ratds'")

rhoCoordinates = r.std.vector('double')()
zCoordinates = r.std.vector('double')()
distance = args.distance
name = args.name

# Creates rho and z coordinate lists for histogram segmentation based on user input
if args.shape == "square":
    squareLength = args.sideLength
    for i in range(-math.floor(squareLength/distance, math.floor(squareLength/distance))):
        for j in range(-math.floor(squareLength/distance, math.floor(squareLength/distance))):
            rhoCoordinates.push_back(i)
            zCoordinates.push_back(j)                

if args.shape == "list":
	for i in args.rhoCoordinates:
		rhoCoordinates.push_back(i)
	for i in args.zCoordinates:
		zCoordinates.push_back(i)
        
xTicks = []
yTicks = []

# Generates a list of x tick locations and y tick locations for histogram appearance
for i in range(0, math.floor(min(rhoCoordinates)-max(rhoCoordinates)/distance)):
    xTicks.append(min(rhoCoordinates) + (distance*i))
                                 
for i in range(0, math.floor(min(zCoordinates)-max(zCoordinates)/distance)):
    yTicks.append(min(zCoordinates) + (distance*i))

ClassifierYoudenArray = []
ValueYoudenArray = []
ClassifierGeneralArray = []
ValueGeneralArray = []

AlphaRejectionYoudenArray = []
BetaAcceptanceYoudenArray = []
AlphaRejectionGeneralArray = []
BetaAcceptanceGeneralArray = []

ClassifierAlphaArray = []
NhitAlphaArray = []
ClassifierBetaArray = []
NhitBetaArray = []

graphs = [ClassifierYoudenArray, ClassifierGeneralArray, ValueYoudenArray, ValueGeneralArray, AlphaRejectionYoudenArray, BetaAcceptanceYoudenArray, AlphaRejectionGeneralArray, BetaAcceptanceGeneralArray, ClassifierAlphaArray, NhitAlphaArray, ClassifierBetaArray, NhitBetaArray]
titles = ["Youden Classifier Cut Value","General Classifier Cut Value", "Youden Cut Value", "General Cut Value", r"Youden $\alpha$ Rejection", r"Youden $\beta$ Acceptance", r"General $\alpha$ Rejection", r"General $\beta$ Acceptance", r"Classification $\alpha$ Summary",r"$N_{\mathrm{hit}}$ $\alpha$ Summary", r"Classification $\beta$ Summary", r"$N_{\mathrm{hit}}$ $\beta$ Summary"]
graphs2 = ["YoudenClassifierCut", "GeneralClassifierCut", "YoudenCutValue", "GeneralCutValue", "YoudenAlphaRejection", "YoudenBetaAcceptance", "GeneralAlphaRejection", "GeneralBetaAcceptance", "ClassifierAlphaArray", "NhitAlphaArray", "ClassifierBetaArray", "NhitBetaArray"]
colorbar = ["Classifier Value", "Classifier Value", "Youden Statistic Value", "General Statistic Value", r"$\alpha$ Rejection", r"$\beta$ Acceptance", r"$\alpha$ Rejection", r"$\beta$ Acceptance", "Classifier Value", r"$N_{\mathrm{hit}}$ Value", "Classifier Value", r"$N_{\mathrm{hit}}$ Value"]

alphaHistograms = []
betaHistograms = []

# Based on filetype, calls a c++ function to sort the events from the given files into seperate histograms for every coordinate
if args.filetype == "ntuple":
    alphaHistograms = r.NhitHistograms(alphaFile, betaFile, np.array(rhoCoordinates), np.array(zCoordinates), ratio, distance, "alphaHist")
    betaHistograms = r.NhitHistograms(betaFile, betaFile, np.array(rhoCoordinates), np.array(zCoordinates), ratio, distance, "betaHist")
elif args.filetype == "ratds":
    alphaHistograms = r.NhitHistograms(alphaFile, betaFile, "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", rhoCoordinates, zCoordinates, distance, "alphaHist")
    betaHistograms = r.NhitHistograms(betaFile, betaFile, "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", rhoCoordinates, zCoordinates, distance, "betaHist")

for i in range(0, len(alphaHistograms)):
    values = rejectionInfo(alphaHistograms[i], betaHistograms[i], ratio)

    ClassifierYoudenArray.append(values[0])
    ValueYoudenArray.append(values[1])
    ClassifierGeneralArray.append(values[2])
    ValueGeneralArray.append(values[3])
    AlphaRejectionYoudenArray.append(values[4])
    BetaAcceptanceYoudenArray.append(values[5])
    AlphaRejectionGeneralArray.append(values[6])
    BetaAcceptanceGeneralArray.append(values[7])
                                          
    ClassifierAlphaArray.append(alphaHistograms[i].GetMean(2))
    NhitAlphaArray.append(alphaHistograms[i].GetMean(1))
    ClassifierBetaArray.append(betaHistograms[i].GetMean(2))
    NhitBetaArray.append(betaHistograms[i].GetMean(1))                              

#for graph in graphs:
#        graph.extend([-100, -100, -100])

#rhoCoordinates.push_back(3)
#rhoCoordinates.push_back(4)
#rhoCoordinates.push_back(4)
#zCoordinates.push_back(4)
#zCoordinates.push_back(4)
#zCoordinates.push_back(4)

for i in range(0,12):
    p.hist2d(rhoCoordinates, zCoordinates, bins=(5,4), range=((-.5,4.5),(1.5,4.5)), weights = graphs[i], cmap=p.cm.viridis, cmin=-10)

    for k in range(len(rhoCoordinates)):
        array = graphs[i]
        p.text(rhoCoordinates[k], zCoordinates[k], s.significantFigures(array[k], 4), ha="center",va="center",color="w")

    # Aesthetic customization
    p.xlabel(r"$\rho$ coordinate")
    p.ylabel(r"$z$ coordinate (m)")
    p.title(titles[i])
    p.xticks(xTicks)
    p.yticks(yTicks)
    p.colorbar(label=colorbar[i])
    p.show()               
    p.savefig("{}{}.pdf".format(name,graphs2[i]))
    p.clf()
