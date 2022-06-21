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

parser = argparse.ArgumentParser()
parser.add_argument('--filetype', '-a', type = str, default = '', help = 'File type of alpha and beta files ("ntuple" OR "ratds")')
parser.add_argument('--alphafile', '-a', type = str, default = '', help = 'Alpha file to use (no ".root" extension)')
parser.add_argument('--betafile', '-b', type = str, default = '', help = 'Beta file to use (no ".root" extension)')
parser.add_argument('--shape', '-s', type = str, default = 'list', help = 'Shape for plot ("list" OR "square")')

args = parser.parse_args()

r.gROOT.SetBatch(1) 

if args.filetype == "ntuple":
    r.gROOT.LoadMacro("./NtupleValues.cpp+")
elif args.filetype == "ratds":
    r.gROOT.LoadMacro("./RATDSValues.cpp+")
else:
    raise Exception("Filetype must be either 'ntuple' OR 'ratds'")

rhoCoordinates = []
zCoordinates = []
distance = 1

#We discussed some ways to reimplement this on Thursday
if args.shape == "square"
        squareLength = int(raw_input("Please enter the side length of your square:"))
        distance = int(raw_input("Please enter the distance between the coordinates:"))
        for i in range(-math.floor(squareLength/distance, math.floor(squareLength/distance)
                for j in range(-math.floor(squareLength/distance, math.floor(squareLength/distance)
                        rhoCoordinates.extend(i)
                        zCoordinates.extend(j)
                       
if args.shape == "list"
        rhoCoordinates = list(map(int, raw_input("Please enter the rho coordinates in the form '3 2 3 4 8':").split()))
        zCoordinates = list(map(int, raw_input("Please enter the z coordinates in the form '4 6 3 8 8':").split()))
        distance = int(raw_input("Please enter the distance between the coordinates:"))
                                 
xTicks = []
yTicks = []

for i in range(0, math.floor(min(rhoCoordinates)-max(rhoCoordinates)/distance))
        xTicks.append(min(rhoCoordinates) + (distance*i))
                                 
for i in range(0, math.floor(min(zCoordinates)-max(zCoordinates)/distance))
        yTicks.append(min(zCoordinates) + (distance*i))
                

#We can make the ratio an argument with argparse
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


graphs = [ClassifierYoudenArray, ClassifierGeneralArray, ValueYoudenArray, ValueGeneralArray, AlphaRejectionYoudenArray, BetaAcceptanceYoudenArray, AlphaRejectionGeneralArray, BetaAcceptanceGeneralArray, ClassifierAlphaArray, NhitAlphaArray, ClassifierBetaArray, NhitBetaArray]
titles = ["Youden Classifier Cut Value","General Classifier Cut Value", "Youden Cut Value", "General Cut Value", r"Youden $\alpha$ Rejection", r"Youden $\beta$ Acceptance", r"General $\alpha$ Rejection", r"General $\beta$ Acceptance", r"Classification $\alpha$ Summary",r"$N_{\mathrm{hit}}$ $\alpha$ Summary", r"Classification $\beta$ Summary", r"$N_{\mathrm{hit}}$ $\beta$ Summary"]
graphs2 = ["YoudenClassifierCut", "GeneralClassifierCut", "YoudenCutValue", "GeneralCutValue", "YoudenAlphaRejection", "YoudenBetaAcceptance", "GeneralAlphaRejection", "GeneralBetaAcceptance", "ClassifierAlphaArray", "NhitAlphaArray", "ClassifierBetaArray", "NhitBetaArray"]
colorbar = ["Classifier Value", "Classifier Value", "Youden Statistic Value", "General Statistic Value", r"$\alpha$ Rejection", r"$\beta$ Acceptance", r"$\alpha$ Rejection", r"$\beta$ Acceptance", "Classifier Value", r"$N_{\mathrm{hit}}$ Value", "Classifier Value", r"$N_{\mathrm{hit}}$ Value"]

for i in range(0, len(rhoCoordinates)):

        if args.filetype == "ntuple"
                alphaHistogram = r.NhitHistogram("{}*.root".format(args.alphafile), rhoCoordinates[i], zCoordinates[i], ratio, distance)
                betaHistogram = r.NhitHistogram("{}*.root".format(args.betafile), rhoCoordinates[i], zCoordinates[i], ratio, distance)
                values = rejectionInfo(alphaHistogram, betaHistogram, ratio)

        
        else if args.filetype == "ratds"
                                           alphaHist, rhoCoordinate, zCoordinate, distance
                alphaHistogram = r.NhitHistogram("{}*.root".format(args.alphafile), "{}*.root".format(args.betafile), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", alphaHist, rhoCoordinates[i], zCoordinates[i], distance)
                betaHistogram = r.NhitHistogram("{}*.root".format(args.alphafile), "{}*.root".format(args.betafile), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", betaHist, rhoCoordinates[i], zCoordinates[i], distance)
                values = = rejectionInfo(alphaHistogram, betaHistogram, ratio)


        ClassifierYoudenArray.append(values[0])
        ValueYoudenArray.append(values[1])
        ClassifierGeneralArray.append(values[2])
        ValueGeneralArray.append(values[3])
        AlphaRejectionYoudenArray.append(values[4])
        BetaAcceptanceYoudenArray.append(values[5])
        AlphaRejectionGeneralArray.append(values[6])
        BetaAcceptanceGeneralArray.append(values[7])
                                          
        ClassifierAlphaArray.append(alphaHistogram.GetMean(2))
        NhitAlphaArray.append(alphaHistogram.GetMean(1))
        ClassifierBetaArray.append(betaHistogram.GetMean(2))
        NhitBetaArray.append(betaHistogram.GetMean(1))                              

for graph in graphs:
        graph.extend([-100, -100, -100])

rhoCoordinates.extend(3, 4, 4)
zCoordinates.extend([4, 4, 4])

for i in range(0,12):
        p.hist2d(rhoCoordinates, zCoordinates, bins=(5,4), range=((-.5,4.5),(1.5,4.5)), weights = graphs[i], cmap=p.cm.viridis, cmin=-10)

        for k in range(len(rhoCoordinates)):
                array = graphs[i]
                p.text(rhoCoordinates[k], zCoordinates[k], s.significantFigures(array[k], 4), ha="center",va="center",color="w")

        p.xlabel(r"$\rho$ coordinate")
        p.ylabel(r"$z$ coordinate (m)")
        p.title(titles[i])
        p.xticks(xTicks)
        p.yticks(yTicks)
        p.colorbar(label=colorbar[i])
        p.show()               
        p.savefig("partialFillPDFs/SummaryPartialFill{}.pdf".format(graphs2[i]))
        p.clf()

def rejectionInfo(alpha_hist, beta_hist, ratio)
                                           
    alpha_histogram = alpha_hist
    beta_histogram = beta_hist                                       
                                           
    meanNhit = (alpha_histogram.GetMean(1)+beta_histogram.GetMean(1))/(alpha_histogram.Integral()+beta_histogram.Integral())

    # cut selection histograms
    youden_histogram = TH1D("Youden's J Statistic", "Youden's J Statistic", alpha_histogram.GetNbinsY(),  alpha_histogram.GetMinimum(), alpha_histogram.GetMaximum())
    general_histogram = TH1D("General Cut Statistic", "General Cut Statistic", alpha_histogram.GetNbinsY(), alpha_histogram.GetMinimum(), alpha_histogram.GetMaximum())
                                           
    alpha_rejection = 0
    alpha_acceptance = 0
    beta_rejection = 0
    beta_acceptance = 0
    x = general_histogram.GetXaxis().GetXmin()

    youden_statistic = 0
    general_statistic = 0

    for k in range(0, youden_histogram.GetNBinsX()):
        alpha_rejection = ratio*alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), k, youden_histogram.GetNbinsY())
        alpha_acceptance = ratio*alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, k)
        beta_acceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, k)
        beta_rejection = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), k, youdenSelection.GetNbinsY())

        if !(beta_acceptance == 0 && alpha_rejection == 0)
            youden_statistic = beta_acceptance/(beta_acceptance+beta_rejection) + alpha_rejection/(alpha_acceptance+alpha_rejection)
            general_statistic = beta_acceptance/sqrt(beta_acceptance+alpha_rejection)

        else
            youden_statistic = 0
            general__statistic = 0

        #fill for cut selection stats
        youden_histogram.Fill(x, youden_statistic)
        general_histogram.Fill(x, general_statistic)

        x += abs(youden_histogram.GetXaxis().GetXmin() - youden_histogram.GetXaxis().GetXmax())/youden_histogram.getNbinsX()
    }

    youdenClassifierBin = youden_histogram.GetMaximumBin()
    youdenClassifierMax = youden_histogram.GetXaxis().GetBinCenter(youdenClassifierBin)
    generalClassifierBin = general_histogram.GetMaximumBin()
    generalClassifierMax = general_histogram.GetXaxis().GetBinCenter(generalClassifierBin)
    youdenNhitMax = youden_histogram.GetMaximum()
    generalNhitMax = general_histogram.GetMaximum()

    allAlphas = alpha_histogram.Integral()
    allBetas = beta_histogram.Integral()

    if allAlphas==0
        allAlphas = 1e-15
   
    if allBetas == 0
        allBetas = 1e-15

    youdenAlphaRejection = alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), youdenClassifierBin, youden_histogram.GetNbinsY()) / allAlphas
    youdenBetaAcceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, youdenClassifierBin) / allBetas
    generalAlphaRejection = alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), generalClassifierBin, youden_histogram.GetNbinsY()) / allAlphas
    generalBetaAcceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, generalClassifierBin) /allBetas

    return [youdenClassifierMax, youdenNhitMax, generalClassifierMax, generalNhitMax, youdenAlphaRejection, youdenBetaAcceptance, generalAlphaRejection, generalBetaAcceptance, meanNhit]
                                           
                                           
                                           
                                           
