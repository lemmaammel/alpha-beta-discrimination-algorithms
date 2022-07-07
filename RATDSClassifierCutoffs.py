# compiles information about the classifier cutoff values across positions

#Useful generic python imports
from __future__ import print_function
from string import Template
import os, sys, time, math, csv, argparse, glob
import numpy as np
from array import array

#Imports particularly important for our purposes
import ROOT as r
import rat

loopvars = ["0_0_2","0_0_3", "0_0_4", "1_0_2","1_0_3","1_0_4","2_0_2","2_0_3","2_0_4","3_0_2","3_0_3","4_0_2"]

parser = argparse.ArgumentParser()
parser.add_argument('--alphafile', '-a', type = str, default = '', help = 'Alpha file to use')
parser.add_argument('--betafile', '-b', type = str, default = '', help = 'Beta file to use')
parser.add_argument('--rhoCoordinates', '-r', type = int, nargs='+', help = 'List of rho coordinates')
parser.add_argument('--zCoordinates', '-z', type = int, nargs='+', help = 'List of z coordinates')
parser.add_argument('--outfile', '-o', type = str, default = 'output.csv', help = 'Name of file to output')

args = parser.parse_args()

r.gROOT.SetBatch(1) 
r.gROOT.LoadMacro("/data/snoplus/home/ammel/rat/example/root/NhitHistogram.cpp+")

fields = ["location", "classifierCutoffGeneral"]
rows = []
filename = args.outfile
rhoCoordinates = args.rhoCoordinates
zCoordinates = args.zCoordinates

for i in range(0, len(rhoCoordinates)):
        
        name1 = args.alphafile 
        name2 = args.betafile
        
        alphaHistogram = r.NhitHistogram(args.alphafile, args.betafile, "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", "alphaHist", rhoCoordinates[i], zCoordinates[i], distance)
        betaHistogram = r.NhitHistogram(args.alphafile, args.betafile, "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood", "betaHist", rhoCoordinates[i], zCoordinates[i], distance)
        values = rejectionInfo(alphaHistogram, betaHistogram, ratio)
        
        currentRow = [loopvar, values[2]]
        rows.append(currentRow)

with open(filename, "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)

def rejectionInfo(alpha_hist, beta_hist, ratio):

    alpha_histogram = alpha_hist
    beta_histogram = beta_hist                                       

    #Like in the other file let's report a separate mean for alphas and betas
    meanNhit = (alpha_histogram.GetMean(1)+beta_histogram.GetMean(1))/(alpha_histogram.Integral()+beta_histogram.Integral())

    # cut selection histograms
    youden_histogram = TH1D("Youden's J Statistic", "Youden's J Statistic", alpha_histogram.GetNbinsY(),  alpha_histogram.GetMinimum(), alpha_histogram.GetMaximum())
    general_histogram = TH1D("General Cut Statistic", "General Cut Statistic", alpha_histogram.GetNbinsY(), alpha_histogram.GetMinimum(), alpha_histogram.GetMaximum())
                                           
    alpha_rejection = 0
    alpha_acceptance = 0
    beta_rejection = 0
    beta_acceptance = 0

    youden_statistic = 0
    general_statistic = 0

    for k in range(0, youden_histogram.GetNBinsX()):
        alpha_acceptance = ratio*alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, k)
        alpha_rejection = ratio*alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), k, youden_histogram.GetNbinsY())
        beta_acceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, k)
        beta_rejection = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), k, youdenSelection.GetNbinsY())

        if not (beta_acceptance == 0 and alpha_rejection == 0):
            youden_statistic = beta_acceptance/(beta_acceptance+beta_rejection) + alpha_rejection/(alpha_acceptance+alpha_rejection)
            general_statistic = beta_acceptance/sqrt(beta_acceptance+alpha_rejection)

        else:
            youden_statistic = 0
            general_statistic = 0

        #fill for cut selection stats
        youden_histogram.Fill(x, youden_statistic)
        general_histogram.Fill(x, general_statistic)

        x += abs(youden_histogram.GetXaxis().GetXmin() - youden_histogram.GetXaxis().GetXmax())/youden_histogram.getNbinsX()


    youdenClassifierBin = youden_histogram.GetMaximumBin()
    youdenClassifierMax = youden_histogram.GetXaxis().GetBinCenter(youdenClassifierBin)
    generalClassifierBin = general_histogram.GetMaximumBin()
    generalClassifierMax = general_histogram.GetXaxis().GetBinCenter(generalClassifierBin)
    youdenNhitMax = youden_histogram.GetMaximum()
    generalNhitMax = general_histogram.GetMaximum()

    allAlphas = alpha_histogram.Integral()
    allBetas = beta_histogram.Integral()

    if allAlphas == 0:
        allAlphas = 1e-15
   
    if allBetas == 0:
        allBetas = 1e-15

    youdenAlphaRejection = alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), youdenClassifierBin, youden_histogram.GetNbinsY()) / allAlphas
    youdenBetaAcceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, youdenClassifierBin) / allBetas
    generalAlphaRejection = alpha_histogram.Integral(1, youden_histogram.GetNbinsX(), generalClassifierBin, youden_histogram.GetNbinsY()) / allAlphas
    generalBetaAcceptance = beta_histogram.Integral(1, youden_histogram.GetNbinsX(), 1, generalClassifierBin) /allBetas

    return [youdenClassifierMax, youdenNhitMax, generalClassifierMax, generalNhitMax, youdenAlphaRejection, youdenBetaAcceptance, generalAlphaRejection, generalBetaAcceptance, meanNhit]
                                           
