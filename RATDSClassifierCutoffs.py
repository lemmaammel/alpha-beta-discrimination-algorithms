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
import AlphaBetaDiscriminationHistograms as a

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
        values = a.rejectionInfo(alphaHistogram, betaHistogram, ratio)
        
        currentRow = [loopvar, values[2]]
        rows.append(currentRow)

with open(filename, "w") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(fields)
        csvwriter.writerows(rows)
