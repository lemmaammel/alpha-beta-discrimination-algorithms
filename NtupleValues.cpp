// contains functions to generate scatterplots of real alpha and beta events, find the optimal cutoff values, and calculate the data relevant to those values
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <string>
#include <cmath>
#include <iostream>

bool isInRange(double posx, double posy, double posz, double rho, double z, double distance) {
    bool inRange = true;
    if (posz < z-distance*1000 || posz > z+distance*1000) {
        inRange = false;
    }
    double posrho = sqrt((posx*posx) + (posy*posy));
    if (posrho < rho-distance*1000 || posrho > rho+distance*1000) {
        inRange = false;
    }
    return inRange;
}

std::vector<TH2D*> NhitHistograms(std::string alphaFile, std::string betaFile, std::vector<double> rho, 
                    const std::string& type = "bothHists", std::vector<double> z, int style, double distance, 
                    const int nhitBins = 100, const double nhitMin = 100, const double nhitMax = 1200,
                    const int classBins = 100, const double classMin = -0.1, const double classMax = 0.1, 
                    const bool full = false) {

    std::vector<TH2D*> alphaHistograms;
    std::vector<TH2D*> betaHistograms;

    // create map of histograms
    std::map<std::string, std::vector<TH2D*>> histFileMap;

    for(size_t i = 0; i < rho.size(); i++) {
        alphaHistograms.push_back(new TH2D("#alpha Events", "N_{hit} Vs. Classification", nhitBins, nhitMin, nhitMax, classBins, classMin, classMax));
        betaHistograms.push_back(new TH2D("#beta Events", "N_{hit} Vs. Classification", nhitBins, nhitMin, nhitMax, classBins, classMin, classMax));
    }

    TCanvas *c2 = new TCanvas("c2", "Classification Histogram", 200, 10, 1300, 1300);
    c2->cd();

    if (type=="alphaHist" || type=="bothHists") {
        // customize aesthetic features and labels
            for(size_t i = 0; i < rho.size(); i++) {
                alphaHistograms[i]->SetMarkerColor(4);
                alphaHistograms[i]->SetMarkerStyle(20);
                alphaHistograms[i]->GetYaxis()->SetTitle("Classification Value/Number of Hits");
                alphaHistograms[i]->GetYaxis()->SetTitleOffset(1.2);
                alphaHistograms[i]->GetXaxis()->SetTitle("Number of Hits");
                alphaHistograms[i]->GetXaxis()->SetLabelSize(0.03);
                alphaHistograms[i]->GetYaxis()->SetLabelSize(0.02);
            }
            histFileMap[alphaFile] = alphaHistograms;
        }
            
        if(type=="betaHist" || type=="bothHists") {
        // customize aesthetic features and labels
            for(size_t i = 0; i < rho.size(); i++) {
                betaHistograms[i]->SetMarkerColor(6);
                betaHistograms[i]->SetMarkerStyle(20);
                betaHistograms[i]->GetYaxis()->SetTitle("Classification Value/Number of Hits");
                betaHistograms[i]->GetYaxis()->SetTitleOffset(1.2);
                betaHistograms[i]->GetXaxis()->SetTitle("Number of Hits");
                betaHistograms[i]->GetXaxis()->SetLabelSize(0.03);
                betaHistograms[i]->GetYaxis()->SetLabelSize(0.02);
            }
            histFileMap[betaFile] = betaHistograms;
        }

    for (std::map<std::string, TH2D*>::iterator it = histFileMap.begin(); it != histFileMap.end(); ++it) {
        TFile *f = TFile::Open(it->first.c_str());
        TTree *t = (TTree*)f->Get("output");

        double posx, posy, posz, berkeleyAlphaBeta;
        int nhits;
        bool fitValid;
        t->SetBranchAddress("posx", &posx);
        t->SetBranchAddress("posy", &posy);
        t->SetBranchAddress("posz", &posz);
        t->SetBranchAddress("berkeleyAlphaBeta", &berkeleyAlphaBeta);
        t->SetBranchAddress("nhits", &nhits);
        size_t k = t->GetEntries();

        for (size_t i = 0; i < k; i++) {
            t->GetEntry(i);
            if (nhits < 0) {
                continue;
            }
            for(size_t i = 0; i < rho.size(); i++) {
                if(full) {
                    it->second[i]->Fill(nhits, berkeleyAlphaBeta/nhits);
                }
                else if(isInRange(posx, posy, posz, rho[i], z[i], distance)) {
                    it->second[i]->Fill(nhits, berkeleyAlphaBeta/nhits);
                }
            }
        }

        f->Close();
    }

    // build a legend
    TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend->SetHeader("Legend");

    // draw relevant histograms on canvas and build legend
    if (alphaHist) {
        alphaHistograms[0]->Draw();
        legend->AddEntry(alphaHistograms[0], "#alpha Events", "p");
    }
    if (betaHist && !alphaHist) {
        betaHistograms[0]->Draw();
        legend->AddEntry(betaHistograms[0], "#beta Events", "p");
    }
    else if(bothHists) {
	    betaHistograms[0]->Draw("same");
	    legend->AddEntry(betaHistograms[0], "#beta Events", "p");
    }

    legend->Draw();

    // store histograms in PDFs
    c1->Print("realDataNhitHistogram.pdf", "pdf");

    // return histograms
    if (alphaHist || bothHists) {
        return alphaHistograms;
    }
    else {
        return betaHistograms;
    }
}

std::vector<double> averageValues(std::string filename) {
    TFile *f = TFile::Open(filename.c_str());
    TTree *t = (TTree*)f->Get("output");

    double posx, posy, posz;
    int nhits;
    bool fitValid;
    t->SetBranchAddress("posx", &posx);
    t->SetBranchAddress("posy", &posy);
    t->SetBranchAddress("posz", &posz);
    t->SetBranchAddress("nhits", &nhits);
    size_t k = t->GetEntries();

    double avgX = 0.0;
    double avgY = 0.0;
    double avgZ = 0.0;
    double avgRho = 0.0;

    for (size_t i = 0; i < k; i++) {
        t->GetEntry(i);

        if(nhits < 0) {
            continue;
        }
        double posrho = sqrt(posx*posx + posy*posy);
        //We should make the valuewe check against flexible
        if (posrho < 0.0 || posrho > 5.0) {
            continue;
        }
        if(posz < -5.0 || posz > 5.0) {
            continue;
        }

        avgX += posx;
        avgY += posy;
        avgZ += posz;
        avgRho += posrho;
    }

    f->Close();

    std::vector<double> values;
    values.push_back(avgX/k);
    values.push_back(avgY/k);
    values.push_back(avgZ/k);
    values.push_back(avgRho/k);

    return values;
}
