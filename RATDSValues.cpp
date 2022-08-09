// contains functions to generate scatterplots of simulated alpha and beta events, find the optimal cutoff values, and calculate the data relevant to those values

#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>

#include <TStyle.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>

#include <cmath>
#include <string>
#include <iostream>

enum histType { alphaHist = 1, betaHist = 2, bothHists = 3 };

bool isInRange(const double posx, const double posy, const double posz, const double rho, const double z, const double distance) {
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

// draw scatterplots of alpha and beta events (nhit and classification value as parameters) for an alpha and beta file
std::vector<TH2D*> NhitHistograms(const std::string& alphaFile, const std::string& betaFile, const std::string& fitName, 
                    const std::string& className, const std::string& classification, const std::string& type = "bothHists", 
                    std::vector<double> rho, std::vector<double> z, const double distance = 1.0,
                    const int nhitBins = 100, const double nhitMin = 100, const double nhitMax = 1000,
                    const int classBins = 100, const double classMin = -0.1, const double classMax = 0.1, const bool full = false) {

    std::vector<TH2D*> alphaHistograms;
    std::vector<TH2D*> betaHistograms;

    // create map of histograms
    std::map<std::string, std::vector<TH2D*>> histFileMap;

    for(size_t i = 0; i < rho.size(); i++) {
        alphaHistograms.push_back(new TH2D("#alpha Events", "N_{hit} Vs. Classification", nhitBins, nhitMin, nhitMax, classBins, classMin, classMax));
        betaHistograms.push_back(new TH2D("#beta Events", "N_{hit} Vs. Classification", nhitBins, nhitMin, nhitMax, classBins, classMin, classMax));

    }

    // create canvas to draw the histogram on
    TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
    c1->cd();

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


    RAT::DB::Get()->SetAirplaneModeStatus(true);
    for (std::map<std::string, TH2D*>::iterator it = histFileMap.begin(); it != histFileMap.end(); ++it) {
        RAT::DU::DSReader currentReader(it->first);
        //loop through all entries in filename
        for (size_t i = 0; i < currentReader.GetEntryCount(); i++) {
            const RAT::DS::Entry& rDS = currentReader.GetEntry(i);
            
            //loop through all events in entries
            for (size_t j = 0; j < rDS.GetEVCount(); j++) {
                //Ignore retriggers from residual detector light
                if (j > 0) { break;}
            
                //get the ev
                const RAT::DS::EV& rEV = rDS.GetEV(j);

                if (!rEV.ClassifierResultExists(className)) { continue;}
                if (!rEV.GetClassifierResult(className).GetValid()) { continue;}
            
                // classifier result
                RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
                double cValue = cResult.GetClassification(classification);

                // nHit value for calibrated hits
                const RAT::DS::EV& numberHits = rDS.GetEV(j);
            
                // test position
                if(!rEV.FitResultExists(fitName)) { continue;}
                if(!rEV.GetFitResult(fitName).GetVertex(0).ContainsPosition()) { continue;}
                if(!rEV.GetFitResult(fitName).GetVertex(0).ValidPosition()) { continue;}
                if(!rEV.GetFitResult(fitName).GetVertex(0).ContainsEnergy()) { continue;}
                if(!rEV.GetFitResult(fitName).GetVertex(0).ValidEnergy()) { continue;}

                RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
                RAT::DS::FitVertex fVertex = fResult.GetVertex(0);

                TVector3 pos = fVertex.GetPosition();

                for(size_t i = 0; i < rho.size(); i++) {
                    if(full) {
                        it->second[i]->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
                    }
                    else if(isInRange(pos.X(), pos.Y(), pos.Z(), rho[i], z[i], distance)) {
                        it->second[i]->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
                    }
                }   
            }
        }
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