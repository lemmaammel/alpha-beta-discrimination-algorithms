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

bool getEventCoordinates(const double posx, const double posy, const double posz, const double rho, const double z, const double distance) {
    bool inRange = true;
    //Should 1000 be hardcoded?
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
TH2D* NhitHistogram(const std::string& alphaFile, const std::string& betaFile, const std::string& fitName, const std::string& className,
                    const std::string& classification, const histType type = bothHists,
                    const double rho = 100.0, const double z = 100.0, const double distance = 1.0,
                    const int nhitBins = 100, const double nhitMin = 100, const double nhitMax = 1000,
                    const int classBins = 100, const double classMin = -0.1, const double classMax = 0.1, const bool full = false) {

    // create canvas to draw the histogram on
    TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
    c1->cd();

    // create map of histograms
    std::map<std::string, TH2D*> histFileMap;

    if (alphaHist || bothHists) {
        // create histogram for alpha events
        TH2D* alphaHistogram = new TH2D("#alpha Events", "N_{hit} Vs. Classification", nhitBins, nhitMin, nhitMax, classBins, classMin, classMax);
        // customize aesthetic features and labels
        alphaHistogram->SetMarkerColor(4);
        alphaHistogram->SetMarkerStyle(20);
        alphaHistogram->GetYaxis()->SetTitle("Classification Value/Number of Hits");
        alphaHistogram->GetYaxis()->SetTitleOffset(1.2);
        alphaHistogram->GetXaxis()->SetTitle("Number of Hits");
        alphaHistogram->GetXaxis()->SetLabelSize(0.03);
        alphaHistogram->GetYaxis()->SetLabelSize(0.02);
        // add histogram to the map
        histFileMap[alphaFile] = alphaHistogram;
    }

    if (betaHist || bothHists) {
        // create histogram for beta events
        TH2D* betaHistogram = new TH2D("#beta Events","N_{hit} Vs. Classification", nhitBins, nhitMin, nhitMax, classBins, classMin, classMax);
        // customize aesthetic features and labels
        betaHistogram->SetMarkerColor(6);
        betaHistogram->SetMarkerStyle(20);
        betaHistogram->GetYaxis()->SetTitle("Classification Value/Number of Hits");
        betaHistogram->GetYaxis()->SetTitleOffset(1.2);
        betaHistogram->GetXaxis()->SetTitle("Number of Hits");
        betaHistogram->GetXaxis()->SetLabelSize(0.03);
        betaHistogram->GetYaxis()->SetLabelSize(0.02);
        // add histogram to the map
        histFileMap[betaFile] = betaHistogram;
    }

    RAT::DB::Get()->SetAirplaneModeStatus(true);

    //loop through all entries in filename
    for (std::map<std::string, TH2D*>::iterator it = histFileMap.begin(); it != histFileMap.end(); ++it) {
        RAT::DU::DSReader currentReader(it->first);

        for (size_t i = 0; i < currentReader.GetEntryCount(); i++) {
            const RAT::DS::Entry& rDS = currentReader.GetEntry(i);

            //loop through all events in entries
            for (size_t j = 0; j < rDS.GetEVCount(); j++) {
                //Ignore retriggers from residual detector light
                if (j > 0) {
                    break;
                }
                //get the ev
                const RAT::DS::EV& rEV = rDS.GetEV(j);

                if (!rEV.ClassifierResultExists(className)) {
                    continue;
                }
                if (!rEV.GetClassifierResult(className).GetValid()) {
                    continue;
                }

                // classifier result
                RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
                double cValue = cResult.GetClassification(classification);

                // nHit value for calibrated hits
                const RAT::DS::EV& numberHits = rDS.GetEV(j);

                // test position
                if(!rEV.FitResultExists(fitname)) {
                    continue;
                }
                if(rEV.GetFitResult(fitname).GetVertexCount() != 0) {
                    continue;
                }
                if(!rEV.GetFitResult(fitname).GetVertex(0).ContainsPosition()) {
                    continue;
                }

                // consider changes about whether to require valid position and energy
                if(!rEV.GetFitResult(fitname).GetVertex(0).ValidPosition()) {
                    continue;
                }
                if(!rEV.GetFitResult(fitname).GetVertex(0).ContainsEnergy()) {
                    continue;
                }
                if(!rEV.GetFitResult(fitname).GetVertex(0).ValidEnergy()) {
                    continue;
                }
                RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
                RAT::DS::FitVertex fVertex = fResult.GetVertex(0);

                TVector3 pos = fVertex.GetPosition();

                if(full) {
                    it->second->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
                }
                else if(getEventCoordinates(sqrt(pos.X()*pos.X() + pos.Y()*pos.Y()), pos.Z(), rho, z, distance) {
                    it->second->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
                }
            }
        }
    }

    // build a legend
    TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
    legend->SetHeader("Legend");

    // draw relevant histograms on canvas and build legend
    if (alphaHist || bothHists) {
        alphaHistogram->Draw();
        legend->AddEntry(alphaHistogram, "#alpha Events", "p");
    }
    else if (betaHist || bothHists) {
        betaaHistogram->Draw("same");
        legend->AddEntry(betaHistogram, "#beta Events", "p");
    }

    legend->Draw();

    // store histograms in PDFs
    c1->Print("realDataNhitHistogram.pdf", "pdf");

    // return histograms
    if (alphaHist || bothHists) {
        return alphaHistogram;
    }
    else {
        return betaHistogram;
    }
}

// calculate and return statistics about the optimal classifier cutoff and the corresponding acceptances and rejections
std::vector<double> rejectionInfo(const std::string& alphaFile, const std::string& betaFile, const std::string& fitname,
                                  const std::string& classname, const std::string& classification,
                                  const double ratio, const double rhoCoordinate, const double zCoordinate, const double distance, const bool printHistogram = false) {

    TH2D* analysisAlphaHistogram = NhitHistogram(alphaFile, betaFile, fitname, classname, classification, alphaHist, rhoCoordinate, zCoordinate, distance);
    TH2D* analysisBetaHistogram = NhitHistogram(alphaFile, betaFile, fitname, classname, classification, betaHist, rhoCoordinate, zCoordinate, distance);
    
    TCanvas *c1 = new TCanvas("c1", "Rejection Histogram", 100, 10, 1300, 1300);
    c1->cd();

    double meanNhit = (analysisBetaHistogram->GetMean(1)+analysisAlphaHistogram->GetMean(1))/(analysisBetaHistogram->Integral()+analysisAlphaHistogram->Integral());

    //cut selection histograms
    TH1D* alphaRejectionHistogram = new TH1D("#alpha and #beta Analysis", "#alpha and #beta Analysis", analysisAlphaHistogram.GetNbinsY(),  analysisAlphaHistogram.GetMinimum(), analysisAlphaHistogram.GetMaximum());
    TH1D* betaAcceptanceHistogram = new TH1D("#beta Acceptance", "#beta Acceptance", analysisAlphaHistogram.GetNbinsY(),  analysisAlphaHistogram.GetMinimum(), analysisAlphaHistogram.GetMaximum());
    TH1D* betaSampleFraction = new TH1D("#beta Sample Fraction", "#beta Sample Fraction", analysisAlphaHistogram.GetNbinsY(),  analysisAlphaHistogram.GetMinimum(), analysisAlphaHistogram.GetMaximum());
    
    TH1D* youdenSelection = new TH1D("Youden's J Statistic", "Youden's J Statistic", analysisAlphaHistogram.GetNbinsY(),  analysisAlphaHistogram.GetMinimum(), analysisAlphaHistogram.GetMaximum());
    TH1D* generalSelection = new TH1D("General Cut Statistic", "General Cut Statistic", analysisAlphaHistogram.GetNbinsY(), analysisAlphaHistogram.GetMinimum(), analysisAlphaHistogram.GetMaximum());
    
    double totalAlphaHits = ratio*analysisAlphaHistogram->Integral();
    double totalBetaHits = analysisBetaHistogram->Integral();
    double currentAlphaHits1 = 0;
    double currentAlphaHits2 = 0;
    double currentBetaHits = 0;
    double currentBetaHits2 = 0;
    double currentX = generalSelection.GetXaxis()->GetXmin();

    double youdenStatistic;
    double generalStatistic;

    for (size_t k = 0; k < youdenSelection.GetNbinsX(); k++) {
        currentAlphaHits1 = ratio*analysisAlphaHistogram->Integral(1, youdenSelection.GetNbinsX(), k, youdenSelection.GetNbinsY());
        currentAlphaHits2 = ratio*analysisAlphaHistogram->Integral(1, youdenSelection.GetNbinsX(), 1, k);
        currentBetaHits = analysisBetaHistogram->Integral(1, youdenSelection.GetNbinsX(), 1, k);
        currentBetaHits2 = analysisBetaHistogram->Integral(1, youdenSelection.GetNbinsX(), k, youdenSelection.GetNbinsY());
        
        alphaRejectionHistogram->SetBinContent(k, currentAlphaHits1/totalAlphaHits);
        betaAcceptanceHistogram->SetBinContent(k, currentBetaHits/totalBetaHits);

        if(!(currentBetaHits == 0 && currentAlphaHits2 == 0)) {
            betaSampleFraction->SetBinContent(k, currentBetaHits/(currentBetaHits+currentAlphaHits2));
            youdenStatistic = currentBetaHits/(currentBetaHits+currentBetaHits2) + currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2);
            generalStatistic = currentBetaHits/sqrt(currentBetaHits+currentAlphaHits2);
        }
        else {
            betaSampleFraction->SetBinContent(k, 1.0);
            youdenStatistic = 0;
            generalStatistic = 0;
        }

        //fill for cut selection stats
        youdenSelection->Fill(currentX, youdenStatistic);
        generalSelection->Fill(currentX, generalStatistic);

        currentX+= (std::abs(youdenSelection.GetXaxis()->GetXmin() - youdenSelection.GetXaxis()->GetXmax()))/youdenSelection.getNbinsX();
    }
    
    if(printHistogram) {
        
         // customize aesthetics for histograms
        alphaRejectionHistogram->SetLineColor(4);
        betaAcceptanceHistogram->SetLineColor(6);
        betaSampleFraction->SetLineColor(3);
        youdenSelection->SetLineColor(2);
        generalSelection->SetLineColor(7);

        alphaRejectionHistogram->SetLineWidth(3);
        betaAcceptanceHistogram->SetLineWidth(3);
        betaSampleFraction->SetLineWidth(3);
        youdenSelection->SetLineWidth(3);
        generalSelection->SetLineWidth(3);
        
        alphaRejectionHistogram->GetYaxis()->SetTitle("#alpha Rejection / #beta Acceptance (%)");
        alphaRejectionHistogram->GetYaxis()->SetTitleOffset(1.2);
        alphaRejectionHistogram->GetXaxis()->SetTitle("Classification Cutoff");
        alphaRejectionHistogram->GetXaxis()->SetLabelSize(0.03);
        alphaRejectionHistogram->GetYaxis()->SetLabelSize(0.02);

        alphaRejectionHistogram->Draw();
        betaAcceptanceHistogram->Draw("same");
        betaSampleFraction->Draw("same");
        youdenSelection->Draw("same");
        generalSelection->Draw("same");

        // build a legend
        TStyle* gStyle = new TStyle("", "");

        TLegend *legend = new TLegend(0.1,0.7,0.3,0.9);
        legend->SetHeader("Legend");
        legend->AddEntry(alphaRejectionHistogram, "#alpha Rejection", "l");
        legend->AddEntry(betaAcceptanceHistogram, "#beta Acceptance", "l");
        legend->AddEntry(betaSampleFraction, "#beta Sample Fraction", "l");
        legend->AddEntry(youdenSelection, "Youden's J Statistic", "l");
        legend->AddEntry(generalSelection, "General Selection Statistic", "l");
        legend->Draw();

        c1->Print("AlphaRejectionHistogram_Revised.pdf", "pdf");
    }

    double youdenClassifierBin = youdenSelection->GetMaximumBin();
    double youdenClassifierMax = youdenSelection->GetXaxis()->GetBinCenter(youdenClassifierBin);
    double generalClassifierBin = generalSelection->GetMaximumBin();
    double generalClassifierMax = generalSelection->GetXaxis()->GetBinCenter(generalClassifierBin);
    double youdenNhitMax = youdenSelection->GetMaximum();
    double generalNhitMax = generalSelection->GetMaximum();

    double allAlphas = analysisAlphaHistogram->Integral();
    double allBetas = analysisBetaHistogram->Integral();

    if (allAlphas==0) {
        allAlphas = 1e-15;
    }
    if (allBetas == 0) {
        allBetas = 1e-15;
    }

    double youdenAlphaRejection = analysisAlphaHistogram->Integral(1, youdenSelection.GetNbinsX(), youdenClassifierBin, youdenSelection.GetNbinsY()) / allAlphas;
    double youdenBetaAcceptance = analysisBetaHistogram->Integral(1, youdenSelection.GetNbinsX(), 1, youdenClassifierBin) / allBetas;
    double generalAlphaRejection = analysisAlphaHistogram->Integral(1, youdenSelection.GetNbinsX(), generalClassifierBin, youdenSelection.GetNbinsY()) / allAlphas;
    double generalBetaAcceptance = analysisBetaHistogram->Integral(1, youdenSelection.GetNbinsX(), 1, generalClassifierBin) /allBetas;

    std::vector<double> values;
    values.push_back(youdenClassifierMax);
    values.push_back(youdenNhitMax);
    values.push_back(generalClassifierMax);
    values.push_back(generalNhitMax);
    values.push_back(youdenAlphaRejection);
    values.push_back(youdenBetaAcceptance);
    values.push_back(generalAlphaRejection);
    values.push_back(generalBetaAcceptance);
    values.push_back(meanNhit);

    return values;
}
