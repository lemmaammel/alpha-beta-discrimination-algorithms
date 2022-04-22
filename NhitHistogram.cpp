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

#include <string>
#include <iostream>

// draw scatterplots of alpha and beta events (nhit and classification value as parameters) for an alpha and beta file
TCanvas* NhitHistogram(const std::string& alphaFile, const std::string& betaFile, const std::string& fitName, const std::string& className, const std::string& classification) {
	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* alphaHistogram = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);
	TH2D* betaHistogram = new TH2D("#alpha Events","N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);

    std::map<std::string, TH2D*> histFileMap;
    histFileMap[alphaFile] = alphaHistogram;
    hsitFileMap[betaFile] = betaHistogram;

	// customize aesthetic of histograms
	alphaHistogram->SetMarkerColor(4);
	betaHistogram->SetMarkerColor(6);
	alphaHistogram->SetMarkerStyle(20);
	betaHistogram->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);

	//loop through all entries in filename
	for (std::map<std::string, TH2D*>::iterator it = histFileMap.begin(); it != histFileMap.end(); ++it) {
		RAT::DU::DSReader currentReader(it->first);

		for (size_t i = 0; i < currentReader.GetEntryCount(); i++) {
			const RAT::DS::Entry& rDS = currentReader.GetEntry(i);

			//loop through all events in entries
			for (size_t j = 0; j < rDS.GetEVCount(); j++) {
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
                //Consider changes about whether to require valid position and energy
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

                //We should think about whether we can remove these events in a more consistent way
				// discard events from residual detector light
				if(numberHits.GetCalPMTs().GetAllCount() > 25) {
					it->second->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
				}
			}
		}
	}

	// customize histogram labeling
	alphaHistogram->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	alphaHistogram->GetYaxis()->SetTitleOffset(1.2);
	alphaHistogram->GetXaxis()->SetTitle("Number of Hits");
	alphaHistogram->GetXaxis()->SetLabelSize(0.03);
	alphaHistogram->GetYaxis()->SetLabelSize(0.02);

	alphaHistogram->Draw();
	betaHistogram->Draw("same");

	// build a legend
	TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);
	legend->SetHeader("Legend");
	legend->AddEntry(alphaHistogram, "#alpha Events", "p");
	legend->AddEntry(betaHistogram, "#beta Events", "p");
	legend->Draw();


	c1->Print("realDataNhitHistogram.pdf", "pdf");
	return c1;
}

// draw scatterplots of alpha events (nhit and classification value as parameters) for an alpha file
TH2D* NhitHistogramAlpha(const std::string& alphaFile, const std::string& fitName, const std::string& className, const std::string& classification) {
	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* alphaHistogram = new TH2D("#alpha Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);

	alphaHistogram->SetMarkerColor(4);
	alphaHistogram->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);
	RAT::DU::DSReader currentReader(alphaFile);

	for (size_t i = 0; i < currentReader.GetEntryCount(); i++) {
		const RAT::DS::Entry& rDS = currentReader.GetEntry(i);

		//loop through all events in entries
		for (size_t j = 0; j < rDS.GetEVCount(); j++) {
			//get the ev
			const RAT::DS::EV& rEV = rDS.GetEV(j);

			if (!rEV.ClassifierResultExists(className)) {
                continue;
            }
            if (!rEV.GetClassifierResult(className).GetValid()) {
                continue;
            }

			//classifier result
			RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
			double cValue = cResult.GetClassification(classification);

			//nHit value for calibrated hits
			const RAT::DS::EV& numberHits = rDS.GetEV(j);

			if (!rEV.FitResultExists(fitname)) {
                continue;
            }
            if (rEV.GetFitResult(fitname).GetVertexCount() != 0) {
                continue;
            }
            if (!rEV.GetFitResult(fitname).GetVertex(0).ContainsPosition()) {
                continue;
            }
            //Consider changes about whether to require valid position and energy
            if (!rEV.GetFitResult(fitname).GetVertex(0).ValidPosition()) {
                continue;
            }
            if (!rEV.GetFitResult(fitname).GetVertex(0).ContainsEnergy()) {
                continue;
            }
            if (!rEV.GetFitResult(fitname).GetVertex(0).ValidEnergy()) {
                continue;
            }
			// test position
			RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
			RAT::DS::FitVertex fVertex = fResult.GetVertex(0);

			if (numberHits.GetCalPMTs().GetAllCount() > 25) {
				alphaHistogram->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
			}
		}
	}

	alphaHistogram->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	alphaHistogram->GetYaxis()->SetTitleOffset(1.2);
	alphaHistogram->GetXaxis()->SetTitle("Number of Hits");
	alphaHistogram->GetXaxis()->SetLabelSize(.03);
	alphaHistogram->GetYaxis()->SetLabelSize(.02);

	alphaHistogram->Draw();
	alphaHistogram->Draw("same");

	// build a legend
	TLegend *legend = new TLegend(0.1,0.7,0.48,.9);
	legend->SetHeader("Legend");
	legend->AddEntry(alphaHistogram, "#alpha Events", "p");
	legend->Draw();

	c1->Print("realDataAlphaNhitHistogram.pdf", "pdf");

	return alphaEvents;
}

// draw scatterplots of beta events (nhit and classification value as parameters) for a beta file
TH2D* NhitHistogramBeta(const std::string& betaFile, const std::string& fitName, const std::string& className, const std::string& classification) {

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* betaHistogram = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);

	betaHistogram->SetMarkerColor(4);
	betaHistogram->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);
	RAT::DU::DSReader currentReader(betaFile);

	for (size_t i = 0; i < currentReader.GetEntryCount(); i++) {
		const RAT::DS::Entry& rDS = currentReader.GetEntry(i);

		//loop through all events in entries
		for (size_t j=0; j<rDS.GetEVCount(); j++) {
			//get the ev
			const RAT::DS::EV& rEV = rDS.GetEV(j);

			if (!rEV.ClassifierResultExists(className)) {
                continue;
            }
            if (!rEV.GetClassifierResult(className).GetValid()) {
                continue;
            }

			//classifier result
			RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
			double cValue = cResult.GetClassification(classification);

			//nHit value for calibrated hits
			const RAT::DS::EV& numberHits = rDS.GetEV(j);

			// test position
			if (!rEV.FitResultExists(fitname)) {
                continue;
            }
            if (rEV.GetFitResult(fitname).GetVertexCount() != 0) {
                continue;
            }
            if (!rEV.GetFitResult(fitname).GetVertex(0).ContainsPosition()) {
                continue;
            }
            //Consider changes about whether to require valid position and energy
            if (!rEV.GetFitResult(fitname).GetVertex(0).ValidPosition()) {
                continue;
            }
            if (!rEV.GetFitResult(fitname).GetVertex(0).ContainsEnergy()) {
                continue;
            }
            if (!rEV.GetFitResult(fitname).GetVertex(0).ValidEnergy()) {
                continue;
            }
			RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
			RAT::DS::FitVertex fVertex = fResult.GetVertex(0);

			if (numberHits.GetCalPMTs().GetAllCount() > 25) {
				betaHistogram->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
			}
		}
	}

	betaHistogram->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	betaHistogram->GetYaxis()->SetTitleOffset(1.2);
	betaHistogram->GetXaxis()->SetTitle("Number of Hits");
	betaHistogram->GetXaxis()->SetLabelSize(.03);
	betaHistogram->GetYaxis()->SetLabelSize(.02);

	betaHistogram->Draw();
	betaHistogram->Draw("same");

	// build a legend
	TLegend *legend = new TLegend(0.1,0.7,0.48,.9);
	legend->SetHeader("Legend");
	legend->AddEntry(betaHistogram, "#beta Events", "p");
	legend->Draw();

	c1->Print("realDataBetaNhitHistogram.pdf", "pdf");

	return betaEvents;
}

// plot the change in alpha rejection, beta acceptance, the beta sample fraction, Youden's J Statistic, and a general statistic relative to classification cutoff
TCanvas* rejectionHistogram(const std::string& alphaFile, const std::string& betaFile, const std::string& fitname, const std::string& classname, const std::string& classification, double ratio) {
	TH2D* analysisAlphaHistogram = NhitHistogramAlpha(alphaFile, fitname, classname, classification);
	TH2D* analysisBetaHistogram = NhitHistogramBeta(betaFile, fitname, classname, classification);

	TCanvas *c1 = new TCanvas("c1", "Rejection Histogram", 100, 10, 1300, 1300);
	c1->cd();

	TH1D* alphaRejectionHistogram = new TH1D("#alpha and #beta Analysis", "#alpha and #beta Analysis", 50, -0.1, 0.1);
	TH1D* betaAcceptanceHistogram = new TH1D("#beta Acceptance", "#beta Acceptance", 50, -0.1, 0.1);
	TH1D* betaSampleFraction = new TH1D("#beta Sample Fraction", "#beta Sample Fraction", 50, -0.1, 0.1);

	//cut selection histograms
	TH1D* youdenSelection = new TH1D("Youden's J Statistic", "Youden's J Statistic", 50, -0.1, 0.1);
	TH1D* generalSelection = new TH1D("General Cut Statistic", "General Cut Statistic", 50, -0.1, 0.1);

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

	double totalAlphaHits = ratio*analysisAlphaHistogram->Integral();
	double totalBetaHits = analysisBetaHistogram->Integral();
	double currentAlphaHits1 = 0;
	double currentAlphaHits2 = 0;
	double currentBetaHits = 0;
	double currentBetaHits2 = 0;
	//Where does this number come from?
    double currentX = -0.1;

	double youdenStatistic;
	double generalStatistic;

	alphaRejectionHistogram->SetBinContent(-0.098, 0);

	for (size_t k = 0;  k< 50; k++) {
		currentAlphaHits1 = ratio*analysisAlphaHistogram->Integral(1, 100, 2*k, 100);
		currentAlphaHits2 = ratio*analysisAlphaHistogram->Integral(1, 100, 1, 2*k);
		currentBetaHits = analysisBetaHistogram->Integral(1, 100, 1, 2*k);
		currentBetaHits2 = analysisBetaHistogram->Integral(1, 100, 2*k, 100);


		alphaRejectionHistogram->SetBinContent(k, currentAlphaHits1/totalAlphaHits);
		betaAcceptanceHistogram->SetBinContent(k, currentBetaHits/totalBetaHits);

		if(!(currentBetaHits==0 && currentAlphaHits2 == 0)) {
			betaSampleFraction->SetBinContent(k, currentBetaHits/(currentBetaHits+currentAlphaHits2));
			youdenStatistic = (currentBetaHits/(currentBetaHits+currentBetaHits2)) + (currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2))-1;
			generalStatistic = (currentBetaHits)/(sqrt(currentBetaHits+currentAlphaHits2)*50);
		}
		else {
			betaSampleFraction->SetBinContent(k, 1.0);
			youdenStatistic = 0;
			generalStatistic = 0;
		}

		//fill for cut selection stats
		youdenSelection->SetBinContent(k, youdenStatistic);
		generalSelection->SetBinContent(k, generalStatistic);

		currentX+=0.004;
	}

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

	return c1;
}

// calculate and return statistics about the optimal classifier cutoff and the corresponding acceptances and rejections
std::vector<double> rejectionInfo(const std::string& alphaFile, const std::string& betaFile, const std::string& fitname, const std::string& classname, const std::string& classification, double ratio) {

	TH2D* analysisAlphaHistogram = NhitHistogramAlpha(alphaFile, fitname, classname, classification);
	TH2D* analysisBetaHistogram = NhitHistogramBeta(betaFile, fitname, classname, classification);

	double meanNhit = (analysisBetaHistogram->GetMean(1)+analysisAlphaHistogram->GetMean(1))/(analysisBetaHistogram->Integral()+analysisAlphaHistogram->Integral());

	//cut selection histograms
	TH1D* youdenSelection = new TH1D("Youden's J Statistic", "Youden's J Statistic", 100, -0.1, 0.1);
	TH1D* generalSelection = new TH1D("General Cut Statistic", "General Cut Statistic", 100, -0.1, 0.1);

	double currentAlphaHits1 = 0;
	double currentAlphaHits2 = 0;
	double currentBetaHits = 0;
	double currentBetaHits2 = 0;
	double currentX = -0.098;

	double youdenStatistic;
	double generalStatistic;

	for (size_t k = 0; k < 100; k++) {
		currentAlphaHits1 = ratio*analysisAlphaHistogram->Integral(1, 100, k, 100);
		currentAlphaHits2 = ratio*analysisAlphaHistogram->Integral(1, 100, 1, k);
		currentBetaHits = analysisBetaHistogram->Integral(1, 100, 1, k);
		currentBetaHits2 = analysisBetaHistogram->Integral(1, 100, k, 100);

		if(!(currentBetaHits == 0 && currentAlphaHits2 == 0)) {
			youdenStatistic = currentBetaHits/(currentBetaHits+currentBetaHits2) + currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2);
			generalStatistic = currentBetaHits/sqrt(currentBetaHits+currentAlphaHits2);
		}
		else {
			youdenStatistic = 0;
			generalStatistic = 0;
		}

		//fill for cut selection stats
		youdenSelection->Fill(currentX, youdenStatistic);
		generalSelection->Fill(currentX, generalStatistic);

		currentX+= .002;
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
        allAlphas = 0.000001;
    }
    if (allBetas == 0) {
        allBetas = 0.000001;
    }

	double youdenAlphaRejection = analysisAlphaHistogram->Integral(1, 100, youdenClassifierBin, 100) / allAlphas;
	double youdenBetaAcceptance = analysisBetaHistogram->Integral(1, 100, 1, youdenClassifierBin) / allBetas;
	double generalAlphaRejection = analysisAlphaHistogram->Integral(1, 100, generalClassifierBin, 100) / allAlphas;
	double generalBetaAcceptance = analysisBetaHistogram->Integral(1, 100, 1, generalClassifierBin) /allBetas;

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
