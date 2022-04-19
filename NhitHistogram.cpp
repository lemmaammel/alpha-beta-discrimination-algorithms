#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DB.hh>
#include <RAT/DU/Utility.hh>
#include <TStreamerInfo.h>

#include <RAT/DS/MC.hh>
#include <TStyle.h>
#include <iostream>

#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
using namespace std;

// draw scatterplots of alpha and beta events (nhit and classification value as parameters) for an alpha and beta file
TCanvas* NhitHistogram(const string& alphaFile, const string& betaFile, const string& fitName, const string& className, const string& classification)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* alphaHistogram = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);
	TH2D* betaHistogram = new TH2D("#alpha Events","N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);
	TH2D* histograms [2] = {alphaHistogram, betaHistogram};
	string files[2] = {alphaFile, betaFile};

	// customize aesthetic of histograms
	alphaHistogram->SetMarkerColor(4);
	betaHistogram->SetMarkerColor(6);
	alphaHistogram->SetMarkerStyle(20);
	betaHistogram->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);

	//loop through all entries in filename
	for(int m = 0; m<2; m++)
	{
		RAT::DU::DSReader currentReader(files[m]);
		
		for(size_t i=0; i<currentReader.GetEntryCount(); i++)
		{
			const RAT::DS::Entry& rDS = currentReader.GetEntry(i);
			
			//loop through all events in entries
			for(size_t j=0; j<rDS.GetEVCount(); j++)
			{
				
				//get the ev
				const RAT::DS::EV& rEV = rDS.GetEV(j);
	
				if(!rEV.ClassifierResultExists(className)) continue;
				if(!rEV.GetClassifierResult(className).GetValid()) continue;

				// classifier result
				RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
				double cValue = cResult.GetClassification(classification);

				// nHit value for calibrated hits
				const RAT::DS::EV& numberHits = rDS.GetEV(j);
				
				// test position
				RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
				RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
				double x = fVertex.GetPosition().Z();

				// discard events from residual detector light
				if(numberHits.GetCalPMTs().GetAllCount() > 25)
				{
					histograms[m]->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
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
TH2D* NhitHistogramAlpha(const string& alphaFile, const string& fitName, const string& className, const string& classification)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* alphaHistogram = new TH2D("#alpha Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);

	alphaHistogram->SetMarkerColor(4);
	alphaHistogram->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);	
	RAT::DU::DSReader currentReader(alphaFile);

	for(size_t i=0; i<currentReader.GetEntryCount(); i++)
	{
		const RAT::DS::Entry& rDS = currentReader.GetEntry(i);
			
		//loop through all events in entries
		for(size_t j=0; j<rDS.GetEVCount(); j++)
		{
				
			//get the ev
			const RAT::DS::EV& rEV = rDS.GetEV(j);
	
			if(!rEV.ClassifierResultExists(className)) continue;
			if(!rEV.GetClassifierResult(className).GetValid()) continue;

			//classifier result
			RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
			double cValue = cResult.GetClassification(classification);

			//nHit value for calibrated hits
			const RAT::DS::EV& numberHits = rDS.GetEV(j);
				
			// test position
			RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
			RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
			double x = fVertex.GetPosition().Z();

			if(numberHits.GetCalPMTs().GetAllCount() > 25)
			{
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
TH2D* NhitHistogramBeta(const string& betaFile, const string& fitName, const string& className, const string& classification)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* betaHistogram = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -0.1, 0.1);

	betaHistogram->SetMarkerColor(4);
	betaHistogram->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);	
	RAT::DU::DSReader currentReader(betaFile);

	for(size_t i=0; i<currentReader.GetEntryCount(); i++)
	{
		const RAT::DS::Entry& rDS = currentReader.GetEntry(i);
			
		//loop through all events in entries
		for(size_t j=0; j<rDS.GetEVCount(); j++)
		{
				
			//get the ev
			const RAT::DS::EV& rEV = rDS.GetEV(j);
	
			if(!rEV.ClassifierResultExists(className)) continue;
			if(!rEV.GetClassifierResult(className).GetValid()) continue;

			//classifier result
			RAT::DS::ClassifierResult cResult = rEV.GetClassifierResult(className);
			double cValue = cResult.GetClassification(classification);

			//nHit value for calibrated hits
			const RAT::DS::EV& numberHits = rDS.GetEV(j);
				
			// test position
			RAT::DS::FitResult fResult = rEV.GetFitResult(fitName);
			RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
			double x = fVertex.GetPosition().Z();

			if(numberHits.GetCalPMTs().GetAllCount() > 25)
			{
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
TCanvas* rejectionHistogram(const string& alphaFile, const string& betaFile, const string& fitname, const string& classname, const string& classification, float ratio)
{

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

	float totalAlphaHits = ratio*analysisAlphaHistogram->Integral(1, 100, 1, 100);
	float totalBetaHits = analysisBetaHistogram->Integral(1, 100, 1, 100);
	float currentAlphaHits1 = 0;
	float currentAlphaHits2 = 0;
	float currentBetaHits = 0;
	float currentBetaHits2 = 0;
	float currentX = -0.1;

	float youdenStatistic;
	float generalStatistic;

	alphaRejectionHistogram->SetBinContent(-0.098, 0);

	for(int k = 0; k<50; k++)
	{
		currentAlphaHits1 = ratio*analysisAlphaHistogram->Integral(1, 100, 2*k, 100);
		currentAlphaHits2 = ratio*analysisAlphaHistogram->Integral(1, 100, 1, 2*k);
		currentBetaHits = analysisBetaHistogram->Integral(1, 100, 1, 2*k);
		currentBetaHits2 = analysisBetaHistogram->Integral(1, 100, 2*k, 100);


		alphaRejectionHistogram->SetBinContent(k, currentAlphaHits1/totalAlphaHits);
		betaAcceptanceHistogram->SetBinContent(k, currentBetaHits/totalBetaHits);

		if(not(currentBetaHits==0 & currentAlphaHits2==0))
		{	
			betaSampleFraction->SetBinContent(k, currentBetaHits/(currentBetaHits+currentAlphaHits2));
			youdenStatistic = (currentBetaHits/(currentBetaHits+currentBetaHits2)) + (currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2))-1;
			generalStatistic = (currentBetaHits)/(sqrt(currentBetaHits+currentAlphaHits2)*50);
		}
		else 
		{
			betaSampleFraction->SetBinContent(k, 1.0);
			youdenStatistic = 0;
			generalStatistic = 0;
		}		

		//fill for cut selection stats
		youdenSelection->SetBinContent(k, youdenStatistic);
		generalSelection->SetBinContent(k, generalStatistic);

		currentX+=0.004;
		cout << currentX;
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
float* rejectionInfo(const string& alphaFile, const string& betaFile, const string& fitname, const string& classname, const string& classification, float ratio)
{

	TH2D* analysisAlphaHistogram = NhitHistogramAlpha(alphaFile, fitname, classname, classification);

	TH2D* analysisBetaHistogram = NhitHistogramBeta(betaFile, fitname, classname, classification);

	float meanNhit = (analysisBetaHistogram->GetMean(1)+analysisAlphaHistogram->GetMean(1))/2;

	//cut selection histograms
	TH1D* youdenSelection = new TH1D("Youden's J Statistic", "Youden's J Statistic", 100, -0.1, 0.1);
	TH1D* generalSelection = new TH1D("General Cut Statistic", "General Cut Statistic", 100, -0.1, 0.1); 
	
	float currentAlphaHits1 = 0;
	float currentAlphaHits2 = 0;
	float currentBetaHits = 0;
	float currentBetaHits2 = 0;
	float currentX = -0.098;

	float youdenStatistic;
	float generalStatistic;

	for(int k = 0; k<100; k++)
	{
		currentAlphaHits1 = float(ratio)*(analysisAlphaHistogram->Integral(1, 100, k, 100));
		currentAlphaHits2 = float(ratio)*(analysisAlphaHistogram->Integral(1, 100, 1, k));
		currentBetaHits = analysisBetaHistogram->Integral(1, 100, 1, k);
		currentBetaHits2 = analysisBetaHistogram->Integral(1, 100, k, 100);

		if(not(currentBetaHits==0 & currentAlphaHits2==0))
		{
			youdenStatistic = (currentBetaHits/(currentBetaHits+currentBetaHits2)) + (currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2))-float(1);
			generalStatistic = (currentBetaHits)/sqrt(currentBetaHits+currentAlphaHits2);
		}
		else 
		{
			youdenStatistic = 0;
			generalStatistic = 0;
		}		

		//fill for cut selection stats
		youdenSelection->Fill(currentX, youdenStatistic);
		generalSelection->Fill(currentX, generalStatistic);

		currentX+= .002;
	}
	
	float youdenClassifierBin = youdenSelection->GetMaximumBin();
	float youdenClassifierMax = youdenSelection->GetXaxis()->GetBinCenter(youdenClassifierBin);
	float generalClassifierBin = generalSelection->GetMaximumBin();
	float generalClassifierMax = generalSelection->GetXaxis()->GetBinCenter(generalClassifierBin);
	float youdenNhitMax = youdenSelection->GetMaximum();
	float generalNhitMax = generalSelection->GetMaximum();
	
	float allAlphas = analysisAlphaHistogram->Integral(1,100,1,100);
	float allBetas = analysisBetaHistogram->Integral(1,100,1,100);	
	
	if(allAlphas==0) allAlphas = 0.000001;
	if(allBetas==0) allBetas = 0.000001;

	float youdenAlphaRejection = analysisAlphaHistogram->Integral(1, 100, youdenClassifierBin, 100) / allAlphas;
	float youdenBetaAcceptance = analysisBetaHistogram->Integral(1, 100, 1, youdenClassifierBin) / allBetas;
	float generalAlphaRejection = analysisAlphaHistogram->Integral(1, 100, generalClassifierBin, 100) / allAlphas;
	float generalBetaAcceptance = analysisBetaHistogram->Integral(1, 100, 1, generalClassifierBin) /allBetas;

	static float values[9];
	values[0] = youdenClassifierMax;
	values[1] = youdenNhitMax;
	values[2] = generalClassifierMax;
	values[3] = generalNhitMax;
	values[4] = youdenAlphaRejection;
	values[5] = youdenBetaAcceptance;
	values[6] = generalAlphaRejection;
	values[7] = generalBetaAcceptance;
	values[8] = meanNhit;

	return values;
}