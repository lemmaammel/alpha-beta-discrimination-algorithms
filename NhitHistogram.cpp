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
//#include <format.h>

#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string> 
using namespace std;

TCanvas* NhitHistogram(const string& filename, const string& filename2, const string& fitName, const string& className, const string& classification)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* histogramClassificationNhit1 = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -.1, 0.1);
	TH2D* histogramClassificationNhit2 = new TH2D("#alpha Events","",100,100,1000,100,-0.1,0.1);
	TH2D* histograms [2] = {histogramClassificationNhit1, histogramClassificationNhit2};

	histogramClassificationNhit1->SetMarkerColor(4);
	histogramClassificationNhit2->SetMarkerColor(6);
	histograms[0]->SetMarkerStyle(20);
	histograms[1]->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);	

	//const string& files[2] = {filename, filename2};

	// for(every entry in filename)
	//    for(every triggered event in entry
	//     for(i=0; i<nhit)
	//      histogramClassificiationNhit->fill(classificationValue)

	//loop through all entries in filename
	for(int m = 0; m<2; m++)
	{
		string files[2] = {filename, filename2};
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
					histograms[m]->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
					//histograms[m]->Fill(200, x/100000);
				}				
			}
		}
	}

	histogramClassificationNhit1->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	histogramClassificationNhit1->GetYaxis()->SetTitleOffset(1.2);
	histogramClassificationNhit1->GetXaxis()->SetTitle("Number of Hits");
	histogramClassificationNhit1->GetXaxis()->SetLabelSize(.03);
	histogramClassificationNhit1->GetYaxis()->SetLabelSize(.02);
	
	histograms[0]->Draw();
	histograms[1]->Draw("same");
	
	// build a legend
	TLegend *legend = new TLegend(0.1,0.7,0.48,.9);
	legend->SetHeader("Legend");
	legend->AddEntry(histograms[0], "#alpha Events", "p");
	legend->AddEntry(histograms[1], "#beta Events", "p");
	legend->Draw();


	c1->Print("partialFillHistogramClassificationNhit.pdf", "pdf");
	return c1;
}
	
TH2D* NhitHistogramAlpha(const string& filename, const string& filename2, const string& fitName, const string& className, const string& classification)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* histogramClassificationNhit1 = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -.1, .1);
	TH2D* histogramClassificationNhit2 = new TH2D("#alpha Events","",100,100,1000,100,-0.1,.1);
	TH2D* histograms [2] = {histogramClassificationNhit1, histogramClassificationNhit2};

	histogramClassificationNhit1->SetMarkerColor(4);
	histogramClassificationNhit2->SetMarkerColor(6);
	histograms[0]->SetMarkerStyle(20);
	histograms[1]->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);	

	//const string& files[2] = {filename, filename2};

	// for(every entry in filename)
	//    for(every triggered event in entry
	//     for(i=0; i<nhit)
	//      histogramClassificiationNhit->fill(classificationValue)

	//loop through all entries in filename
	for(int m = 0; m<2; m++)
	{
		string files[2] = {filename, filename2};
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
					histograms[m]->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
					//histograms[m]->Fill(200, x/100000);
				}				
			}
		}
	}

	histogramClassificationNhit1->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	histogramClassificationNhit1->GetYaxis()->SetTitleOffset(1.2);
	histogramClassificationNhit1->GetXaxis()->SetTitle("Number of Hits");
	histogramClassificationNhit1->GetXaxis()->SetLabelSize(.03);
	histogramClassificationNhit1->GetYaxis()->SetLabelSize(.02);
	
	histograms[0]->Draw();
	histograms[1]->Draw("same");
	
	// build a legend
	TLegend *legend = new TLegend(0.1,0.7,0.48,.9);
	legend->SetHeader("Legend");
	legend->AddEntry(histograms[0], "#alpha Events", "p");
	legend->AddEntry(histograms[1], "#beta Events", "p");
	legend->Draw();


	//c1->Print("histogramClassificationNhit.pdf", "pdf");
	return histogramClassificationNhit2;
}	



TH2D* NhitHistogramBeta(const string& filename, const string& filename2, const string& fitName, const string& className, const string& classification)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* histogramClassificationNhit1 = new TH2D("#beta Events", "N_{hit} Vs. Classification", 100, 100, 1000, 100, -.1, .1);
	TH2D* histogramClassificationNhit2 = new TH2D("#alpha Events","",100,100,1000,100,-0.1,.1);
	TH2D* histograms [2] = {histogramClassificationNhit1, histogramClassificationNhit2};

	histogramClassificationNhit1->SetMarkerColor(4);
	histogramClassificationNhit2->SetMarkerColor(6);
	histograms[0]->SetMarkerStyle(20);
	histograms[1]->SetMarkerStyle(20);

	RAT::DB::Get()->SetAirplaneModeStatus(true);	

	//loop through all entries in filename
	for(int m = 0; m<2; m++)
	{
		string files[2] = {filename, filename2};
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
					histograms[m]->Fill(numberHits.GetCalPMTs().GetAllCount(), cValue/(numberHits.GetCalPMTs().GetAllCount()));
					//histograms[m]->Fill(200, x/100000);
				}				
			}
		}
	}

	histogramClassificationNhit1->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	histogramClassificationNhit1->GetYaxis()->SetTitleOffset(1.2);
	histogramClassificationNhit1->GetXaxis()->SetTitle("Number of Hits");
	histogramClassificationNhit1->GetXaxis()->SetLabelSize(.03);
	histogramClassificationNhit1->GetYaxis()->SetLabelSize(.02);
	
	histograms[0]->Draw();
	histograms[1]->Draw("same");
	
	// build a legend
	TLegend *legend = new TLegend(0.2,0.2,0.8,.8);
	legend->SetHeader("Legend");
	legend->AddEntry(histograms[0], "#alpha Events", "p");
	legend->AddEntry(histograms[1], "#beta Events", "p");
	legend->Draw();


	//c1->Print("histogramClassificationNhit.pdf", "pdf");
	return histogramClassificationNhit1;
}

TCanvas* AlphaRejectionHistogram(const string& filename, const string& filename2, const string& fitname, const string& classname, const string& classification, float ratio)
{

	TH2D* analysisAlphaHistogram = NhitHistogramAlpha(filename, filename2, fitname, classname, classification);

	TH2D* analysisBetaHistogram = NhitHistogramBeta(filename,filename2, fitname, classname, classification);

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

	//	alphaRejectionHistogram->Fill(currentX, currentAlphaHits1/totalAlphaHits);
	//	betaAcceptanceHistogram->Fill(currentX, currentBetaHits/totalBetaHits);
		alphaRejectionHistogram->SetBinContent(k, currentAlphaHits1/totalAlphaHits);
		betaAcceptanceHistogram->SetBinContent(k, currentBetaHits/totalBetaHits);

		if(not(currentBetaHits==0 & currentAlphaHits2==0))
		{
			//betaSampleFraction->Fill(currentX, currentBetaHits/(currentBetaHits+currentAlphaHits2));	
			betaSampleFraction->SetBinContent(k, currentBetaHits/(currentBetaHits+currentAlphaHits2));
			youdenStatistic = (currentBetaHits/(currentBetaHits+currentBetaHits2)) + (currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2))-1;
			generalStatistic = (currentBetaHits)/(sqrt(currentBetaHits+currentAlphaHits2)*50);
		}
		else 
		{
			//betaSampleFraction->Fill(currentX, 1.0);
			betaSampleFraction->SetBinContent(k, 1.0);
			youdenStatistic = 0;
			generalStatistic = 0;
		}		

		//fill for cut selection stats
		//youdenSelection->Fill(currentX, youdenStatistic);
		//generalSelection->Fill(currentX, generalStatistic);
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


float* AlphaRejectionInfo(const string& filename, const string& filename2, const string& fitname, const string& classname, const string& classification, float ratio)
{

	TH2D* analysisAlphaHistogram = NhitHistogramAlpha(filename, filename2, fitname, classname, classification);

	TH2D* analysisBetaHistogram = NhitHistogramBeta(filename,filename2, fitname, classname, classification);

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
			cout<"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
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


/*
void compileHistogram()
{

	int copies = 3;
	string loopvars[9] = ["0_0_3", "0_0_4", "1_0_2","1_0_4","2_0_2","2_0_3","3_0_2","3_0_3","4_0_2"];

	for(int x=0; x<9; x++)
	{
    
    		TCanvas* c1 = new TCanvas("c1", "", 200, 10, 1300, 1300);
		c1->cd();

    		string name1 =  format("/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_e-*", loopvars[x]);
    		string name2 = format("/data/snoplus/home/ammel/projects/alphabeta_test/{}_dec2020_recoord_alpha*", loopvars[x]);

    		TH2D* htot = NhitHistogram(format("{}.root", name2), format("{}.root", name1), "partialFitter", "BerkeleyAlphaBeta:partialFitter", "likelihood");

    		htot->Draw("AL");

    		c1->Print(format("compiledHistogram{}.pdf", loopvars[x]);
	}
}
*/
