#include <RAT/DU/DSReader.hh>
#include <RAT/DU/Utility.hh>
#include <RAT/DB.hh>
#include <RAT/DU/PMTInfo.hh>
#include <RAT/DU/LightPathCalculator.hh>
#include <RAT/DU/GroupVelocity.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/PMT.hh>
#include <RAT/DS/FitResult.hh>
#include <RAT/DS/MC.hh>
#include <RAT/DS/MCParticle.hh>
#include <cmath>
#include <TStyle.h>
#include <iostream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <RAT/DS/Meta.hh>

#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string> 
using namespace std;

bool getEventCoordinates(float posx, float posy, float posz, float rho, float z)
{
	bool inRange = true;
	if(posz<z-500 || posz>z+500) inRange = false;
	float posrho = sqrt((posx*posx)+(posy*posy));
	if(posrho*posrho/36000000<rho-0.056 || posrho*posrho/36000000>rho+0.056) inRange = false;
	cout << inRange;
	return inRange;
}

TH2D* nHitHistogramRealData(char* filename, float rho, float z, int style)
{

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* histogram = new TH2D("Events", "N_{hit} Vs. Classification", 100, 100, 1200, 100, -0.1, 0.1);
	histogram->SetMarkerColor(style);
	histogram->SetMarkerStyle(20);

	const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
	const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();
	
	TFile *f = TFile::Open(filename);
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
	
	for(size_t i=0; i<k; i++)
	{	
		t->GetEntry(i);
	
		if(nhits<0) continue;

		/*if(getEventCoordinates(posx, posy,posz, rho, z))*/ histogram->Fill(nhits, berkeleyAlphaBeta/nhits);
	}

	f->Close();

	histogram->GetYaxis()->SetTitle("Classification Value/Number of Hits");
	histogram->GetYaxis()->SetTitleOffset(1.2);
	histogram->GetXaxis()->SetTitle("Number of Hits");
	histogram->GetXaxis()->SetLabelSize(.03);
	histogram->GetYaxis()->SetLabelSize(.02);
	
	histogram->Draw();
	
	c1->Print("realDataHistogramClassificationNhit2.pdf", "pdf");

	return histogram;
}

void nHitHistogramRealDataNew(char* filename, char* filename2, float rho, float z, int style1, int style2)
{
	TH2D* analysisAlphaHistogram = nHitHistogramRealData(filename, rho, z, style1);
	TH2D* analysisBetaHistogram = nHitHistogramRealData(filename2, rho, z, style2);
	
	TCanvas *c2 = new TCanvas("c2", "Classification Histogram", 200, 10, 1300, 1300);
	c2->cd();
	
	analysisAlphaHistogram->Draw();
	analysisBetaHistogram->Draw("SAME");
	
	// build a legend
	TLegend *legend = new TLegend(0.1,0.7,0.48,.9);
	legend->SetHeader("Legend");
	legend->AddEntry(analysisAlphaHistogram, "#alpha Events", "p");
	legend->AddEntry(analysisBetaHistogram, "#beta Events", "p");
	legend->Draw();
	
	c2->Print("realDataHistogramClassificationNhit4.pdf", "pdf");
}


float* rejectionInfo(char* filename, char* filename2, float rho, float z, float ratio)
{

	TH2D* analysisAlphaHistogram = nHitHistogramRealData(filename, rho, z,1);
	TH2D* analysisBetaHistogram = nHitHistogramRealData(filename2, rho, z,1);

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

	static float values[10];
	values[0] = youdenClassifierMax;
	values[1] = youdenNhitMax;
	values[2] = generalClassifierMax;
	values[3] = generalNhitMax;
	values[4] = youdenAlphaRejection;
	values[5] = youdenBetaAcceptance;
	values[6] = generalAlphaRejection;
	values[7] = generalBetaAcceptance;
	values[8] = allAlphas;
	values[9] = allBetas;
	
	return values;
}

float* cutPerformanceValues(char* filename, char* filename2, float rho, float z, float ratio, float cutValue1)
{
	
	float cutValue = abs(cutValue1+0.1)/0.002;
	if(cutValue > 100) cutValue = 100;
	if(cutValue < 0) cutValue = 0;
	
	TH2D* analysisAlphaHistogram = nHitHistogramRealData(filename, rho, z,ratio);
	TH2D* analysisBetaHistogram = nHitHistogramRealData(filename2, rho, z,ratio);

	float alphaHits = float(ratio)*(analysisAlphaHistogram->Integral(1, 100, cutValue, 100));
	float alphaHits2 = float(ratio)*(analysisAlphaHistogram->Integral(1, 100, 1, cutValue));
	float betaHits = analysisBetaHistogram->Integral(1, 100, 1, cutValue);
	
	float betaHits2 = analysisBetaHistogram->Integral(1, 100, cutValue, 100);
	
	float youdenCutValue = (betaHits/(betaHits+betaHits2)) + (alphaHits/(alphaHits+alphaHits2)) - float(1);
	float generalCutValue = (betaHits)/sqrt(betaHits+alphaHits2);
	
	float allAlphas = analysisAlphaHistogram->Integral(1,100,1,100);
	float allBetas = analysisBetaHistogram->Integral(1,100,1,100);	
	
	if(allAlphas==0) allAlphas = 0.000001;
	if(allBetas==0) allBetas = 0.000001;
	
	float alphaRejection = analysisAlphaHistogram->Integral(1, 100, cutValue, 100) / allAlphas;
	float betaAcceptance = analysisBetaHistogram->Integral(1, 100, 1, cutValue) / allBetas;

	static float values [] = { youdenCutValue, generalCutValue, alphaRejection, betaRejection };
	
	return values;
}

float* averageValues(char* filename)
{
	const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
	const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();
	
	TFile *f = TFile::Open(filename);
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

	for(size_t i=0; i<k; i++)
	{	
		t->GetEntry(i);
	
		if(nhits<0) continue;
		if(sqrt(pow(posX,2) + pow(posY,2)) < 0.0 || sqrt(pow(posX,2) + pow(posY,2)) > 5.0) continue;
		if(posZ < -5.0 || posZ > 5.0) continue;

		avgX = avgX + posx;
		avgY = avgY + posy;
		avgZ = avgZ + posz;	
	}

	f->Close();

	static float values [] = { avgX/k, avgY/k, avgZ/k, sqrt(pow(values[0],2) + pow(values[1],2)) };
	
	return values;
}


