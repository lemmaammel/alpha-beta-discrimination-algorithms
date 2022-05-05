// contains functions to generate scatterplots of real alpha and beta events, find the optimal cutoff values, and calculate the data relevant to those values

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <string>
#include <cmath>
#include <iostream>

bool getEventCoordinates(double posx, double posy, double posz, double rho, double z) {
    bool inRange = true;
	if (posz < z-500 || posz > z+500) {
        inRange = false;
    }
    double posrho = sqrt((posx*posx) + (posy*posy));
	if (posrho*posrho/36000000 < rho-0.056 || posrho*posrho/36000000 > rho+0.056) {
        inRange = false;
    }
	return inRange;
}

TH2D* nHitHistogramRealData(std::string filename, double rho, double z, int style) {

	TCanvas *c1 = new TCanvas("c1", "Classification Histogram", 200, 10, 1300, 1300);
	c1->cd();

	TH2D* histogram = new TH2D("Events", "N_{hit} Vs. Classification", 100, 100, 1200, 100, -0.1, 0.1);
	histogram->SetMarkerColor(style);
	histogram->SetMarkerStyle(20);

	TFile *f = TFile::Open(filename.c_str());
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

		if(getEventCoordinates(posx, posy,posz, rho, z)) {
			histogram->Fill(nhits, berkeleyAlphaBeta/nhits);
		}
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

void nHitHistogramRealDataNew(std::string filename, std::string filename2, double rho, double z, int style1, int style2) {
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


double* rejectionInfo(char* filename, char* filename2, double rho, double z, double ratio) {
	TH2D* analysisAlphaHistogram = nHitHistogramRealData(filename, rho, z, 1);
	TH2D* analysisBetaHistogram = nHitHistogramRealData(filename2, rho, z, 1);

	//cut selection histograms
	TH1D* youdenSelection = new TH1D("Youden's J Statistic", "Youden's J Statistic", 100, -0.1, 0.1);
	TH1D* generalSelection = new TH1D("General Cut Statistic", "General Cut Statistic", 100, -0.1, 0.1);

	double currentAlphaHits1 = 0;
	double currentAlphaHits2 = 0;
	double currentBetaHits = 0;
	double currentBetaHits2 = 0;
    //Where does this value come from?
	double currentX = -0.098;

	double youdenStatistic;
	double generalStatistic;

    //We should think about how to have consistent bin sizes instead of having to hard code 100 in places
    for (int k = 0; k < 100; k++) {
        //We should pick more descriptive names than 1/2, maybe Acc/Rej
		currentAlphaHits1 = ratio*analysisAlphaHistogram->Integral(1, 100, k, 100);
		currentAlphaHits2 = ratio*analysisAlphaHistogram->Integral(1, 100, 1, k);
		currentBetaHits = analysisBetaHistogram->Integral(1, 100, 1, k);
		currentBetaHits2 = analysisBetaHistogram->Integral(1, 100, k, 100);

		if (!(currentBetaHits == 0 && currentAlphaHits2 == 0)) {
            //The denominators for Youden should be constants, so we can likely avoid some computations
			youdenStatistic = (currentBetaHits/(currentBetaHits+currentBetaHits2)) + (currentAlphaHits1/(currentAlphaHits1+currentAlphaHits2)) - 1.0;
			generalStatistic = (currentBetaHits)/sqrt(currentBetaHits+currentAlphaHits2);
		}
		else  {
			youdenStatistic = 0;
			generalStatistic = 0;
		}

		//fill for cut selection stats
		youdenSelection->Fill(currentX, youdenStatistic);
		generalSelection->Fill(currentX, generalStatistic);

        //Where does this value come from?
		currentX+= .002;
	}

	double youdenClassifierBin = youdenSelection->GetMaximumBin();
	double youdenClassifierMax = youdenSelection->GetXaxis()->GetBinCenter(youdenClassifierBin);
	double generalClassifierBin = generalSelection->GetMaximumBin();
	double generalClassifierMax = generalSelection->GetXaxis()->GetBinCenter(generalClassifierBin);
	double youdenNhitMax = youdenSelection->GetMaximum();
	double generalNhitMax = generalSelection->GetMaximum();

	double allAlphas = analysisAlphaHistogram->Integral(1,100,1,100);
	double allBetas = analysisBetaHistogram->Integral(1,100,1,100);

	if(allAlphas==0) allAlphas = 0.000001;
	if(allBetas==0) allBetas = 0.000001;

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
	values.push_back(allAlphas);
	values.push_back(allBetas);

	return values;
}

std::vector<double> cutPerformanceValues(std::string filename, std::string filename2, double rho, double z, double ratio, double cutValue1) {
	//Where do these numbers come from
    double cutValue = abs(cutValue1+0.1)/0.002;
	if(cutValue > 100) cutValue = 100;
	if(cutValue < 0) cutValue = 0;

	TH2D* analysisAlphaHistogram = nHitHistogramRealData(filename, rho, z, ratio);
	TH2D* analysisBetaHistogram = nHitHistogramRealData(filename2, rho, z, ratio);

	double alphaHits = double(ratio)*(analysisAlphaHistogram->Integral(1, 100, cutValue, 100));
	double alphaHits2 = double(ratio)*(analysisAlphaHistogram->Integral(1, 100, 1, cutValue));
	double betaHits = analysisBetaHistogram->Integral(1, 100, 1, cutValue);

	double betaHits2 = analysisBetaHistogram->Integral(1, 100, cutValue, 100);

	double youdenCutValue = (betaHits/(betaHits+betaHits2)) + (alphaHits/(alphaHits+alphaHits2)) - double(1);
	double generalCutValue = (betaHits)/sqrt(betaHits+alphaHits2);

	double allAlphas = analysisAlphaHistogram->Integral(1,100,1,100);
	double allBetas = analysisBetaHistogram->Integral(1,100,1,100);

	if(allAlphas==0) allAlphas = 0.000001;
	if(allBetas==0) allBetas = 0.000001;

	double alphaRejection = analysisAlphaHistogram->Integral(1, 100, cutValue, 100) / allAlphas;
	double betaAcceptance = analysisBetaHistogram->Integral(1, 100, 1, cutValue) / allBetas;

    std::vector<double> values;
    values.push_back(youdenCutValue);
    values.push_back(generalCutValue);
    values.push_back(alphaRejection)
    values.push_back(betaRejection);

	return values;
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


