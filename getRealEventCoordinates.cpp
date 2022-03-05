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
#include <string>
#include <iostream>

int* getEventCoordinates(const std::string& filename, const std::string& fitname, float rho, float z)
{
	const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
	const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();
	
	TFile *f = TFile::Open(filename);
	TTree *t = (TTree*)f->Get("output");

	double posx, posy, posz;
	bool fitValid;
	t->SetBranchAddress("posx", &posx);
	t->SetBranchAddress("posy", &posy);
	t->SetBranchAddress("posz", &posz);
	t->SetBranchAddress("fitValid", &fitValid);

	float averageXCoordinate = 0;
	float averageYCoordinate = 0;
	float averageZCoordinate = 0;
	float averageRhoCoordinate = 0;
	float counter = 0;
	
	for(size_t i=0; i<t->GetEntries(); i++)
	{	
		t->GetEntry(i);

		if(fitValid) continue;

		averageXCoordinate = averageXCoordinate + posX;
		averageYCoordinate = averageYCoordinate + posY;
		averageZCoordinate = averageZCoordinate + posZ;
		counter++;
	}

	averageXCoordinate = averageXCoordinate / counter;
	averageYCoordinate = averageYCoordinate / counter;
	averageZCoordinate = averageZCoordinate / counter;
	averageRhoCoordinate = sqrt((averageXCoordinate*averageXCoordinate)+(averageYCoordinate*averageYCoordinate));

	static int values[2];
	values[0] = (int)averageZCoordinate;
	values[1] = (int)averageRhoCoordinate;
	return values;
}
	
