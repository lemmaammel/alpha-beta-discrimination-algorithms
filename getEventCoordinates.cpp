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

int* getEventCoordinates(const std::string& filename, const std::string& fitname)
{
	const RAT::DU::ReconCorrector &eCorr = RAT::DU::Utility::Get()->GetReconCorrector();
	const RAT::DU::ReconCalibrator &eCalib = RAT::DU::Utility::Get()->GetReconCalibrator();
	
	RAT::DU::DSReader dsReader(filename);

	float averageXCoordinate = 0;
	float averageYCoordinate = 0;
	float averageZCoordinate = 0;
	float averageRhoCoordinate = 0;
	float counter = 0;
	
	for(size_t i=0; i<dsReader.GetEntryCount(); i++)
	{	
	
		const RAT::DS::Entry& rDS = dsReader.GetEntry(i);

		for(size_t iEV=0; iEV<rDS.GetEVCount(); iEV++)
		{
			const RAT::DS::EV& rEV = rDS.GetEV(iEV);

			if(!rEV.FitResultExists(fitname)) continue;
			if(!rEV.GetFitResult(fitname).GetVertex(0).ContainsEnergy()) continue;
			if(!rEV.GetFitResult(fitname).GetVertex(0).ValidEnergy()) continue;

			cout<<"test3";
			RAT::DS::FitResult fResult = rEV.GetFitResult(fitname);
			RAT::DS::FitVertex fVertex = fResult.GetVertex(0);

			averageXCoordinate = averageXCoordinate + fVertex.GetPosition().X();
			averageYCoordinate = averageYCoordinate + fVertex.GetPosition().Y();
			averageZCoordinate = averageZCoordinate + fVertex.GetPosition().Z();
			counter++;
		}
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
	
