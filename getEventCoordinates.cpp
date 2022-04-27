#include <RAT/DU/DSReader.hh>
#include <RAT/DS/Entry.hh>
#include <RAT/DS/EV.hh>
#include <RAT/DS/FitResult.hh>
#include <cmath>
#include <string>

std::vector<double> getEventCoordinates(const std::string& filename, const std::string& fitname) {
	RAT::DU::DSReader dsReader(filename);

	double averageXCoordinate = 0.0;
    double averageYCoordinate = 0.0;
    double averageZCoordinate = 0.0;
    double averageRhoCoordinate = 0.0;
    size_t counter = 0;

	for (size_t i = 0; i < dsReader.GetEntryCount(); i++) {
		const RAT::DS::Entry& rDS = dsReader.GetEntry(i);

		for (size_t iEV = 0; iEV < rDS.GetEVCount(); iEV++) {
			const RAT::DS::EV& rEV = rDS.GetEV(iEV);

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

			RAT::DS::FitResult fResult = rEV.GetFitResult(fitname);
			RAT::DS::FitVertex fVertex = fResult.GetVertex(0);
            TVector3 pos = fVertex.GetPosition();

			averageXCoordinate += pos.X();
			averageYCoordinate += pos.Y();
			averageZCoordinate += pos.Z();
            //We want the average rho coordinate, not the rho coordinate corresponding
            //to the average x and average y
			averageRhoCoordinate += pos.Perp();
            counter++;
		}
	}

	averageXCoordinate /= counter;
	averageYCoordinate /= counter;
	averageZCoordinate /= counter;
	averageRhoCoordinate /= counter;

    std::vector<double> values;
	values.push_back(averageZCoordinate);
	values.push_back(averageRhoCoordinate);

	return values;
}

