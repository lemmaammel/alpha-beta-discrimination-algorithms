// returns average z and rho coordinates from real data files
#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <string>

std::vector<double> getRealEventCoordinates(const std::string& filename, const std::string& fitname) {
    TFile *f = TFile::Open(filename.c_str());
    TTree *t = (TTree*)f->Get("output");

    double posx, posy, posz;
    bool fitValid;
    t->SetBranchAddress("posx", &posx);
    t->SetBranchAddress("posy", &posy);
    t->SetBranchAddress("posz", &posz);
    t->SetBranchAddress("fitValid", &fitValid);

    double averageXCoordinate = 0;
    double averageYCoordinate = 0;
    double averageZCoordinate = 0;
    double averageRhoCoordinate = 0;
    size_t counter = 0;

    for (size_t i = 0; i < t->GetEntries(); i++) {
        t->GetEntry(i);

        if (fitValid) {
            continue;
        }

        averageXCoordinate += posx;
        averageYCoordinate += posy;
        averageZCoordinate += posz;
        averageRhoCoordinate += sqrt(posx*posx + posy*posy);
        counter++;
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

std::vector<double> getSimulatedEventCoordinates(const std::string& filename, const std::string& fitname) {
    RAT::DU::DSReader dsReader(filename);

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

            averageZCoordinate += pos.Z();
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


