#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <string>

std::vector<double> getEventCoordinates(const std::string& filename, const std::string& fitname) {
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

		if(fitValid) continue;

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

