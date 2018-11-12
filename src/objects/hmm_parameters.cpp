//$Id: hmm_parameters.cpp 613 2012-07-26 15:06:09Z koskos $

#include <objects/hmm_parameters.h>

hmm_parameters::hmm_parameters() {
	ed = 0.0001;
	ee = 0.9999;
}

hmm_parameters::~hmm_parameters() {
	t.clear();
	nt.clear();
	tfreq.clear();
}

void hmm_parameters::initialise(variant_map & mapG, int Neff, int Nhap, int Nstate, bool gmap) {
	t = vector < double > (mapG.size() - 1, 0.0);
	nt = vector < double > (mapG.size() - 1, 0.0);
	tfreq = vector < double > (mapG.size() - 1, 0.0);
	for (int l = 1 ; l < mapG.size() ; l ++) {
		double rho;
		if (gmap) rho = 0.04 * Neff * (mapG.vec_pos[l]->cm - mapG.vec_pos[l-1]->cm);
		else 0.0004 *  (mapG.vec_pos[l]->bp - mapG.vec_pos[l-1]->bp);
		if (rho == 0.0) rho = 0.00001;
		nt[l-1] = exp(-1.0 * rho / Nhap);
		t[l-1] = 1-nt[l-1];
		tfreq[l-1] = t[l-1] * 1.0 / Nstate;
	}
}

