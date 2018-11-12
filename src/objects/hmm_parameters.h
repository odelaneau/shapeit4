//$Id: hmm_parameters.h 597 2012-07-18 14:00:06Z koskos $

#ifndef _HMM_PARAMETERS_H
#define _HMM_PARAMETERS_H

#include <utils/otools.h>
#include <containers/variant_map.h>

class hmm_parameters {
public :
	//DATA
	vector < double > t;
	vector < double > tfreq;
	vector < double > nt;
	double ee;
	double ed;
	double efreq;
	double dfreq;

	//CONSTRUCTOR/DESTRUCTOR
	hmm_parameters();
	~hmm_parameters();

	//METHODS
	void initialise(variant_map &, int, int, int, bool);
};

#endif
