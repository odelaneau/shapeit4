#ifndef _GENOTYPE_SET_H
#define _GENOTYPE_SET_H

#include <utils/otools.h>

#include <objects/genotype/genotype_header.h>
#include <containers/variant_map.h>

class genotype_set {
public:
	//DATA
	int n_site, n_ind;
	vector < genotype * > vecG;

	//CONSTRUCTOR/DESTRUCTOR
	genotype_set();
	~genotype_set();

	//METHODS
	void imputeMonomorphic(variant_map &);
	unsigned int largestNumberOfTransitions();
	unsigned long numberOfSegments();
	void masking();
	void solve();
};

#endif
