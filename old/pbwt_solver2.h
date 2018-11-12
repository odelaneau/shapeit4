#ifndef _PBWT_SOLVER2_H
#define _PBWT_SOLVER2_H

#include <utils/otools.h>
#include <containers/haplotype_set.h>

class pbwt_solver2 {
private:
	bitmatrix & H;
	unsigned int n_site, n_main_hap, n_ref_hap, n_total_hap, n_total, n_resolved;
	vector < vector < int > > pbwt_clusters, pbwt_indexes, pbwt_divergences;
	vector < char > Guess;
	vector < bool > Het, Mis;
	vector < float > scoreBit;
	vector < vector < bool > > Amb;

public:
	pbwt_solver2(haplotype_set &);
	~pbwt_solver2();
	void free();

	void forward(genotype_set &);
	void backward(genotype_set &);
};

#endif
