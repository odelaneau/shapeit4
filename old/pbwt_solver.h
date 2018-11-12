#ifndef _PBWT_SOLVER_H
#define _PBWT_SOLVER_H

#include <utils/otools.h>
#include <containers/haplotype_set.h>

class pbwt_solver {
private:
	bitmatrix & H;
	unsigned int n_site, n_main_hap, n_ref_hap, n_total_hap, n_total, n_resolved;
	double MP[4], MUT[3];
	vector < vector < int > > pbwt_clusters, pbwt_indexes;
	vector < vector < bool > > vecA;
	vector < int > B, GuessPrev, GuessNext;
	vector < bool > Guess, Het, Mis;
	vector < float > Prob;

public:
	pbwt_solver(haplotype_set &);
	~pbwt_solver();
	void free();

	void forward(genotype_set &, int, int);
	void backward(genotype_set &, int, int);
};

#endif
