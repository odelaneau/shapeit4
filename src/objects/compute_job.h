#ifndef _COMPUTE_THREAD_H
#define _COMPUTE_THREAD_H

#include <utils/otools.h>

#include <containers/haplotype_set.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>

class coordinates {
public:
	int start_locus;
	int start_segment;
	int start_ambiguous;
	int start_transition;
	int stop_locus;
	int stop_segment;
	int stop_ambiguous;
	int stop_transition;

	coordinates() {
		start_locus = 0;
		start_segment = 0;
		start_ambiguous = 0;
		start_transition = 0;
		stop_locus = 0;
		stop_segment = 0;
		stop_ambiguous = 0;
		stop_transition = 0;
	}

	~coordinates() {
		start_locus = 0;
		start_segment = 0;
		start_ambiguous = 0;
		start_transition = 0;
		stop_locus = 0;
		stop_segment = 0;
		stop_ambiguous = 0;
		stop_transition = 0;
	}
};

class compute_job {
public:
	variant_map & V;
	genotype_set & G;
	haplotype_set & H;
	vector < double > T;
	vector < coordinates > C;
	vector < vector < unsigned int > > Kvec;

	compute_job(variant_map & , genotype_set & , haplotype_set & , unsigned int n_max_transitions);
	~compute_job();

	void free();
	void reset();
	void make(unsigned int, double);
	unsigned int size();
	void maskingTransitions(unsigned int, double);
};

inline
unsigned int compute_job::size() {
	 return C.size();
}

#endif
