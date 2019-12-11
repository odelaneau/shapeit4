/*******************************************************************************
 * Copyright (C) 2018 Olivier Delaneau, University of Lausanne
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/
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
	int start_missing;
	int start_transition;
	int stop_locus;
	int stop_segment;
	int stop_ambiguous;
	int stop_missing;
	int stop_transition;

	coordinates() {
		start_locus = 0;
		start_segment = 0;
		start_ambiguous = 0;
		start_missing = 0;
		start_transition = 0;
		stop_locus = 0;
		stop_segment = 0;
		stop_ambiguous = 0;
		stop_missing = 0;
		stop_transition = 0;
	}

	~coordinates() {
		start_locus = 0;
		start_segment = 0;
		start_ambiguous = 0;
		start_missing = 0;
		start_transition = 0;
		stop_locus = 0;
		stop_segment = 0;
		stop_ambiguous = 0;
		stop_missing = 0;
		stop_transition = 0;
	}

	string toString() {
		string str="";
		str += "L=[" + stb.str(start_locus) + "->" + stb.str(stop_locus) + "]";
		return str;
	}
};

class compute_job {
public:
	variant_map & V;
	genotype_set & G;
	haplotype_set & H;
	vector < double > T;
	vector < float > M;
	vector < coordinates > C;
	vector < vector < unsigned int > > Kvec;

	compute_job(variant_map & , genotype_set & , haplotype_set & , unsigned int n_max_transitions , unsigned int n_max_missing);
	~compute_job();

	void free();
	void reset();
	void make(unsigned int, double);
	unsigned int size();
	void maskingTransitions(unsigned int, double);
	bool reccursive_window_splitting(double, int, int, vector < int > &, vector < int > &, vector < double > &, vector < double > &, vector < int > &);
};

inline
unsigned int compute_job::size() {
	 return C.size();
}

#endif
