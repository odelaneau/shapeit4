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
#ifndef _HAPLOTYPE_SET_H
#define _HAPLOTYPE_SET_H

#include <utils/otools.h>

#include <containers/bitmatrix.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>

struct IBD2track {
	int ind, from, to;

	IBD2track(int _ind, int _from, int _to) {
		ind = _ind;
		from = _from;
		to = _to;
	}

	bool operator<(const IBD2track & rhs) const {
		if (ind < rhs.ind) return true;
		if (ind > rhs.ind) return false;
		return (from < rhs.from);
	}

	bool overlap (const IBD2track & rhs) const {
		return ((ind==rhs.ind) && (rhs.to >= from) && (rhs.from <= to));
	}

	bool merge (const IBD2track & rhs) {
		from = min(from, rhs.from);
		to = max(to, rhs.to);
	}


} ;

class haplotype_set {
public:
	//Haplotype Data
	bitmatrix H_opt_hap;		// Bit matrix of haplotypes (haplotype first). Transposed version of H_opt_var.
	bitmatrix H_opt_var;		// Bit matrix of haplotypes (variant first). Transposed version of H_opt_hap.
	unsigned long n_site;		// #variants
	unsigned long n_hap;		// #haplotypes
	unsigned long n_ind;		// #individuals

	//PBWT parameters
	double pbwt_modulo;				// Modulo used to store PBWT indexes (--pbwt-modulo)
	unsigned long pbwt_depth;		// #neighbours in the PBWT to use for conditioning (--pbwt-depth)
	unsigned long pbwt_mac;			// Minor Allele Count to consider in PBWT pass (--pbwt-mac)
	double pbwt_mdr;				// Missinga Data Rate to consider in PBWT pass (--pbwt-mdr)
	unsigned int nthreads;			// Number of threads (--thread)

	//PBWT indexes & arrays
	unsigned long pbwt_nstored;		//#variants with PBWT indexes stored
	vector < double > pbwt_cm;		//Variants at which PBWT is evaluated
	vector < int > pbwt_grp;		//Variant groups based on cm positions
	vector < int > pbwt_evaluated;	//Variants at which PBWT is evaluated
	vector < int > pbwt_stored;		//Variants at which PBWT is stored
	vector < int > pbwt_parray;		//PBWT prefix array
	vector < int > pbwt_darray;		//PBWT divergence array
	vector < int > pbwt_neighbours; //Closest neighbours

	//PBWT IBD2 protect
	vector < vector < IBD2track > > bannedPairs;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();
	void clear();

	//PBWT routines
	void parametrizePBWT(int, double, int, double, int);
	void initializePBWTmapping(variant_map &);
	void updatePBWTmapping();
	void allocatePBWTarrays();
	void selectPBWTarrays();
	void transposePBWTarrays();

	//IBD2 routines
	//void searchIBD2matching(genotype_set & G, variant_map & V, double minLengthIBDtrack, double windowSize, double ibd2_maf, double ibd2_mdr, int ibd2_count);
	//void writeIBD2matching(genotype_set & G, string);
	void mergeIBD2constraints();
	bool checkIBD2matching(int, int, int);

	//Haplotype routines
	void updateHaplotypes(genotype_set & G, bool first_time = false);
	void transposeHaplotypes_H2V(bool full);
	void transposeHaplotypes_V2H(bool full);
};

inline
bool haplotype_set::checkIBD2matching(int mh, int ch, int idx) {
	int mi = min(mh/2,ch/2);
	int ci = max(mh/2,ch/2);
	// Prevents self copying, who knows ?
	if (mi == ci) return false;
	// Prevents copying for IBD2 individuals
	for (int i = 0 ; i < bannedPairs[mi].size() && bannedPairs[mi][i].ind <= ci; i ++)
		if (bannedPairs[mi][i].ind==ci && bannedPairs[mi][i].from<=idx && bannedPairs[mi][i].to>=idx) return false;
	return true;
}

#endif
