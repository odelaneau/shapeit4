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

class haplotype_set {
public:
	//DATA
	unsigned long n_site, n_hap, n_ind, n_save, mod, lengthIBD2, depth;
	bitmatrix H_opt_hap;
	bitmatrix H_opt_var;
	vector < int > abs_indexes, rel_indexes;
	vector < int > curr_clusters, dist_clusters, save_clusters;

	//IBD2
	vector < vector < bool > > flagIBD2;
	vector < vector < pair < int, int > > > idxIBD2;

	//CONSTRUCTOR/DESTRUCTOR/INITIALIZATION
	haplotype_set();
	~haplotype_set();

	//ROUTINES
	void allocate(variant_map &, int, int);
	void update(genotype_set & G, bool first_time = false);
	void select();
	void transposeH2V(bool full);
	void transposeV2H(bool full);
	void transposeC2H();
	void updateMapping();

	void searchIBD2(int);
	bool banned(int, int, int);
};

inline
bool haplotype_set::banned(int b, int _i0, int _i1) {
	int i0 = min(_i0, _i1);
	int i1 = max(_i0, _i1);
	for (int c = 0 ; c < idxIBD2[b].size() ; c++) if (idxIBD2[b][c].first == i0 && idxIBD2[b][c].second == i1) return true;
	return false;
}

#endif
