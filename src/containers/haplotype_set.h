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
