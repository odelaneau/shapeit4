////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////
#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	clear();
}

haplotype_set::~haplotype_set() {
	clear();
}

void haplotype_set::clear() {
	n_site = 0;
	n_hap = 0;
	n_ind = 0;
	pbwt_modulo = 0.0;
	pbwt_depth = 0;
	pbwt_mac = 0;
	pbwt_mdr = 0.0;
	pbwt_nstored = 0;
	nthreads = 0;
	pbwt_evaluated.clear();
	pbwt_stored.clear();
	pbwt_parray.clear();
	pbwt_darray.clear();
	pbwt_neighbours.clear();
}

void haplotype_set::parametrizePBWT(int _pbwt_depth, double _pbwt_modulo, int _pbwt_mac, double _pbwt_mdr, int _nthreads) {
	pbwt_modulo = _pbwt_modulo;
	pbwt_depth = _pbwt_depth;
	pbwt_mac = _pbwt_mac;
	pbwt_mdr = _pbwt_mdr;
	nthreads = _nthreads;
}

void haplotype_set::initializePBWTmapping(variant_map & V) {
	tac.clock();
	for (int l = 0 ; l < n_site ; l ++) {
		if (V.vec_pos[l]->getMAC() >= pbwt_mac && V.vec_pos[l]->getMDR() <= pbwt_mdr) {
			pbwt_evaluated.push_back(l);
			pbwt_cm.push_back(V.vec_pos[l]->cm);
			pbwt_grp.push_back((int)round(V.vec_pos[l]->cm / pbwt_modulo));
		}
	}
	for (int l = 0, src = -1, tar = -1 ; l < pbwt_grp.size() ; l ++) {
		if (src == pbwt_grp[l]) pbwt_grp[l] = tar;
		else {
			src = pbwt_grp[l];
			pbwt_grp[l] = ++tar;
		}
	}
	pbwt_nstored = pbwt_grp.back() + 1;
	updatePBWTmapping();
	vrb.bullet("PBWT indexing [l=" + stb.str(pbwt_nstored) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::updatePBWTmapping() {
	pbwt_stored = vector < int > (pbwt_grp.size() , -1);
	for (int idx = 0, loffset = 0, liter = 0 ; idx <= pbwt_grp.back() ; idx ++) {
		for (liter = 0; (loffset+liter) < pbwt_grp.size() && (pbwt_grp[loffset+liter]==idx) ; liter++);
		pbwt_stored[loffset + rng.getInt(liter)] = idx;
		loffset += liter;
	}
}

void haplotype_set::allocatePBWTarrays() {
	assert(pbwt_evaluated.size() > 0);
	pbwt_neighbours = vector < int > ((pbwt_depth+1) * pbwt_nstored * n_ind * 2UL, 0);	//The depth+1 is there for transposing the matrix
	pbwt_parray = vector < int > (n_hap, 0);
	pbwt_darray = vector < int > (n_hap, 0);
}

void haplotype_set::updateHaplotypes(genotype_set & G, bool first_time) {
	tac.clock();
	for (unsigned int i = 0 ; i < G.n_ind ; i ++) {
		for (unsigned int v = 0 ; v < n_site ; v ++) {
			if (first_time || (VAR_GET_HET(MOD2(v), G.vecG[i]->Variants[DIV2(v)])) || (VAR_GET_MIS(MOD2(v), G.vecG[i]->Variants[DIV2(v)]))) {
				bool a0 = VAR_GET_HAP0(MOD2(v), G.vecG[i]->Variants[DIV2(v)]);
				bool a1 = VAR_GET_HAP1(MOD2(v), G.vecG[i]->Variants[DIV2(v)]);
				H_opt_hap.set(2*i+0, v, a0);
				H_opt_hap.set(2*i+1, v, a1);
			}
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_H2V(bool full) {
	tac.clock();
	if (!full) H_opt_hap.transpose(H_opt_var, 2*n_ind, n_site);
	else H_opt_hap.transpose(H_opt_var, n_hap, n_site);
	vrb.bullet("H2V transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeHaplotypes_V2H(bool full) {
	tac.clock();
	if (!full) H_opt_var.transpose(H_opt_hap, n_site, 2*n_ind);
	else H_opt_var.transpose(H_opt_hap, n_site, n_hap);
	vrb.bullet("V2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposePBWTarrays() {
	tac.clock();
	int block = 32;
	unsigned long addr_tar, addr_src;
	unsigned long addr_offset = pbwt_nstored * n_ind * 2UL;
	for (int d = 0; d < pbwt_depth ; d ++) {
		for (int s = 0; s < pbwt_nstored ; s += block) {
			for(int h = 0; h < n_ind * 2; ++h) {
				for(int b = 0; b < block && s + b < pbwt_nstored; ++b) {
					addr_tar = pbwt_depth * addr_offset + h*pbwt_nstored + s + b;
					addr_src = d * addr_offset + (s + b)*n_ind*2UL + h;
					pbwt_neighbours[addr_tar] = pbwt_neighbours[addr_src];
				}
			}
		}
		std::copy(pbwt_neighbours.begin() + pbwt_depth * addr_offset , pbwt_neighbours.end(), pbwt_neighbours.begin() + d * addr_offset );
	}
	vrb.bullet("C2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::selectPBWTarrays() {
	tac.clock();
	if (bannedPairs.size() == 0) bannedPairs = vector < vector < IBD2track > > (n_ind);
	vector < int > B = vector < int > (n_hap, 0);
	vector < int > D = vector < int > (n_hap, 0);
	for (int h = 0 ; h < n_hap ; h ++) pbwt_parray[h] = h;
	fill(pbwt_darray.begin(), pbwt_darray.end(), 0);
	for (int l = 0 ; l < pbwt_evaluated.size() ; l ++) {
		int u = 0, v = 0, p = l, q = l;

		//PBWT PASS
		for (int h = 0 ; h < n_hap ; h ++) {
			int alookup = pbwt_parray[h];
			int dlookup = pbwt_darray[h];
			if (dlookup > p) p = dlookup;
			if (dlookup > q) q = dlookup;
			if (!H_opt_var.get(pbwt_evaluated[l], alookup)) {
				pbwt_parray[u] = alookup;
				pbwt_darray[u] = p;
				p = 0;
				u++;
			} else {
				B[v] = alookup;
				D[v] = q;
				q = 0;
				v++;
			}
		}
		std::copy(B.begin(), B.begin()+v, pbwt_parray.begin()+u);
		std::copy(D.begin(), D.begin()+v, pbwt_darray.begin()+u);

		//PBWT STORAGE
		if (pbwt_stored[l] >= 0) {
			unsigned long addr_offset = pbwt_nstored * n_ind * 2UL;
			for (int h = 0 ; h < n_hap ; h ++) {
				int chap = pbwt_parray[h];
				int cind = chap / 2;
				if (cind < n_ind) {
					int add_guess0 = 0, add_guess1 = 0, offset0 = 1, offset1 = 1, hap_guess0 = -1, hap_guess1 = -1, div_guess0 = -1, div_guess1 = -1;
					unsigned long tar_idx = pbwt_stored[l] * 2UL * n_ind + chap;
					for (int n_added = 0 ; n_added < pbwt_depth ; ) {
						if ((h-offset0)>=0) {
							hap_guess0 = pbwt_parray[h-offset0];
							div_guess0 = max(pbwt_darray[h-offset0+1], div_guess0);
							add_guess0 = checkIBD2matching(chap, hap_guess0, pbwt_evaluated[l]);
							//cout << "AG0= " << hap_guess0 << " " << div_guess0 << " " << add_guess0 << endl;
						} else { add_guess0 = 0; div_guess0 = l+1; }
						if ((h+offset1)<n_hap) {
							hap_guess1 = pbwt_parray[h+offset1];
							div_guess1 = max(pbwt_darray[h+offset1], div_guess1);
							add_guess1 = checkIBD2matching(chap, hap_guess1, pbwt_evaluated[l]);
							//cout << "AG1= " << hap_guess0 << " " << div_guess0 << " " << add_guess0 << endl;
						} else { add_guess1 = 0; div_guess1 = l+1; }
						if (add_guess0 && add_guess1) {
							if (div_guess0 < div_guess1) {
								pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess0;
								offset0++; n_added++;
							} else {
								pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess1;
								offset1++; n_added++;
							}
						} else if (add_guess0) {
							pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess0;
							offset0++; n_added++;
						} else if (add_guess1) {
							pbwt_neighbours[n_added*addr_offset+tar_idx] = hap_guess1;
							offset1++; n_added++;
						} else {
							offset0++;
							offset1++;
						}
					}
				}
			}
		}
		vrb.progress("  * PBWT selection", (l+1)*1.0/pbwt_evaluated.size());
	}
	vrb.bullet("PBWT selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::mergeIBD2constraints() {
	unsigned int n_inds_with_ibd2 = 0;
	unsigned int n_ibd2_blocs = 0;
	unsigned int n_ibd2_merged = 0;
	for (int i = 0 ; i < n_ind ; i ++) {
		sort(bannedPairs[i].begin(), bannedPairs[i].end());
		if (bannedPairs[i].size() > 1) {
			for(vector < IBD2track > :: iterator it = bannedPairs[i].begin() + 1 ; it != bannedPairs[i].end() ; ) {
				if (it->overlap(*(it-1))) {
					(it-1)->merge(*it);
					it = bannedPairs[i].erase(it);
					n_ibd2_merged++;
				} else ++it;
			}
		}
		n_ibd2_blocs += bannedPairs[i].size();
		n_inds_with_ibd2+= (bannedPairs[i].size() > 0);
	}
	vrb.bullet("IBD2 constraints [#inds=" + stb.str(n_inds_with_ibd2) + " / #contraints=" + stb.str(n_ibd2_blocs) + " / #merged = " + stb.str(n_ibd2_merged) + "]");
}
