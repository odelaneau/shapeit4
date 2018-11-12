#include <containers/haplotype_set.h>

haplotype_set::haplotype_set() {
	n_site = 0;
	n_hap = 0;
	n_ind = 0;
	n_save = 0;
}

haplotype_set::~haplotype_set() {
	n_site = 0;
	n_hap = 0;
	n_ind = 0;
	n_save = 0;
	curr_clusters.clear();
	save_clusters.clear();
	dist_clusters.clear();
	abs_indexes.clear();
	rel_indexes.clear();
}

void haplotype_set::updateMapping() {
	rel_indexes.clear();
	int rint = rng.getInt(mod);
	for (int l = 0 ; l < abs_indexes.size() ; l ++) rel_indexes.push_back(((l%mod == rint)?(l/mod):(-1)));
}

void haplotype_set::allocate(variant_map & V, int _mod, int _depth) {
	mod = _mod;
	depth = _depth;
	for (int al = 0, rl = 0 ; al < n_site ; al ++) if (V.vec_pos[al]->getMAC() >= 2) {
		abs_indexes.push_back(al);
		rel_indexes.push_back(((rl%mod)?(-1):(rl/mod)));
		rl++;
	}
	n_save = (abs_indexes.size() / mod) + (abs_indexes.size() % mod != 0);
	save_clusters = vector < int > ((depth+1) * (unsigned long)n_save * n_ind * 2UL, 0);
	curr_clusters = vector < int > (2 * n_hap, 0);
	dist_clusters = vector < int > (2 * n_hap, 0);
}

void haplotype_set::update(genotype_set & G, bool first_time) {
	tac.clock();
	for (unsigned int i = 0 ; i < G.n_ind ; i ++) {
		for (unsigned int v = 0 ; v < n_site ; v ++) {
			if (first_time || (VAR_GET_AMB(MOD2(v), G.vecG[i]->Variants[DIV2(v)]))) {
				bool a0 = VAR_GET_HAP0(MOD2(v), G.vecG[i]->Variants[DIV2(v)]);
				bool a1 = VAR_GET_HAP1(MOD2(v), G.vecG[i]->Variants[DIV2(v)]);
				H_opt_hap.set(2*i+0, v, a0);
				H_opt_hap.set(2*i+1, v, a1);
			}
		}
	}
	vrb.bullet("HAP update (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeH2V(bool full) {
	tac.clock();
	if (!full) H_opt_hap.transpose(H_opt_var, 2*n_ind, n_site);
	else H_opt_hap.transpose(H_opt_var, n_hap, n_site);
	vrb.bullet("H2V transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeV2H(bool full) {
	tac.clock();
	if (!full) H_opt_var.transpose(H_opt_hap, n_site, 2*n_ind);
	else H_opt_var.transpose(H_opt_hap, n_site, n_hap);
	vrb.bullet("V2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::transposeC2H() {
	tac.clock();
	int block = 32;
	unsigned long addr_tar, addr_src;
	unsigned long addr_offset = n_save * (unsigned long)n_ind * 2UL;
	for (int d = 0; d < depth ; d ++) {
		for (int s = 0; s < n_save ; s += block) {
			for(int h = 0; h < n_ind * 2; ++h) {
				for(int b = 0; b < block && s + b < n_save; ++b) {
					addr_tar = depth * addr_offset + ((unsigned long)h)*n_save + s + b;
					addr_src = d * addr_offset + (s + b)*((unsigned long)n_ind)*2 + h;
					save_clusters[addr_tar] = save_clusters[addr_src];
				}
			}
		}
		std::copy(save_clusters.begin() + depth * addr_offset , save_clusters.end(), save_clusters.begin() + d * addr_offset );
	}
	vrb.bullet("C2H transpose (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::select() {
	tac.clock();
	updateMapping();
	int iprev = 1, inext = 0, i_added = 0;
	unsigned long addr_offset = n_save * (unsigned long)n_ind * 2UL;
	vector < int > B = vector < int >(n_hap, 0);
	vector < int > D = vector < int >(n_hap, 0);
	for (int l = 0 ; l < abs_indexes.size() ; l ++) {
		int u = 0, v = 0, p = l, q = l;
		int curr_block = abs_indexes[l] / lengthIBD2;
		int poffset = n_hap*(l%2), coffset = n_hap-poffset;
		iprev = inext;
		inext = 1 - iprev;
		for (int h = 0 ; h < n_hap ; h ++) {
			int alookup = l?curr_clusters[poffset+h]:h;
			int dlookup = l?dist_clusters[poffset+h]:0;
			if (dlookup > p) p = dlookup;
			if (dlookup > q) q = dlookup;
			if (!H_opt_var.get(abs_indexes[l], alookup)) {
				curr_clusters[coffset+u] = alookup;
				dist_clusters[coffset+u] = p;
				p = 0;
				u++;
			} else {
				B[v] = alookup;
				D[v] = q;
				q = 0;
				v++;
			}
		}
		std::copy(B.begin(), B.begin()+v,curr_clusters.begin()+coffset+u);
		std::copy(D.begin(), D.begin()+v,dist_clusters.begin()+coffset+u);
		if (rel_indexes[l] >= 0) {
			for (int h = 0 ; h < n_hap ; h ++) {
				int chap = curr_clusters[coffset+h];
				int cind = chap / 2;
				if (cind < n_ind) {
					unsigned long tar_idx = ((unsigned long)rel_indexes[l])*2*n_ind + chap;
					int offset0 = 1, offset1 = 1, lmatch0 = -1, lmatch1 = -1;
					bool add0 = false, add1 = false;
					for (int n_added = 0 ; n_added < depth ; ) {
						if ((h-offset0)>=0) {
							lmatch0 = max(dist_clusters[coffset+h-offset0+1], lmatch0);
							add0 = (curr_clusters[coffset+h-offset0]/2!=cind);
							if (flagIBD2[curr_block][cind]) add0 = (add0 && !banned(curr_block, cind, curr_clusters[coffset+h-offset0]/2));
						} else { add0 = false; lmatch0 = l; }
						if ((h+offset1)<n_hap) {
							lmatch1 = max(dist_clusters[coffset+h+offset1], lmatch1);
							add1 = (curr_clusters[coffset+h+offset1]/2!=cind);
							if (flagIBD2[curr_block][cind]) add1 = (add1 && !banned(curr_block, cind, curr_clusters[coffset+h+offset1]/2));
						} else { add1 = false; lmatch1 = l; }
						if (add0 && add1) {
							if (lmatch0 < lmatch1) {
								save_clusters[n_added*addr_offset+tar_idx] = curr_clusters[coffset+h-offset0];
								offset0++; n_added++;
							} else if (lmatch0 > lmatch1) {
								save_clusters[n_added*addr_offset+tar_idx] = curr_clusters[coffset+h+offset1];
								offset1++; n_added++;
							} else if (rng.flipCoin()) {
								save_clusters[n_added*addr_offset+tar_idx] = curr_clusters[coffset+h-offset0];
								offset0++; n_added++;
							} else {
								save_clusters[n_added*addr_offset+tar_idx] = curr_clusters[coffset+h+offset1];
								offset1++; n_added++;
							}
						} else if (add0) {
							save_clusters[n_added*addr_offset+tar_idx] = curr_clusters[coffset+h-offset0];
							offset0++; n_added++;
						} else if (add1) {
							save_clusters[n_added*addr_offset+tar_idx] = curr_clusters[coffset+h+offset1];
							offset1++; n_added++;
						} else {
							offset0++;
							offset1++;
						}
					}
				}
			}
		}
		vrb.progress("  * PBWT selection", (l+1)*1.0/abs_indexes.size());
	}
	vrb.bullet("PBWT selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void haplotype_set::searchIBD2(int _lengthIBD2) {
	tac.clock();
	lengthIBD2 = _lengthIBD2;
	flagIBD2 = vector < vector < bool > > (n_site / lengthIBD2 + 1, vector < bool > (n_ind, false));
	idxIBD2 = vector < vector < pair < int, int > > > (n_site / lengthIBD2 + 1, vector < pair < int, int > > ());
	int nBan = 0;
	vector < int > a0 = vector < int >(n_ind, 0);
	vector < int > a1 = vector < int >(n_ind, 0);
	vector < int > a2 = vector < int >(n_ind, 0);
	vector < int > d0 = vector < int >(n_ind, 0);
	vector < int > d1 = vector < int >(n_ind, 0);
	vector < int > d2 = vector < int >(n_ind, 0);
	vector < int > ban = vector < int >(n_ind, -1);
	for (int l = 0 ; l < n_site ; l ++) {
		int curr_block = l / lengthIBD2;
		int u = 0, v = 0, w = 0, p = l+1, q = l+1, r = l+1;
		for (int i = 0 ; i < n_ind ; i ++) {
			int alookup = l?a0[i]:i;
			int dlookup = l?d0[i]:0;
			if (dlookup > p) p = dlookup;
			if (dlookup > q) q = dlookup;
			if (dlookup > r) r = dlookup;
			unsigned char genotype = H_opt_var.get(l, 2*alookup+0) + H_opt_var.get(l, 2*alookup+1);
			switch (genotype) {
			case 0: a0[u] = alookup;
					d0[u] = p;
					p = 0;
					u++;
					break;
			case 1: a1[v] = alookup;
					d1[v] = q;
					q = 0;
					v++;
					break;
			case 2: a2[w] = alookup;
					d2[w] = r;
					r = 0;
					w++;
					break;
			}
		}
		copy(a1.begin(), a1.begin()+v,a0.begin()+u);
		copy(a2.begin(), a2.begin()+w,a0.begin()+u+v);
		copy(d1.begin(), d1.begin()+v,d0.begin()+u);
		copy(d2.begin(), d2.begin()+w,d0.begin()+u+v);

		if ((l%lengthIBD2) == (lengthIBD2 - 1)) {
			for (int i = 1 ; i < n_ind ; i ++) {
				int offset = 1;
				while (((i-offset)>=0) && ((l-d0[i-offset+1]+1)>=lengthIBD2)) {
					int ind0 = a0[i];
					int ind1 = a0[i-offset];
					flagIBD2[curr_block][ind0] = true;
					flagIBD2[curr_block][ind1] = true;
					idxIBD2[curr_block].push_back(pair <int, int > (min(ind0, ind1), max(ind0, ind1)));
					offset ++;
				}
			}
			sort(idxIBD2[curr_block].begin(), idxIBD2[curr_block].end());
			idxIBD2[curr_block].erase(unique(idxIBD2[curr_block].begin(), idxIBD2[curr_block].end()), idxIBD2[curr_block].end());
			nBan += idxIBD2[curr_block].size();
		}

		vrb.progress("  * IBD2 mask", (l+1)*1.0/n_site);
	}
	vrb.bullet("IBD2 mask [l=" + stb.str(lengthIBD2) + " / n=" + stb.str(nBan) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}


