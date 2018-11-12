#include <modules/pbwt_solver2.h>

pbwt_solver2::pbwt_solver2(haplotype_set & _H) : H(_H.H_opt_var) {
	n_total = 0;
	n_resolved = 0;
	n_site = _H.n_site;
	n_main_hap = 2 * _H.n_ind;
	n_total_hap = _H.n_hap;
	pbwt_clusters = vector < vector < int > > (2, vector < int > (n_total_hap, 0));
	pbwt_indexes = vector < vector < int > > (2, vector < int > (n_total_hap, 0));
	pbwt_divergences = vector < vector < int > > (2, vector < int > (n_total_hap, 0));
	Guess = vector < char > (n_total_hap, 0);
	Het = vector < bool > (n_main_hap/2, 0);
	Mis = vector < bool > (n_main_hap/2, 0);
	scoreBit = vector < float > (n_site, 0.0);
	for (int l = 0 ; l < n_site ; ++l) scoreBit[l] = log (l + 1.0);
	Amb = vector < vector < bool > > (n_site, vector < bool > (n_main_hap/2, false));
}

pbwt_solver2::~pbwt_solver2() {
	free();
}

void pbwt_solver2::free() {
	vector < vector < int > > ().swap(pbwt_clusters);
	vector < vector < int > > ().swap(pbwt_indexes);
	vector < vector < int > > ().swap(pbwt_divergences);
	vector < char > ().swap(Guess);
	vector < bool > ().swap(Het);
	vector < bool > ().swap(Mis);
}

void pbwt_solver2::backward(genotype_set & G) {
	tac.clock();
	vector < int > B = vector < int >(n_total_hap, 0);

	for (int l = n_site - 1 ; l >= 0 ; l--) {
		int idx_next = (l%2 == 1);
		int idx_prev = (l%2 == 0);
		if (l < n_site - 1) {
			double thresh = 2.5, s, s0, s1;
			unsigned int nh = 0;
			for (int h = 0 ; h < n_total_hap ; h++) Guess[h] = (H.get(l, h)?1:-1);
			for (int i = 0 ; i < n_main_hap/2 ; i ++) {
				Mis[i] = VAR_GET_MIS(MOD2(l), G.vecG[i]->Variants[DIV2(l)]);
				Het[i] = VAR_GET_HET(MOD2(l), G.vecG[i]->Variants[DIV2(l)]);
				Amb[l][i] = (Mis[i] || Het[i]);
				if (Amb[l][i]) {
					Guess[2*i+0] = 0;
					Guess[2*i+1] = 0;
				}
				nh+=Het[i];
			}
			//cout << na << endl;
			while (nh && thresh > 1.0) {
				int nhOld = nh;
				nh = 0;
				for (int i = 0, h = 0 ; i < n_main_hap/2 ; ++i, h += 2) {
					if (Amb[l][i]) {
						if (Het[i]) {
							unsigned int hidx0 = pbwt_indexes[idx_prev][h+0];
							unsigned int hidx1 = pbwt_indexes[idx_prev][h+1];
							s = 0.0;
							if (hidx0>0) s += Guess[pbwt_clusters[idx_prev][hidx0-1]];
							if (hidx0<(n_total_hap-1)) s += Guess[pbwt_clusters[idx_prev][hidx0+1]];
							if (hidx1>0) s -= Guess[pbwt_clusters[idx_prev][hidx1-1]];
							if (hidx1<(n_total_hap-1)) s -= Guess[pbwt_clusters[idx_prev][hidx1+1]];
							if (s > thresh) { Guess[h+0] = 1.0; Guess[h+1] = -1.0; Amb[l][i] = false; }
							else if (s < -thresh) { Guess[h+0] = -1.0; Guess[h+1] = 1.0; Amb[l][i] = false; }
							else ++nh;
						} else {
							unsigned int hidx0 = pbwt_indexes[idx_prev][h+0];
							unsigned int hidx1 = pbwt_indexes[idx_prev][h+1];
							if (hidx0>0) s0 = Guess[pbwt_clusters[idx_prev][hidx0-1]];
							if (hidx0<(n_total_hap-1)) s0 += Guess[pbwt_clusters[idx_prev][hidx0+1]];
							if (hidx1>0) s1 = Guess[pbwt_clusters[idx_prev][hidx1-1]];
							if (hidx1<(n_total_hap-1)) s1 += Guess[pbwt_clusters[idx_prev][hidx1+1]];
							if (s0 == -2 && s1 == -2) { Guess[h+0] = -1.0; Guess[h+1] = -1.0; Amb[l][i] = false; }
							else if (s0 == -2 && s1 == 2) { Guess[h+0] = -1.0; Guess[h+1] = 1.0; Amb[l][i] = false; }
							else if (s0 == 2 && s1 == -2) { Guess[h+0] = 1.0; Guess[h+1] = -1.0; Amb[l][i] = false; }
							else if (s0 == 2 && s1 == 2) { Guess[h+0] = 1.0; Guess[h+1] = 1.0; Amb[l][i] = false; }
						}
					}
				}
				if (nh == nhOld) thresh -= 1.0 ;
			}
			for (int h = 0 ; h < n_main_hap ; h++) if ((!Amb[l][h/2]) && (Mis[h/2]||Het[h/2])) H.set(l, h, Guess[h] > 0);
		}

		int u = 0, v = 0, p = l, q = l;
		for (int h = 0 ; h < n_total_hap ; h ++) {
			int alookup = (l<(n_site-1))?pbwt_clusters[idx_prev][h]:h;
			if (!H.get(l, alookup)) pbwt_clusters[idx_next][u++] = alookup;
			else B[v++] = alookup;
		}
		copy(B.begin(), B.begin()+v,pbwt_clusters[idx_next].begin()+u);
		for (int h = 0 ; h < n_total_hap ; h ++) pbwt_indexes[idx_next][pbwt_clusters[idx_next][h]] = h;
		vrb.progress("  * PBWT phase backward", (n_site-l+1)*1.0/n_site);
	}
	vrb.bullet("PBWT phase backward (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}


void pbwt_solver2::forward(genotype_set & G) {
	tac.clock();
	vector < int > B = vector < int >(n_total_hap, 0);
	vector < int > E = vector < int >(n_total_hap, 0);

	for (int l = 0 ; l < n_site ; l++) {
		int idx_next = (l%2 == 1);
		int idx_prev = (l%2 == 0);
		if (l > 0) {
			double thresh = 2.5, s, s0, s1;
			unsigned int nm = 0, nh = 0;
			for (int h = 0 ; h < n_total_hap ; h++) Guess[h] = (H.get(l, h)?1:-1);
			for (int i = 0 ; i < n_main_hap/2 ; i ++) {
				Mis[i] = VAR_GET_MIS(MOD2(l), G.vecG[i]->Variants[DIV2(l)]) && Amb[l][i];
				Het[i] = VAR_GET_HET(MOD2(l), G.vecG[i]->Variants[DIV2(l)]) && Amb[l][i];
				if (Het[i] || Mis[i]) { Guess[2*i+0] = 0; Guess[2*i+1] = 0;}
				nh+=Het[i];
				nm+=Mis[i];
			}
			while (nh && thresh > 1.0) {
				int nhOld = nh; nh = 0, nm = 0 ;
				for (int i = 0, h = 0 ; i < n_main_hap/2 ; ++i, h += 2) {
					if (Het[i]) {
						unsigned int hidx0 = pbwt_indexes[idx_prev][h+0];
						unsigned int hidx1 = pbwt_indexes[idx_prev][h+1];
						s = 0.0;
						if (hidx0>0) s += Guess[pbwt_clusters[idx_prev][hidx0-1]];
						if (hidx0<(n_total_hap-1)) s += Guess[pbwt_clusters[idx_prev][hidx0+1]];
						if (hidx1>0) s -= Guess[pbwt_clusters[idx_prev][hidx1-1]];
						if (hidx1<(n_total_hap-1)) s -= Guess[pbwt_clusters[idx_prev][hidx1+1]];
						if (s > thresh) { Guess[h+0] = 1.0; Guess[h+1] = -1.0; Amb[l][i] = false; }
						else if (s < -thresh) { Guess[h+0] = -1.0; Guess[h+1] = 1.0; Amb[l][i] = false; }
						else ++nh;
					} else if (Mis[i]) {
						unsigned int hidx0 = pbwt_indexes[idx_prev][h+0];
						unsigned int hidx1 = pbwt_indexes[idx_prev][h+1];
						if (hidx0>0) s0 = Guess[pbwt_clusters[idx_prev][hidx0-1]];
						if (hidx0<(n_total_hap-1)) s0 += Guess[pbwt_clusters[idx_prev][hidx0+1]];
						if (hidx1>0) s1 = Guess[pbwt_clusters[idx_prev][hidx1-1]];
						if (hidx1<(n_total_hap-1)) s1 += Guess[pbwt_clusters[idx_prev][hidx1+1]];
						if (s0 == -2 && s1 == -2) { Guess[h+0] = -1.0; Guess[h+1] = -1.0; Amb[l][i] = false; }
						else if (s0 == -2 && s1 == 2) { Guess[h+0] = -1.0; Guess[h+1] = 1.0; Amb[l][i] = false; }
						else if (s0 == 2 && s1 == -2) { Guess[h+0] = 1.0; Guess[h+1] = -1.0; Amb[l][i] = false; }
						else if (s0 == 2 && s1 == 2) { Guess[h+0] = 1.0; Guess[h+1] = 1.0; Amb[l][i] = false; }
						else ++nm;
					}
				}
				if (nh == nhOld) thresh -= 1.0 ;
			}
			if (nh || nm) {
				for (int i = 0, h = 0 ; i < n_main_hap/2 ; ++i, h += 2) {
					if (Amb[l][i]) {
						if (Het[i]) {
							unsigned int hidx0 = pbwt_indexes[idx_prev][h+0];
							unsigned int hidx1 = pbwt_indexes[idx_prev][h+1];
							s = 0.0;
							if (hidx0>0) s += Guess[pbwt_clusters[idx_prev][hidx0-1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx0] + 1];
							if (hidx0<(n_total_hap-1)) s += Guess[pbwt_clusters[idx_prev][hidx0+1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx0+1]+1];
							if (hidx1>0) s -= Guess[pbwt_clusters[idx_prev][hidx1-1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx1] + 1];
							if (hidx1<(n_total_hap-1)) s -= Guess[pbwt_clusters[idx_prev][hidx1+1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx1+1] + 1];
							if (s > 0) { Guess[h+0] = 1 ; Guess[h+1] = -1 ; }
							else { Guess[h+0] = -1 ; Guess[h+1] = 1 ; }
						}  else if (Mis[i]) {
							unsigned int hidx0 = pbwt_indexes[idx_prev][h+0];
							unsigned int hidx1 = pbwt_indexes[idx_prev][h+1];
							if (hidx0>0) s0 = Guess[pbwt_clusters[idx_prev][hidx0-1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx0] + 1];
							if (hidx0<(n_total_hap-1)) s0 += Guess[pbwt_clusters[idx_prev][hidx0+1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx0+1]+1];
							if (hidx1>0) s1 = Guess[pbwt_clusters[idx_prev][hidx1-1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx1] + 1];
							if (hidx1<(n_total_hap-1)) s1 += Guess[pbwt_clusters[idx_prev][hidx1+1]] * scoreBit[l - pbwt_divergences[idx_prev][hidx1+1] + 1];
							if (s0 > 0) Guess[h+0] = 1.0;
							else Guess[h+0] = -1.0;
							if (s1 > 0) Guess[h+1] = 1.0;
							else Guess[h+1] = -1.0;
						}
					}
				}
			}
			for (int h = 0 ; h < n_main_hap ; h++) if (Het[h/2] || Mis[h/2]) H.set(l, h, Guess[h] > 0);
		}

		int u = 0, v = 0, p = l, q = l;
		for (int h = 0 ; h < n_total_hap ; h ++) {
			int alookup = (l>0)?pbwt_clusters[idx_prev][h]:h;
			int dlookup = (l>0)?pbwt_divergences[idx_prev][h]:0;
			if (dlookup > p) p = dlookup;
			if (dlookup > q) q = dlookup;
			if (!H.get(l, alookup)) {
				pbwt_clusters[idx_next][u] = alookup;
				pbwt_divergences[idx_next][u] = p;
				p = 0;
				u++;
			} else {
				B[v] = alookup;
				E[v] = q;
				q = 0;
				v++;
			}
		}
		copy(B.begin(), B.begin()+v,pbwt_clusters[idx_next].begin()+u);
		copy(E.begin(), E.begin()+v,pbwt_divergences[idx_next].begin()+u);
		for (int h = 0 ; h < n_total_hap ; h ++) pbwt_indexes[idx_next][pbwt_clusters[idx_next][h]] = h;
		vrb.progress("  * PBWT phase forward", (l+1)*1.0/n_site);
	}
	vrb.bullet("PBWT phase forward (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
