#include <modules/pbwt_solver.h>

pbwt_solver::pbwt_solver(haplotype_set & _H) : H(_H.H_opt_var) {
	n_total = 0;
	n_resolved = 0;
	n_site = _H.n_site;
	n_main_hap = 2 * _H.n_ind;
	n_total_hap = _H.n_hap;
	pbwt_clusters = vector < vector < int > > (2, vector < int > (n_total_hap, 0));
	pbwt_indexes = vector < vector < int > > (2, vector < int > (n_total_hap, 0));
	vecA = vector < vector < bool > > (n_site, vector < bool > (n_main_hap/2, false));
	B = vector < int > (n_total_hap, 0);
	Guess = vector < bool > (n_total_hap, 0);
	GuessPrev = vector < int > (n_total_hap, 0);
	GuessNext = vector < int > (n_total_hap, 0);
	Het = vector < bool > (n_main_hap/2, 0);
	Mis = vector < bool > (n_main_hap/2, 0);
	Prob = vector < float > (n_main_hap*2, 0.0);
	double theta = 0.0;
	for (int h = 1 ; h < (n_total_hap - 1) ; h ++) theta += 1.0/h;
	theta = 1.0/theta;
	double ed = (theta/(n_total_hap+theta))*0.5;
	double ee = ed+(n_total_hap/(n_total_hap+theta));
	MUT[0] = ee; MUT[1] = (ee+ed)*0.5; MUT[2] = ed;
}

pbwt_solver::~pbwt_solver() {
	free();
}

void pbwt_solver::free() {
	vector < vector < int > > ().swap(pbwt_clusters);
	vector < vector < int > > ().swap(pbwt_indexes);
	vector < vector < bool > > ().swap(vecA);
	vector < int > ().swap(B);
	vector < bool > ().swap(Guess);
	vector < int > ().swap(GuessPrev);
	vector < int > ().swap(GuessNext);
	vector < bool > ().swap(Het);
	vector < bool > ().swap(Mis);
	vector < float > ().swap(Prob);
}

void pbwt_solver::forward(genotype_set & G, int PBWT_NBURN, int PBWT_NSAVE) {
	tac.clock();
	for (int l = 0 ; l < n_site ; l++) {
		int idx_next = (l%2 == 1);
		int idx_prev = (l%2 == 0);
		if (l > 0) {
			std::fill(Prob.begin(), Prob.end(), 0.0);
			for (int h = 0 ; h < n_total_hap ; h++) {
				int idx = pbwt_indexes[idx_prev][h];
				Guess[h] = (H.get(l, h)?1:0);
				GuessPrev[h] = idx?pbwt_clusters[idx_prev][idx-1]:-1;
				GuessNext[h] = (idx<(n_total_hap-1))?pbwt_clusters[idx_prev][idx+1]:-1;
			}
			for (int i = 0 ; i < n_main_hap/2 ; i ++) {
				Mis[i] = VAR_GET_MIS(MOD2(l), G.vecG[i]->Variants[DIV2(l)]);
				Het[i] = VAR_GET_HET(MOD2(l), G.vecG[i]->Variants[DIV2(l)]);
				vecA[l][i] = (Mis[i] || Het[i]);
				n_total += vecA[l][i];
			}
			for (int iter = 0 ; iter < (PBWT_NBURN + PBWT_NSAVE) ; iter++) {
				for (int h = 0, i = 0 ; h < n_main_hap ; h += 2, i++) {
					if (Het[i] || Mis[i]) {
						int hap0_prev = GuessPrev[h+0], hap0_next = GuessNext[h+0];
						int hap1_prev = GuessPrev[h+1], hap1_next = GuessNext[h+1];
						int hap0_allele = ((hap0_prev>=0)?Guess[hap0_prev]:0) + ((hap0_next>=0)?Guess[hap0_next]:0);
						int hap1_allele = ((hap1_prev>=0)?Guess[hap1_prev]:0) + ((hap1_next>=0)?Guess[hap1_next]:0);
						double hap0_probAis0 = MUT[hap0_allele];
						double hap1_probAis0 = MUT[hap1_allele];
						double hap0_probAis1 = MUT[2-hap0_allele];
						double hap1_probAis1 = MUT[2-hap1_allele];
						if (Mis[i]) {
							MP[0] = hap0_probAis0 * hap1_probAis0;
							MP[1] = hap0_probAis0 * hap1_probAis1;
							MP[2] = hap0_probAis1 * hap1_probAis0;
							MP[3] = hap0_probAis1 * hap1_probAis1;
							double sum = MP[0] + MP[1] + MP[2] + MP[3];
							switch(rng.sample4(MP, sum)) {
							case 0: Guess[h+0] = 0; Guess[h+1] = 0; break;
							case 1: Guess[h+0] = 0; Guess[h+1] = 1; break;
							case 2: Guess[h+0] = 1; Guess[h+1] = 0; break;
							case 3: Guess[h+0] = 1; Guess[h+1] = 1; break;
							}
							if (iter >= PBWT_NBURN) {
								Prob[h*2+0] += MP[0] / sum;
								Prob[h*2+1] += MP[1] / sum;
								Prob[h*2+2] += MP[2] / sum;
								Prob[h*2+3] += MP[3] / sum;
							}
						} else {
							double p01 = hap0_probAis0 * hap1_probAis1;
							double p10 = hap0_probAis1 * hap1_probAis0;
							double sum = p01 + p10;
							if (rng.getDouble() <= p01/sum) { Guess[h+0] = 0; Guess[h+1] = 1; }
							else { Guess[h+0] = 1; Guess[h+1] = 0; }
							if (iter >= PBWT_NBURN) {
								Prob[h*2+1] += p01 / sum;
								Prob[h*2+2] += p10 / sum;
							}
						}
					}
				}
			}
			for (int h = 0, i = 0 ; h < n_main_hap ; h+=2, i++) {
				if (Het[i]) {
					if (Prob[2*h+1] >= Prob[2*h+2]) { H.set(l,h,false); H.set(l,h+1,true); }
					else { H.set(l,h,true); H.set(l,h+1,false); }
					if (Prob[2*h+1] >= (0.99*PBWT_NSAVE) || Prob[2*h+2] >= (0.99*PBWT_NSAVE)) { vecA[l][i] = false; n_resolved ++; }
				} else if (Mis[i]) {
					int idxMax = distance(Prob.begin()+2*h, max_element(Prob.begin()+2*h, Prob.begin() + 2*h + 4));
					switch (idxMax) {
					case 0: H.set(l,h,false); H.set(l,h+1,false); break;
					case 1: H.set(l,h,false); H.set(l,h+1,true); break;
					case 2: H.set(l,h,true); H.set(l,h+1,false); break;
					case 3: H.set(l,h,true); H.set(l,h+1,true); break;
					}
					if (Prob[2*h+idxMax] >= (0.99*PBWT_NSAVE)) { vecA[l][i] = false; n_resolved++; }
				}
			}
		}
		int u = 0, v = 0;
		for (int h = 0 ; h < n_total_hap ; h ++) {
			int alookup = (l==0)?h:pbwt_clusters[idx_prev][h];
			if (!H.get(l, alookup)) pbwt_clusters[idx_next][u++] = alookup;
			else B[v++] = alookup;
		}
		copy(B.begin(), B.begin()+v,pbwt_clusters[idx_next].begin()+u);
		for (int h = 0 ; h < n_total_hap ; h ++) pbwt_indexes[idx_next][pbwt_clusters[idx_next][h]] = h;
		vrb.progress("  * PBWT phase forward", (l+1)*1.0/n_site);
	}
	vrb.bullet("PBWT phase forward [resolved=" + stb.str(n_resolved * 100.0 / n_total, 1) + "%] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void pbwt_solver::backward(genotype_set & G, int PBWT_NBURN, int PBWT_NSAVE) {
	tac.clock();
	for (int l = n_site - 1 ; l >= 0 ; l--) {
		int idx_next = (l%2 == 1);
		int idx_prev = (l%2 == 0);
		if (l < n_site - 1) {
			std::fill(Prob.begin(), Prob.end(), 0.0);
			for (int h = 0 ; h < n_total_hap ; h++) {
				int idx = pbwt_indexes[idx_next][h];
				Guess[h] = (H.get(l, h)?1:0);
				GuessPrev[h] = idx?pbwt_clusters[idx_next][idx-1]:-1;
				GuessNext[h] = (idx<(n_total_hap-1))?pbwt_clusters[idx_next][idx+1]:-1;
			}
			for (int i = 0 ; i < n_main_hap/2 ; i ++) {
				Mis[i] = VAR_GET_MIS(MOD2(l), G.vecG[i]->Variants[DIV2(l)]) && vecA[l][i];
				Het[i] = VAR_GET_HET(MOD2(l), G.vecG[i]->Variants[DIV2(l)]) && vecA[l][i];
			}
			for (int iter = 0 ; iter < (PBWT_NBURN + PBWT_NSAVE) ; iter++) {
				for (int h = 0, i = 0 ; h < n_main_hap ; h += 2, ++i) {
					if (Het[i] || Mis[i]) {
						int hap0_prev = GuessPrev[h+0], hap0_next = GuessNext[h+0];
						int hap1_prev = GuessPrev[h+1], hap1_next = GuessNext[h+1];
						int hap0_allele = ((hap0_prev>=0)?Guess[hap0_prev]:0) + ((hap0_next>=0)?Guess[hap0_next]:0);
						int hap1_allele = ((hap1_prev>=0)?Guess[hap1_prev]:0) + ((hap1_next>=0)?Guess[hap1_next]:0);
						double hap0_probAis0 = MUT[hap0_allele];
						double hap1_probAis0 = MUT[hap1_allele];
						double hap0_probAis1 = MUT[2-hap0_allele];
						double hap1_probAis1 = MUT[2-hap1_allele];
						if (Mis[i]) {
							MP[0] = hap0_probAis0 * hap1_probAis0;
							MP[1] = hap0_probAis0 * hap1_probAis1;
							MP[2] = hap0_probAis1 * hap1_probAis0;
							MP[3] = hap0_probAis1 * hap1_probAis1;
							double sum = MP[0] + MP[1] + MP[2] + MP[3];
							switch(rng.sample4(MP, sum)) {
							case 0: Guess[h+0] = 0; Guess[h+1] = 0; break;
							case 1: Guess[h+0] = 0; Guess[h+1] = 1; break;
							case 2: Guess[h+0] = 1; Guess[h+1] = 0; break;
							case 3: Guess[h+0] = 1; Guess[h+1] = 1; break;
							}
							if (iter >= PBWT_NBURN) {
								Prob[h*2+0] += MP[0] / sum;
								Prob[h*2+1] += MP[1] / sum;
								Prob[h*2+2] += MP[2] / sum;
								Prob[h*2+3] += MP[3] / sum;
							}
						} else {
							double p01 = hap0_probAis0 * hap1_probAis1;
							double p10 = hap0_probAis1 * hap1_probAis0;
							double sum = p01 + p10;
							if (rng.getDouble() <= p01/sum) { Guess[h+0] = 0; Guess[h+1] = 1; }
							else { Guess[h+0] = 1; Guess[h+1] = 0; }
							if (iter >= PBWT_NBURN) {
								Prob[h*2+1] += p01 / sum;
								Prob[h*2+2] += p10 / sum;
							}
						}
					}
				}
			}
			for (int h = 0, i = 0 ; h < n_main_hap ; h+=2, ++i) {
				if (Het[i]) {
					if (Prob[2*h+1] >= Prob[2*h+2]) { H.set(l,h,false); H.set(l,h+1,true); }
					else { H.set(l,h,true); H.set(l,h+1,false); }
				} else if (Mis[i]) {
					int idxMax = distance(Prob.begin()+2*h, max_element(Prob.begin()+2*h, Prob.begin() + 2*h + 4));
					switch (idxMax) {
					case 0: H.set(l,h,false); H.set(l,h+1,false); break;
					case 1: H.set(l,h,false); H.set(l,h+1,true); break;
					case 2: H.set(l,h,true); H.set(l,h+1,false); break;
					case 3: H.set(l,h,true); H.set(l,h+1,true); break;
					}
				}
			}
		}

		int u = 0, v = 0;
		for (int h = 0 ; h < n_total_hap ; h ++) {
			int alookup = (l==(n_site-1))?h:pbwt_clusters[idx_next][h];
			if (!H.get(l, alookup)) pbwt_clusters[idx_prev][u++] = alookup;
			else B[v++] = alookup;
		}
		copy(B.begin(), B.begin()+v,pbwt_clusters[idx_prev].begin()+u);
		for (int h = 0 ; h < n_total_hap ; h ++) pbwt_indexes[idx_prev][pbwt_clusters[idx_prev][h]] = h;
		vrb.progress("  * PBWT phase backward", (n_site-l+1)*1.0/n_site);
	}
	vrb.bullet("PBWT phase backward (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}


