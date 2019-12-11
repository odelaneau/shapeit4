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
#include <objects/compute_job.h>

bool compute_job::reccursive_window_splitting(double min_length_cm, int left_index, int right_index, vector < int > & idx_sta, vector < int > & idx_sto, vector < double > & ccm_sta, vector < double > & ccm_sto, vector < int > & output) {
	int number_of_segments = right_index-left_index+1;
	int number_of_variants = idx_sto[right_index] - idx_sta[left_index] + 1;
	double length_of_region = ccm_sto[right_index] - ccm_sta[left_index];

	//A phasing window must (i) span >=4 segments, (ii) contain >= 100 variants and (iii) span more than "min_length_cm" cM
	if (number_of_segments < 4 || number_of_variants < 100 || length_of_region < min_length_cm) return false;
	else {
		int split_point = rng.getInt(number_of_segments/2) + number_of_segments/4 + 1;
		vector <  int > left_output, right_output;
		bool ret1 = reccursive_window_splitting(min_length_cm, left_index, left_index + split_point, idx_sta, idx_sto, ccm_sta, ccm_sto, left_output);
		bool ret2 = reccursive_window_splitting(min_length_cm, left_index + split_point, right_index, idx_sta, idx_sto, ccm_sta, ccm_sto, right_output);

		if (ret1 && ret2) {
			//succesful split, so operate it
			output = vector < int >(left_output.size() + right_output.size());
			std::copy(left_output.begin(), left_output.end(), output.begin());
			std::copy(right_output.begin(), right_output.end(), output.begin() + left_output.size());
		} else {
			//unsuccesful split, so return current coordinates
			output.clear();
			output.push_back(left_index);
			output.push_back(right_index);
		}
		return true;
	}
}

compute_job::compute_job(variant_map & _V, genotype_set & _G, haplotype_set & _H, unsigned int n_max_transitions, unsigned int n_max_missing) : V(_V), G(_G), H(_H) {
	T = vector < double > (n_max_transitions, 0.0);
	M = vector < float > (n_max_missing , 0.0);
}

compute_job::~compute_job() {
	free();
}

void compute_job::free () {
	vector < double > ().swap(T);
	vector < coordinates > ().swap(C);
	vector < vector < unsigned int > > ().swap(Kvec);
}

void compute_job::reset() {
	C.clear();
	Kvec.clear();
}

void compute_job::make(unsigned int ind, double min_window_size) {
	//1. Mapping coordinates of each segment
	vector < unsigned int > loc_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > loc_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > amb_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > amb_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > mis_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > mis_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > tra_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > tra_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < double > ccm_sta = vector < double >(G.vecG[ind]->n_segments, 0);
	vector < double > ccm_sto = vector < double >(G.vecG[ind]->n_segments, 0);
	vector < int > idx_sta = vector < int >(G.vecG[ind]->n_segments, 0);
	vector < int > idx_sto = vector < int >(G.vecG[ind]->n_segments, 0);

	unsigned int prev_dipcounts = 1, curr_dipcounts = 0;
	for (unsigned int s = 0, a = 0, t = 0, v = 0, m = 0 ; s < G.vecG[ind]->n_segments ; s ++) {
		//update a
		amb_idx[s] = a;
		for (unsigned int vrel = 0 ; vrel < G.vecG[ind]->Lengths[s] ; vrel ++) amb_siz[s] += VAR_GET_AMB(MOD2(v+vrel), G.vecG[ind]->Variants[DIV2(v+vrel)]);
		a += amb_siz[s];
		//update m
		mis_idx[s] = m;
		for (unsigned int vrel = 0 ; vrel < G.vecG[ind]->Lengths[s] ; vrel ++) mis_siz[s] += VAR_GET_MIS(MOD2(v+vrel), G.vecG[ind]->Variants[DIV2(v+vrel)]);
		m += mis_siz[s];
		//update v
		loc_idx[s] = v;
		loc_siz[s] = G.vecG[ind]->Lengths[s];
		v += loc_siz[s];
		//update idx
		idx_sta[s] = loc_idx[s];
		idx_sto[s] = loc_idx[s]+loc_siz[s]-1;
		//update ccm
		ccm_sta[s] = V.vec_pos[idx_sta[s]]->cm;
		ccm_sto[s] = V.vec_pos[idx_sto[s]]->cm;
		//update t
		tra_idx[s] = t;
		curr_dipcounts = G.vecG[ind]->countDiplotypes(G.vecG[ind]->Diplotypes[s]);
		tra_siz[s] = prev_dipcounts * curr_dipcounts;
		t += tra_siz[s];
		prev_dipcounts = curr_dipcounts;
	}

	//2. Reccursive split
	vector < int > output;
	output.push_back(0);
	output.push_back(G.vecG[ind]->n_segments-1);
	reccursive_window_splitting(min_window_size, 0, G.vecG[ind]->n_segments-1, idx_sta, idx_sto, ccm_sta, ccm_sto, output);
/*
	for (int w =0 ; w < output.size() ; w += 2) {
		cout << w << " " << output[w] << " " << output[w+1] << endl;
	}
*/
	assert(output[0] == 0);
	assert(output.back() == G.vecG[ind]->n_segments-1);
	for (int w = 1 ; w < output.size()-1 ; w += 2) assert(output[w+0] == output[w+1]);
	int n_windows = output.size()/2;
	//cout << "N windows = " << n_windows << endl;

	//3. Update coordinates
	C = vector < coordinates > (n_windows);
	for (unsigned int w = 0 ; w < n_windows ; w ++) {
		C[w].start_segment = output[2*w+0];
		C[w].stop_segment = output[2*w+1];
		C[w].start_ambiguous = amb_idx[C[w].start_segment];
		C[w].stop_ambiguous = amb_idx[C[w].stop_segment] + amb_siz[C[w].stop_segment] - 1;

		C[w].start_missing = mis_idx[C[w].start_segment];
		C[w].stop_missing = mis_idx[C[w].stop_segment] + mis_siz[C[w].stop_segment] - 1;

		C[w].start_locus = loc_idx[C[w].start_segment];
		C[w].stop_locus = loc_idx[C[w].stop_segment] + loc_siz[C[w].stop_segment] - 1;
		C[w].start_transition = tra_idx[C[w].start_segment] + tra_siz[C[w].start_segment];
		C[w].stop_transition = tra_idx[C[w].stop_segment] + tra_siz[C[w].stop_segment] - 1;

		//cout << w << " " << C[w].toString() << endl;

	}
	assert(C.back().stop_ambiguous == G.vecG[ind]->n_ambiguous - 1);
	assert(C.back().stop_missing == G.vecG[ind]->n_missing - 1);
	assert(C.back().stop_segment == G.vecG[ind]->n_segments - 1);
	assert(C.back().stop_locus == G.vecG[ind]->n_variants - 1);
	assert(C.back().stop_transition == G.vecG[ind]->n_transitions - 1);
	//cout << "Done coordinates"<< endl;

	//4. Update conditional haps
	unsigned long addr_offset = H.pbwt_nstored * H.n_ind * 2UL;
	Kvec = vector < vector < unsigned int > > (n_windows);
	vector < int > phap = vector < int > (2 * H.pbwt_depth, -1);
	for (int l = 0, w = 0 ; l < H.pbwt_evaluated.size() ; l ++) {
		int abs_idx = H.pbwt_evaluated[l], rel_idx = H.pbwt_stored[l];
		if (abs_idx > C[w].stop_locus) { std::fill(phap.begin(), phap.end(), -1); w++; }
		if (rel_idx >= 0) {
			unsigned long curr_hap0 = 2*ind+0, curr_hap1 = 2*ind+1;
			bool addToNext = ((w+1)<n_windows && abs_idx>=C[w+1].start_locus);
			for (int s = 0 ; s < H.pbwt_depth ; s ++) {
				int cond_hap0 = H.pbwt_neighbours[s * addr_offset + curr_hap0*H.pbwt_nstored + rel_idx];
				int cond_hap1 = H.pbwt_neighbours[s * addr_offset + curr_hap1*H.pbwt_nstored + rel_idx];
				if (cond_hap0 != phap[2*s+0]) { Kvec[w].push_back(cond_hap0); phap[2*s+0] = cond_hap0; };
				if (cond_hap1 != phap[2*s+1]) { Kvec[w].push_back(cond_hap1); phap[2*s+1] = cond_hap1; };
				if (addToNext) { Kvec[w+1].push_back(cond_hap0); Kvec[w+1].push_back(cond_hap1); }
			}
		}
/*
		if (w == n_windows - 1) {
			cout << w << " " << l << " " << Kvec[w].size() << endl;
		}
*/
	}
	/*
	unsigned long block_size = H.pbwt_nstored * H.n_ind * 2UL;
	Kvec = vector < vector < unsigned int > > (n_windows);
	for (int d = 0 ; d < H.pbwt_depth ; d ++) {
		for (int h = 0; h < 2; h++) {
			unsigned long curr_hap = 2*ind+h;
			for (int l = 0, w = 0 ; l < H.pbwt_evaluated.size() ; l ++) {
				int abs_idx = H.pbwt_evaluated[l];
				int rel_idx = H.pbwt_stored[l];
				if (abs_idx > C[w].stop_locus) w++;
				if (rel_idx>=0) {
					int curr_cond_hap = H.pbwt_neighbours[d * block_size + curr_hap*H.pbwt_nstored + rel_idx];
					if (Kvec[w].empty() || curr_cond_hap != Kvec[w].back()) Kvec[w].push_back(curr_cond_hap);
					if ((w+1)<n_windows && abs_idx>=C[w+1].start_locus && (Kvec[w+1].empty() || curr_cond_hap != Kvec[w+1].back())) Kvec[w+1].push_back(curr_cond_hap);
				}
			}
		}
	}
	*/
	for (int w = 0 ; w < n_windows; w++) {
		sort(Kvec[w].begin(), Kvec[w].end());
		Kvec[w].erase(unique(Kvec[w].begin(), Kvec[w].end()), Kvec[w].end());
		/*
		cout << ind << " " << w << " " << Kvec[w].size()  << "======================================================="<< endl;
		if (w == n_windows - 1) {
			for (int k  = 0 ; k < Kvec[w].size() ; k ++) cout << " " << Kvec[w][k];
			cout << endl;
		}
		*/
	}
	//cout << "Done selection"<< endl;
}

void compute_job::maskingTransitions(unsigned int ind, double error_rate) {
	vector < double > curr_transitions = vector < double > (4096, 0.0);
	unsigned int prev_dipcount = 1, curr_dipcount = 0, curr_transcount = 0;
	for (unsigned int s = 0, t = 0 ; s < G.vecG[ind]->n_segments ; s ++) {
		curr_dipcount = G.vecG[ind]->countDiplotypes(G.vecG[ind]->Diplotypes[s]);
		curr_transcount = prev_dipcount * curr_dipcount;

		double sumT = 0.0;
		for (unsigned int trel = 0 ; trel < curr_transcount ; trel ++) {
			curr_transitions[trel] = T[t+trel] * (G.vecG[ind]->ProbabilityMask[t+trel]?(1.0-error_rate):(error_rate));
			sumT += curr_transitions[trel];
		}

		if (sumT > numeric_limits<double>::min() && !isnan(sumT))
			for (unsigned int trel = 0 ; trel < curr_transcount ; trel ++)
				T[t+trel] = curr_transitions[trel] / sumT;

		t += curr_transcount;
		prev_dipcount = curr_dipcount;
	}
}
