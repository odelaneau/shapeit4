#include <objects/compute_job.h>

void split(unsigned int min_segments, unsigned int leftB, unsigned int rightB, vector < unsigned int > & output) {
	output.clear();
	unsigned int n_curr_segments = rightB-leftB+1;
	if (n_curr_segments < 2 * min_segments) {
		output.push_back(leftB);
		output.push_back(rightB);
	} else {
		unsigned int split_point = rng.getInt(n_curr_segments - 2 * min_segments + 1) + min_segments;
		vector < unsigned int > left_output, right_output;
		split(min_segments, leftB, leftB + split_point, left_output);
		split(min_segments, leftB + split_point, rightB, right_output);
		output = vector < unsigned int >(left_output.size() + right_output.size());
		std::copy(left_output.begin(), left_output.end(), output.begin());
		std::copy(right_output.begin(), right_output.end(), output.begin() + left_output.size());
	}
}

compute_job::compute_job(variant_map & _V, genotype_set & _G, haplotype_set & _H, unsigned int n_max_transitions) : V(_V), G(_G), H(_H) {
	T = vector < double > (n_max_transitions, 0.0);
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
	unsigned int n_segments_per_window, n_windows;
	unsigned int n_splits = (unsigned int)round(V.length() * 1.0 / min_window_size);
	if (!n_splits) n_splits = 1;
	n_segments_per_window = G.vecG[ind]->n_segments/n_splits;

	//Recursive split into overalaping windows
	vector < unsigned int > output = vector < unsigned int > (2, 0); output[1] = G.vecG[ind]->n_segments -1;
	if (n_segments_per_window >= 2) split(n_segments_per_window, 0, G.vecG[ind]->n_segments-1, output);
	n_windows = output.size()/2;

	//Map coordinates of each segment
	vector < unsigned int > loc_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > loc_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > amb_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > amb_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > tra_idx = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	vector < unsigned int > tra_siz = vector < unsigned int >(G.vecG[ind]->n_segments, 0);
	unsigned int prev_dipcounts = 1, curr_dipcounts = 0;
	for (unsigned int s = 0, a = 0, t = 0, v = 0 ; s < G.vecG[ind]->n_segments ; s ++) {
		//update a
		amb_idx[s] = a;
		for (unsigned int vrel = 0 ; vrel < G.vecG[ind]->Lengths[s] ; vrel ++) amb_siz[s] += VAR_GET_AMB(MOD2(v+vrel), G.vecG[ind]->Variants[DIV2(v+vrel)]);
		a += amb_siz[s];
		//update v
		loc_idx[s] = v;
		loc_siz[s] = G.vecG[ind]->Lengths[s];
		v += loc_siz[s];
		//update t
		tra_idx[s] = t;
		curr_dipcounts = G.vecG[ind]->countDiplotypes(G.vecG[ind]->Diplotypes[s]);
		tra_siz[s] = prev_dipcounts * curr_dipcounts;
		t += tra_siz[s];
		prev_dipcounts = curr_dipcounts;
	}

	//Update coordinates
	C = vector < coordinates > (n_windows);
	for (unsigned int w = 0 ; w < n_windows ; w ++) {
		C[w].start_segment = output[2*w+0];
		C[w].stop_segment = output[2*w+1];
		C[w].start_ambiguous = amb_idx[C[w].start_segment];
		C[w].stop_ambiguous = amb_idx[C[w].stop_segment] + amb_siz[C[w].stop_segment] - 1;
		C[w].start_locus = loc_idx[C[w].start_segment];
		C[w].stop_locus = loc_idx[C[w].stop_segment] + loc_siz[C[w].stop_segment] - 1;
		C[w].start_transition = tra_idx[C[w].start_segment] + tra_siz[C[w].start_segment];
		C[w].stop_transition = tra_idx[C[w].stop_segment] + tra_siz[C[w].stop_segment] - 1;
	}

	assert(C.back().stop_ambiguous == G.vecG[ind]->n_ambiguous - 1);
	assert(C.back().stop_segment == G.vecG[ind]->n_segments - 1);
	assert(C.back().stop_locus == G.vecG[ind]->n_variants - 1);
	assert(C.back().stop_transition == G.vecG[ind]->n_transitions - 1);

	//Update conditional haps
	int n_state = H.depth;
	unsigned long addr_offset = H.n_save * (unsigned long)H.n_ind * 2UL;
	Kvec = vector < vector < unsigned int > > (n_windows);
	vector < int > phap = vector < int > (2 * H.depth, -1);
	for (int l = 0, w = 0 ; l < H.abs_indexes.size() ; l ++) {
		int abs_idx = H.abs_indexes[l];
		int rel_idx = H.rel_indexes[l];
		if (abs_idx > C[w].stop_locus) { std::fill(phap.begin(), phap.end(), -1); w++; }
		bool addToNext = ((w+1)<n_windows && abs_idx>=C[w+1].start_locus);
		if (rel_idx >= 0) {
			unsigned long curr_hap0 = 2*ind+0;
			unsigned long curr_hap1 = 2*ind+1;
			for (int s = 0 ; s < H.depth ; s ++) {
				int cond_hap0 = H.save_clusters[s * addr_offset + curr_hap0*H.n_save + rel_idx];
				int cond_hap1 = H.save_clusters[s * addr_offset + curr_hap1*H.n_save + rel_idx];
				if (cond_hap0 != phap[2*s+0]) { Kvec[w].push_back(cond_hap0); phap[2*s+0] = cond_hap0; };
				if (cond_hap1 != phap[2*s+1]) { Kvec[w].push_back(cond_hap1); phap[2*s+1] = cond_hap1; };
				if (addToNext) { Kvec[w+1].push_back(cond_hap0); Kvec[w+1].push_back(cond_hap1); }
			}
		}
	}
	for (int w = 0 ; w < n_windows; w++) {
		sort(Kvec[w].begin(), Kvec[w].end());
		Kvec[w].erase(unique(Kvec[w].begin(), Kvec[w].end()), Kvec[w].end());
	}
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
