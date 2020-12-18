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
#include <models/haplotype_segment_double.h>


haplotype_segment_double::haplotype_segment_double(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, coordinates & C, hmm_parameters & _M) : G(_G), M(_M){
	segment_first = C.start_segment;
	segment_last = C.stop_segment;
	locus_first = C.start_locus;
	locus_last = C.stop_locus;
	ambiguous_first = C.start_ambiguous;
	ambiguous_last = C.stop_ambiguous;
	missing_first = C.start_missing;
	missing_last = C.stop_missing;
	transition_first = C.start_transition;
	transition_last = C.stop_transition;
	n_cond_haps = idxH.size();
	n_missing = missing_last - missing_first + 1;

	probSumT = 0.0f;
	prob = vector < double > (HAP_NUMBER * n_cond_haps, 0.0f);
	probSumH = vector < double > (HAP_NUMBER, 0.0f);
	probSumK = vector < double > (n_cond_haps, 0.0f);
	Alpha = vector < vector < double > > (segment_last - segment_first + 1, vector < double > (HAP_NUMBER * n_cond_haps, 0.0f));
	AlphaSum = vector < vector < double > > (segment_last - segment_first + 1, vector < double > (HAP_NUMBER, 0.0f));
	AlphaSumSum = vector < double > (segment_last - segment_first + 1, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < vector < double > > (n_missing, vector < double > (HAP_NUMBER * n_cond_haps, 0.0f));
		AlphaSumMissing = vector < vector < double > > (n_missing, vector < double > (HAP_NUMBER, 0.0f));
	}

	//Cache efficient data transfer for conditioning haplotypes
	curr_rel_locus_offset = Hhap.subset(H, idxH, locus_first, locus_last);
	Hvar.allocateFast(Hhap.n_cols, Hhap.n_rows);
	Hhap.transpose(Hvar);
}

haplotype_segment_double::~haplotype_segment_double() {
	G = NULL;
	segment_first = 0;
	segment_last = 0;
	locus_first = 0;
	locus_last = 0;
	ambiguous_first = 0;
	ambiguous_last = 0;
	transition_first = 0;
	n_cond_haps = 0;
	curr_segment_index = 0;
	curr_segment_locus = 0;
	curr_abs_locus = 0;
	curr_rel_locus = 0;
	curr_abs_ambiguous = 0;
	curr_abs_transition = 0;
	probSumT = 0.0;
	prob.clear();
	probSumK.clear();
	probSumH.clear();
	Alpha.clear();
	AlphaSum.clear();
	AlphaSumSum.clear();
}

void haplotype_segment_double::forward() {
	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;

	for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool hom = !(amb || mis);

		if (curr_rel_locus == 0) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != 0) {
			if (hom) RUN_HOM(true);
			else if (amb) RUN_AMB(true);
			else RUN_MIS(true);
		} else {
			if (hom) COLLAPSE_HOM(true);
			else if (amb) COLLAPSE_AMB(true);
			else COLLAPSE_MIS(true);
		}

		if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK();
		if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) {
			Alpha[curr_segment_index - segment_first] = prob;
			AlphaSum[curr_segment_index - segment_first] = probSumH;
			AlphaSumSum[curr_segment_index - segment_first] = probSumT;
		}
		if (mis) {
			AlphaMissing[curr_rel_missing] = prob;
			AlphaSumMissing[curr_rel_missing] = probSumH;
			curr_abs_missing ++;
		}

		curr_segment_locus ++;
		curr_abs_ambiguous += amb;
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
}

int haplotype_segment_double::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities) {
	int n_underflow_recovered = 0;
	curr_segment_index = segment_last;
	curr_segment_locus = G->Lengths[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	curr_abs_transition = transition_last;

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool hom = !(amb || mis);

		if (curr_abs_locus == locus_last) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != G->Lengths[curr_segment_index] - 1) {
			if (hom) RUN_HOM(false);
			else if (amb) RUN_AMB(false);
			else RUN_MIS(false);
		} else {
			if (hom) COLLAPSE_HOM(false);
			else if (amb) COLLAPSE_AMB(false);
			else COLLAPSE_MIS(false);
		}
		if (curr_segment_locus == 0) SUMK();

		if (curr_abs_locus == 0) SET_FIRST_TRANS(transition_probabilities);
		if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (ret < 0) return ret;
			else n_underflow_recovered += ret;
		}

		if (mis) {
			IMPUTE(missing_probabilities);
			curr_abs_missing--;
		}


		curr_segment_locus--;
		curr_abs_ambiguous -= amb;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
	return n_underflow_recovered;
}

void haplotype_segment_double::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	double scale = 1.0f / probSumT, scaleDip = 0.0f;
	unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	vector < double > cprobs = vector < double > (n_transitions, 0.0);
	for (unsigned int d = 0, t = 0 ; d < 64 ; ++d) {
		if (DIP_GET(G->Diplotypes[0], d)) {
			cprobs[t] = (double)(probSumH[DIP_HAP0(d)]*scale) * (double)(probSumH[DIP_HAP1(d)]*scale);
			scaleDip += cprobs[t++];
		}
	}
	scaleDip = 1.0f / scaleDip;
	for (unsigned int t = 0 ; t < n_transitions ; t ++) transition_probabilities[t] = cprobs[t] * scaleDip;
}

int haplotype_segment_double::SET_OTHER_TRANS(vector < double > & transition_probabilities) {
	int underflow_recovered = 0;
	if (TRANS_HAP()) return -1;
	if (TRANS_DIP_MULT()) {
		if (TRANS_DIP_ADD()) return -2;
		else underflow_recovered = 1;
	}
	unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
	unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
	unsigned int n_transitions = curr_dipcount * prev_dipcount;
	double scaleDip = 1.0 / sumDProbs;
	curr_abs_transition -= (n_transitions - 1);
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaleDip;
	curr_abs_transition --;
	return underflow_recovered;
}
