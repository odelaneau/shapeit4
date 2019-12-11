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
#include <models/haplotype_segment.h>


haplotype_segment::haplotype_segment(genotype * _G, bitmatrix & _H, vector < unsigned int > & _idxH, coordinates & C, hmm_parameters & _M) : G(_G), H(_H), idxH(_idxH), M(_M) {
	segment_first = C.start_segment;
	segment_last = C.stop_segment;
	locus_first = C.start_locus;
	locus_last = C.stop_locus;
	ambiguous_first = C.start_ambiguous;
	ambiguous_last = C.stop_ambiguous;
	missing_first = C.start_missing;
	missing_last = C.stop_missing;
	transition_first = C.start_transition;
	n_cond_haps = idxH.size();
	prob1 = aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 1.0);
	prob2 = aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 1.0);
	probSumH1 = aligned_vector32 < float > (HAP_NUMBER, 1.0);
	probSumH2 = aligned_vector32 < float > (HAP_NUMBER, 1.0);
	probSumK1 = aligned_vector32 < float > (n_cond_haps, 1.0);
	probSumK2 = aligned_vector32 < float > (n_cond_haps, 1.0);
	probSumT1 = 1.0;
	probSumT2 = 1.0;
	Alpha = vector < aligned_vector32 < float > > (segment_last - segment_first + 1, aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0));
	Beta = vector < aligned_vector32 < float > > (segment_last - segment_first + 1, aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0));
	AlphaSum = vector < aligned_vector32 < float > > (segment_last - segment_first + 1, aligned_vector32 < float > (HAP_NUMBER, 0.0));
	AlphaSumSum = aligned_vector32 < float > (segment_last - segment_first + 1, 0.0);
	BetaSum = aligned_vector32 < float > (HAP_NUMBER, 0.0);
	n_missing = missing_last - missing_first + 1;
	//cout << G->name << " " << missing_first << " " << missing_last << " " << n_missing << endl;
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0));
		AlphaSumMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER, 0.0));
		ProbM1sums = aligned_vector32 < float >(HAP_NUMBER, 0.0);
		ProbM0sums = aligned_vector32 < float >(HAP_NUMBER, 0.0);
		ProbM = aligned_vector32 < float >(HAP_NUMBER, 0.0);
	}
}

haplotype_segment::~haplotype_segment() {
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
	curr_rel_segment_index = 0;
	curr_abs_ambiguous = 0;
	curr_abs_transition = 0;
	probSumT1 = 0.0;
	probSumT2 = 0.0;
	prob1.clear();
	prob2.clear();
	probSumK1.clear();
	probSumK2.clear();
	probSumH1.clear();
	probSumH2.clear();
	Alpha.clear();
	Beta.clear();
	AlphaSum.clear();
	AlphaSumSum.clear();
	BetaSum.clear();
}

void haplotype_segment::forward() {
	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;

	for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool paired = (curr_rel_locus % 2 == 0);
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);


		if (mis) MIS(paired);
		else if (amb) AMB(paired);
		else HOM(paired);

		if (curr_rel_locus != 0) {
			if (curr_segment_locus == 0) COLLAPSE(true, paired);
			else RUN(true, paired);
		}
		SUM(paired);
		if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK(paired);
		if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) {
			//if (paired) copy(prob2.begin(), prob2.end(), Alpha[curr_segment_index - segment_first].begin());
			//else copy(prob1.begin(), prob1.end(), Alpha[curr_segment_index - segment_first].begin());
			Alpha[curr_segment_index - segment_first] = (paired?prob2:prob1);		//does not compile with aligned vectors (needs to define operator = from std vector)
			AlphaSum[curr_segment_index - segment_first] = (paired?probSumH2:probSumH1);
			AlphaSumSum[curr_segment_index - segment_first] = (paired?probSumT2:probSumT1);
		}
		if (mis) {
			AlphaMissing[curr_rel_missing] = (paired?prob2:prob1);
			AlphaSumMissing[curr_rel_missing] = (paired?probSumH2:probSumH1);
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

void haplotype_segment::backward(vector < float > & missing_probabilities) {
	curr_segment_index = segment_last;
	curr_segment_locus = G->Lengths[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		//cout << curr_abs_missing << " " << missing_first << endl;
		//assert(curr_rel_missing >= 0);
		bool paired = (curr_rel_locus % 2 == 0);
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);

		if (mis) MIS(paired);
		else if (amb) AMB(paired);
		else HOM(paired);

		if (curr_abs_locus != locus_last) {
			if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) COLLAPSE(false, paired);
			else RUN(false, paired);
		}
		SUM(paired);
		if (curr_segment_locus == 0) SUMK(paired);
		if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
			Beta[curr_segment_index - segment_first] = (paired?prob2:prob1);	//does not compile with aligned vectors (needs to define operator = from std vector)
			//if (paired) copy(prob2.begin(), prob2.end(), Beta[curr_segment_index - segment_first].begin());
			//else copy(prob1.begin(), prob1.end(), Beta[curr_segment_index - segment_first].begin());
		}
		if (mis) {
			//Impute missing for 8 possible haplotypes
			fill(ProbM1sums.begin(), ProbM1sums.end(), 0.0);
			fill(ProbM0sums.begin(), ProbM0sums.end(), 0.0);
			for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
				for (int h = 0 ; h < HAP_NUMBER ; h ++) ProbM[h] = (AlphaMissing[curr_rel_missing][i + h] / AlphaSumMissing[curr_rel_missing][h]) * (paired?prob2[i + h]:prob1[i + h]);
				if (H.get(idxH[k], curr_abs_locus)) for (int h = 0 ; h < HAP_NUMBER ; h ++) ProbM1sums[h] += ProbM[h];
				else for (int h = 0 ; h < HAP_NUMBER ; h ++) ProbM0sums[h] += ProbM[h];
			}
			//cout << curr_abs_missing;
			for (int h = 0 ; h < HAP_NUMBER ; h ++) {
				missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = ProbM1sums[h] / (ProbM0sums[h]+ProbM1sums[h]);
				//cout << " " << stb.str(missing_probabilities[curr_abs_missing * HAP_NUMBER + h], 3);
			}
			//cout << endl;
			curr_abs_missing--;
		}

		if (curr_abs_locus == 0) BetaSum=(paired?probSumH2:probSumH1);
		curr_segment_locus--;
		curr_abs_ambiguous -= amb;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
}

int haplotype_segment::expectation(vector < double > & transition_probabilities, vector < float > & missing_probabilities) {
	//cout << "ok1 " << n_cond_haps << endl;
	forward();
	//cout << "ok2" << endl;
	backward(missing_probabilities);
	//cout << "ok3" << endl;

	unsigned int n_transitions = 0;
	if (!segment_first) {
		double sumHap = 0.0, sumDip = 0.0;
		n_transitions = G->countDiplotypes(G->Diplotypes[0]);
		for (int h = 0 ; h < HAP_NUMBER ; h ++) sumHap += BetaSum[h];
		vector < double > cprobs = vector < double > (n_transitions, 0.0);
		for (unsigned int d = 0, t = 0 ; d < 64 ; ++d) {
			if (DIP_GET(G->Diplotypes[0], d)) {
				cprobs[t] = (double)(BetaSum[DIP_HAP0(d)]/sumHap) * (double)(BetaSum[DIP_HAP1(d)]/sumHap);
				sumDip += cprobs[t];
				t++;
			}
		}
		for (unsigned int t = 0 ; t < n_transitions ; t ++) transition_probabilities[t] = (cprobs[t] / sumDip);
		curr_abs_transition += n_transitions;
	}

	unsigned int curr_abs_transition = transition_first;
	unsigned int curr_dipcount = 0, prev_dipcount = G->countDiplotypes(G->Diplotypes[segment_first]);

	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	int n_underflow_recovered = 0;
	for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus ++) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_segment_index = curr_segment_index - segment_first;

		if (curr_rel_locus != 0 && curr_segment_locus == 0) {
			if (TRANSH()) return -1;
			if (TRANSD(n_underflow_recovered)) return -1;
			curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
			n_transitions = curr_dipcount * prev_dipcount;
			double scaling = 1.0 / sumDProbs;

			//Unrolling of: // for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaling;
			int t = 0, repeat = (n_transitions / 4), left = (n_transitions % 4);
			while (repeat --) {
				transition_probabilities[curr_abs_transition + t + 0] = DProbs[t+0] * scaling;
				transition_probabilities[curr_abs_transition + t + 1] = DProbs[t+1] * scaling;
				transition_probabilities[curr_abs_transition + t + 2] = DProbs[t+2] * scaling;
				transition_probabilities[curr_abs_transition + t + 3] = DProbs[t+3] * scaling;
				t += 4;
			}
			switch (left) {
			case 3:	transition_probabilities[curr_abs_transition + t + 2] = DProbs[t+2] * scaling;
			case 2:	transition_probabilities[curr_abs_transition + t + 1] = DProbs[t+1] * scaling;
			case 1:	transition_probabilities[curr_abs_transition + t + 0] = DProbs[t+0] * scaling;
			}
			//Fin unrolling


			curr_abs_transition += n_transitions;
			prev_dipcount = curr_dipcount;
		}

		curr_segment_locus ++;
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
	//cout << "ok4" << endl;
	return n_underflow_recovered;
}



