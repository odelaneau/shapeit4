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
#ifndef _HAPLOTYPE_SEGMENT_DOUBLE_H
#define _HAPLOTYPE_SEGMENT_DOUBLE_H

#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>

class haplotype_segment_double {
private:
	//EXTERNAL DATA
	hmm_parameters & M;
	genotype * G;
	bitmatrix Hhap, Hvar;

	//COORDINATES & CONSTANTS
	int segment_first;
	int segment_last;
	int locus_first;
	int locus_last;
	int ambiguous_first;
	int ambiguous_last;
	int missing_first;
	int missing_last;
	int transition_first;
	int transition_last;
	unsigned int n_cond_haps;
	unsigned int n_missing;

	//CURSORS
	int curr_segment_index;
	int curr_segment_locus;
	int curr_abs_locus;
	int curr_rel_locus;
	int curr_rel_locus_offset;
	int curr_abs_ambiguous;
	int curr_abs_transition;
	int curr_abs_missing;
	int curr_rel_missing;

	//DYNAMIC ARRAYS
	double probSumT;
	vector < double > prob;
	vector < double > probSumK;
	vector < double > probSumH;
	vector < vector < double > > Alpha;
	vector < vector < double > > AlphaSum;
	vector < double > AlphaSumSum;
	vector < vector < double > > AlphaMissing;
	vector < vector < double > > AlphaSumMissing;
	double HProbs [HAP_NUMBER * HAP_NUMBER];
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER];

	//STATIC ARRAYS
	double sumHProbs;
	double sumDProbs;
	double g0[HAP_NUMBER], g1[HAP_NUMBER];

	//INLINED AND UNROLLED ROUTINES
	void INIT_HOM();
	void INIT_AMB();
	void INIT_MIS();
	void RUN_HOM(bool);
	void RUN_AMB(bool);
	void RUN_MIS(bool);
	void COLLAPSE_HOM(bool);
	void COLLAPSE_AMB(bool);
	void COLLAPSE_MIS(bool);
	void SUMK();
	void IMPUTE(vector < float > & );
	bool TRANS_HAP();
	bool TRANS_DIP_MULT();
	bool TRANS_DIP_ADD();
	void SET_FIRST_TRANS(vector < double > & );
	int SET_OTHER_TRANS(vector < double > & );

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment_double(genotype *, bitmatrix &, vector < unsigned int > &, coordinates &, hmm_parameters &);
	~haplotype_segment_double();

	//void fetch();
	void forward();
	int backward(vector < double > &, vector < float > &);
};

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE			************************/
/*******************************************************************************/

inline
void haplotype_segment_double::INIT_HOM() {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (ag==ah)?1.0f:M.ed/M.ee);
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::RUN_HOM(bool forward) {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	vector < double > _tFreq = vector < double >(HAP_NUMBER, M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	for (int h = 0 ; h < HAP_NUMBER ; h++) _tFreq[h] *= probSumH[h];
	double _nt = M.nt[curr_abs_locus-forward] / probSumT;
	double _mismatch = M.ed/M.ee;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] = (prob[i+h]*_nt)+_tFreq[h];
		if (ag!=ah) for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] *= _mismatch;
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::COLLAPSE_HOM(bool forward) {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	double _tFreq = M.t[curr_abs_locus-forward] / n_cond_haps;
	double _nt = M.nt[curr_abs_locus-forward] / probSumT;
	double _mismatch = M.ed/M.ee;
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (ag==ah)?((probSumK[k]*_nt)+_tFreq):(((probSumK[k]*_nt)+_tFreq)*_mismatch));
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE			********************/
/*******************************************************************************/

inline
void haplotype_segment_double::INIT_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		if (ah) memcpy(&prob[i], &g1[0], HAP_NUMBER*sizeof(double));
		else memcpy(&prob[i], &g0[0], HAP_NUMBER*sizeof(double));
		for (int h = 0 ; h < HAP_NUMBER ; h ++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
inline
void haplotype_segment_double::RUN_AMB(bool forward) {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	vector < double > _tFreq = vector < double >(HAP_NUMBER, M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	for (int h = 0 ; h < HAP_NUMBER ; h++) _tFreq[h] *= probSumH[h];
	double _nt = M.nt[curr_abs_locus-forward] / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] = (prob[i+h]*_nt)+_tFreq[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] *= ah?g1[h]:g0[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::COLLAPSE_AMB(bool forward) {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	double _tFreq = M.t[curr_abs_locus-forward] / n_cond_haps;
	double _nt = M.nt[curr_abs_locus-forward] / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (probSumK[k]*_nt)+_tFreq);
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] *= ah?g1[h]:g0[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			MISSING GENOTYPE			************************/
/*******************************************************************************/

inline
void haplotype_segment_double::INIT_MIS() {
	fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * n_cond_haps));
	fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
	probSumT = 1.0f;
}

inline
void haplotype_segment_double::RUN_MIS(bool forward) {
	vector < double > _tFreq = vector < double >(HAP_NUMBER, M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT));
	for (int h = 0 ; h < HAP_NUMBER ; h++) _tFreq[h] *= probSumH[h];
	double _nt = M.nt[curr_abs_locus-forward] / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		for (int h = 0 ; h < HAP_NUMBER ; h++) prob[i+h] = (prob[i+h]*_nt)+_tFreq[h];
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::COLLAPSE_MIS(bool forward) {
	double _tFreq = M.t[curr_abs_locus-forward] / n_cond_haps;
	double _nt = M.nt[curr_abs_locus-forward] / probSumT;
	fill(probSumH.begin(), probSumH.begin()+HAP_NUMBER, 0.0f);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		fill(prob.begin()+i, prob.begin()+i+HAP_NUMBER, (probSumK[k]*_nt)+_tFreq);
		for (int h = 0 ; h < HAP_NUMBER ; h++) probSumH[h] += prob[i+h];
	}
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************					SUM Ks				************************/
/*******************************************************************************/

inline
void haplotype_segment_double::SUMK() {
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		probSumK[k] = prob[i+0] + prob[i+1] + prob[i+2] + prob[i+3] + prob[i+4] + prob[i+5] + prob[i+6] + prob[i+7];
	}
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS			************************/
/*******************************************************************************/

inline
bool haplotype_segment_double::TRANS_HAP() {
	sumHProbs = 0.0;
	unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
	double fact1 = M.nt[curr_abs_locus-1] / AlphaSumSum[curr_rel_segment_index - 1];
	fill_n(HProbs, HAP_NUMBER*HAP_NUMBER, 0.0f);
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		//__m256 _sum = _mm256_set1_ps(0.0f);
		double fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * M.t[curr_abs_locus - 1] / n_cond_haps;
		for (int k = 0 ; k < n_cond_haps ; k ++) {
			for (int h2 = 0 ; h2 < HAP_NUMBER ; h2++) HProbs[h1*HAP_NUMBER+h2]+=((Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1]*fact1 + fact2)*prob[k*HAP_NUMBER+h2]);
		}
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (isnan(sumHProbs) || isinf(sumHProbs) || sumHProbs < numeric_limits<double>::min());
}

inline
bool haplotype_segment_double::TRANS_DIP_MULT() {
	sumDProbs= 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = (double)(HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)] * scaling) * (double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)] * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (isnan(sumDProbs) || isinf(sumDProbs) || sumDProbs < numeric_limits<double>::min());
}

inline
bool haplotype_segment_double::TRANS_DIP_ADD() {
	sumDProbs = 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = (double)(HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)] * scaling) + (double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)] * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (isnan(sumDProbs) || isinf(sumDProbs) || sumDProbs < numeric_limits<double>::min());
}

inline
void haplotype_segment_double::IMPUTE(vector < float > & missing_probabilities) {
	vector < vector < double > > _sumA = vector < vector < double > > (2, vector < double > (HAP_NUMBER, 0.0f));
	vector < double > _scale = AlphaSumMissing[curr_rel_missing];
	for (int h = 0 ; h < HAP_NUMBER ; h++) _scale[h] = 1.0f / _scale[h];
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		for (int h = 0 ; h < HAP_NUMBER ; h++) _sumA[ah][h] += (AlphaMissing[curr_rel_missing][i+h]*_scale[h])*prob[i+h];
	}
	for (int h = 0 ; h < HAP_NUMBER ; h ++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = _sumA[1][h] / (_sumA[0][h]+_sumA[1][h]);
}

#endif
