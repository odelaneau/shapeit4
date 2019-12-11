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
#ifndef _HAPLOTYPE_SEGMENT_H
#define _HAPLOTYPE_SEGMENT_H

#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>

#ifdef __AVX2__
	#include <immintrin.h>
#endif

#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator < T, 32 > >;

class haplotype_segment {
private:
	//EXTERNAL DATA
	bitmatrix & H;
	vector < unsigned int > & idxH;
	hmm_parameters & M;
	genotype * G;

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
	unsigned int n_cond_haps;
	unsigned int n_missing;


	//CURSORS
	int curr_segment_index;
	int curr_segment_locus;
	int curr_abs_locus;
	int curr_rel_locus;
	int curr_rel_segment_index;
	int curr_abs_ambiguous;
	int curr_abs_transition;
	int curr_abs_missing;
	int curr_rel_missing;

	//DYNAMIC ARRAYS
	float probSumT1;
	float probSumT2;
	aligned_vector32 < float > prob1;
	aligned_vector32 < float > prob2;

	aligned_vector32 < float > probSumK1;
	aligned_vector32 < float > probSumK2;
	aligned_vector32 < float > probSumH1;
	aligned_vector32 < float > probSumH2;
	vector < aligned_vector32 < float > > Alpha;
	vector < aligned_vector32 < float > > Beta;
	vector < aligned_vector32 < float > > AlphaSum;
	aligned_vector32 < float > AlphaSumSum;
	aligned_vector32 < float > BetaSum;

	//IMPUTED MISSING DATA
	vector < aligned_vector32 < float > > AlphaMissing;
	vector < aligned_vector32 < float > > AlphaSumMissing;
	aligned_vector32 < float > ProbM0sums, ProbM1sums, ProbM;

	//STATIC ARRAYS
	float sumHProbs;
	float HProbs [HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
	double sumDProbs;
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));


	//INLINED AND UNROLLED ROUTINES
	void HOM(bool);
	void AMB(bool);
	void MIS(bool);

	void SUM(bool);
	void SUMK(bool);
	void COLLAPSE(bool, bool);
	void RUN(bool, bool);
	bool TRANSH();
	bool TRANSD(int &);

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment(genotype *, bitmatrix &, vector < unsigned int > &, coordinates &, hmm_parameters &);
	~haplotype_segment();

	void forward();
	void backward(vector < float > &);
	int expectation(vector < double > &, vector < float > &);
};

inline
void haplotype_segment::HOM(bool paired) {
	bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
	if (paired) {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = H.get(idxH[k], curr_abs_locus);
			if (ag != ah) fill(prob2.begin() + i, prob2.begin() + i + HAP_NUMBER, M.ed);
			else fill(prob2.begin() + i, prob2.begin() + i + HAP_NUMBER, M.ee);
		}
	} else {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = H.get(idxH[k], curr_abs_locus);
			if (ag != ah) fill(prob1.begin() + i, prob1.begin() + i + HAP_NUMBER, M.ed);
			else fill(prob1.begin() + i, prob1.begin() + i + HAP_NUMBER, M.ee);
		}
	}
}

inline
void haplotype_segment::AMB(bool paired) {
	float galleles0[HAP_NUMBER], galleles1[HAP_NUMBER];
	galleles0[0] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],0)?M.ed:M.ee;
	galleles0[1] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],1)?M.ed:M.ee;
	galleles0[2] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],2)?M.ed:M.ee;
	galleles0[3] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],3)?M.ed:M.ee;
	galleles0[4] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],4)?M.ed:M.ee;
	galleles0[5] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],5)?M.ed:M.ee;
	galleles0[6] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],6)?M.ed:M.ee;
	galleles0[7] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],7)?M.ed:M.ee;
	galleles1[0] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],0)?M.ee:M.ed;
	galleles1[1] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],1)?M.ee:M.ed;
	galleles1[2] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],2)?M.ee:M.ed;
	galleles1[3] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],3)?M.ee:M.ed;
	galleles1[4] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],4)?M.ee:M.ed;
	galleles1[5] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],5)?M.ee:M.ed;
	galleles1[6] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],6)?M.ee:M.ed;
	galleles1[7] = HAP_GET(G->Ambiguous[curr_abs_ambiguous],7)?M.ee:M.ed;
	if (paired) {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool a = H.get(idxH[k], curr_abs_locus);
			if (a) memcpy(&prob2[i], &galleles1[0], HAP_NUMBER*sizeof(float));
			else memcpy(&prob2[i], &galleles0[0], HAP_NUMBER*sizeof(float));
		}
	} else {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool a = H.get(idxH[k], curr_abs_locus);
			if (a) memcpy(&prob1[i], &galleles1[0], HAP_NUMBER*sizeof(float));
			else memcpy(&prob1[i], &galleles0[0], HAP_NUMBER*sizeof(float));
		}
	}
}

inline
void haplotype_segment::MIS(bool paired) {
	if (paired) fill(prob2.begin(), prob2.end(), 1.0);
	else fill(prob1.begin(), prob1.end(), 1.0);
}

inline
void haplotype_segment::SUMK(bool paired) {
	if (paired) {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			probSumK2[k] = prob2[i+0] + prob2[i+1] + prob2[i+2] + prob2[i+3] + prob2[i+4] + prob2[i+5] + prob2[i+6] + prob2[i+7];
		}
	} else {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			probSumK1[k] = prob1[i+0] + prob1[i+1] + prob1[i+2] + prob1[i+3] + prob1[i+4] + prob1[i+5] + prob1[i+6] + prob1[i+7];
		}
	}
}

#ifdef __AVX2__

inline
void haplotype_segment::RUN(bool forward, bool paired) {
	if (paired) {
		float nt = M.nt[curr_abs_locus-forward] / probSumT1;
		float tfreq = M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT1);
		__m256 _tFreq = _mm256_set_ps(probSumH1[7] * tfreq, probSumH1[6] * tfreq, probSumH1[5] * tfreq, probSumH1[4] * tfreq, probSumH1[3] * tfreq, probSumH1[2] * tfreq, probSumH1[1] * tfreq, probSumH1[0] * tfreq);
		__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT1);
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			__m256 _prob1 = _mm256_load_ps(&prob1[i]);
			__m256 _prob2 = _mm256_load_ps(&prob2[i]);
			_prob1 = _mm256_fmadd_ps(_prob1, _nt, _tFreq);
			_prob2 = _mm256_mul_ps(_prob1, _prob2);
			_mm256_store_ps(&prob2[i], _prob2);
		}
	} else {
		float nt = M.nt[curr_abs_locus-forward] / probSumT2;
		float tfreq = M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT2);
		__m256 _tFreq = _mm256_set_ps(probSumH2[7] * tfreq, probSumH2[6] * tfreq, probSumH2[5] * tfreq, probSumH2[4] * tfreq, probSumH2[3] * tfreq, probSumH2[2] * tfreq, probSumH2[1] * tfreq, probSumH2[0] * tfreq);
		__m256 _nt = _mm256_set1_ps(M.nt[curr_abs_locus-forward] / probSumT2);
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			__m256 _prob1 = _mm256_load_ps(&prob1[i]);
			__m256 _prob2 = _mm256_load_ps(&prob2[i]);
			_prob2 = _mm256_fmadd_ps(_prob2, _nt, _tFreq);
			_prob1 = _mm256_mul_ps(_prob1, _prob2);
			_mm256_store_ps(&prob1[i], _prob1);
		}
	}
}

inline
void haplotype_segment::SUM(bool paired) {
	__m256 _sum = _mm256_set1_ps(0.0f);
	if (paired) {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			__m256 _prob2 = _mm256_load_ps(&prob2[i]);
			_sum = _mm256_add_ps(_sum, _prob2);
		}
		_mm256_store_ps(&probSumH2[0], _sum);
		probSumT2 = probSumH2[0] + probSumH2[1] + probSumH2[2] + probSumH2[3] + probSumH2[4] + probSumH2[5] + probSumH2[6] + probSumH2[7];
	} else {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			__m256 _prob1 = _mm256_load_ps(&prob1[i]);
			_sum = _mm256_add_ps(_sum, _prob1);
		}
		_mm256_store_ps(&probSumH1[0], _sum);
		probSumT1 = probSumH1[0] + probSumH1[1] + probSumH1[2] + probSumH1[3] + probSumH1[4] + probSumH1[5] + probSumH1[6] + probSumH1[7];
	}
}

inline
void haplotype_segment::COLLAPSE(bool forward, bool paired) {
	if (paired) {
		float tmp_prob0 = M.nt[curr_abs_locus-forward] / probSumT1;
		float tmp_prob1 = M.t[curr_abs_locus-forward] / n_cond_haps;
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			__m256 _factor = _mm256_set1_ps(probSumK1[k] * tmp_prob0 + tmp_prob1);
			__m256 _prob2 = _mm256_load_ps(&prob2[i]);
			_prob2 = _mm256_mul_ps(_factor, _prob2);
			_mm256_store_ps(&prob2[i], _prob2);
		}
	} else {
		float tmp_prob0 = M.nt[curr_abs_locus-forward] / probSumT2;
		float tmp_prob1 = M.t[curr_abs_locus-forward] / n_cond_haps;
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			__m256 _factor = _mm256_set1_ps(probSumK2[k] * tmp_prob0 + tmp_prob1);
			__m256 _prob1 = _mm256_load_ps(&prob1[i]);
			_prob1 = _mm256_mul_ps(_factor, _prob1);
			_mm256_store_ps(&prob1[i], _prob1);
		}
	}
}


inline
bool haplotype_segment::TRANSH() {
	sumHProbs = 0.0;
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		__m256 _sum = _mm256_set1_ps(0.0f);

		float fact1 = M.nt[curr_abs_locus-1] / AlphaSumSum[curr_rel_segment_index - 1];
		float fact2 = (AlphaSum[curr_rel_segment_index - 1][h1]/AlphaSumSum[curr_rel_segment_index - 1]) * M.t[curr_abs_locus - 1] / n_cond_haps;

		for (int k = 0 ; k < n_cond_haps ; k ++) {
			float alpha = Alpha[curr_rel_segment_index - 1][k*HAP_NUMBER + h1] * fact1 + fact2;
			__m256 _alpha = _mm256_set1_ps(alpha);
			__m256 _beta = _mm256_load_ps(&Beta[curr_rel_segment_index][k*HAP_NUMBER]);
			_sum = _mm256_add_ps(_sum, _mm256_mul_ps(_alpha, _beta));
		}
		_mm256_store_ps(&HProbs[h1*HAP_NUMBER], _sum);
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (isnan(sumHProbs) || sumHProbs < numeric_limits<float>::min());
}

#else

inline
void haplotype_segment::SUM(bool paired) {
	float sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;
	if (paired) {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			sum0 += prob2[i + 0];
			sum1 += prob2[i + 1];
			sum2 += prob2[i + 2];
			sum3 += prob2[i + 3];
			sum4 += prob2[i + 4];
			sum5 += prob2[i + 5];
			sum6 += prob2[i + 6];
			sum7 += prob2[i + 7];
		}
		probSumH2[0] = sum0;
		probSumH2[1] = sum1;
		probSumH2[2] = sum2;
		probSumH2[3] = sum3;
		probSumH2[4] = sum4;
		probSumH2[5] = sum5;
		probSumH2[6] = sum6;
		probSumH2[7] = sum7;
		probSumT2 = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
	} else {
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			sum0 += prob1[i + 0];
			sum1 += prob1[i + 1];
			sum2 += prob1[i + 2];
			sum3 += prob1[i + 3];
			sum4 += prob1[i + 4];
			sum5 += prob1[i + 5];
			sum6 += prob1[i + 6];
			sum7 += prob1[i + 7];
		}
		probSumH1[0] = sum0;
		probSumH1[1] = sum1;
		probSumH1[2] = sum2;
		probSumH1[3] = sum3;
		probSumH1[4] = sum4;
		probSumH1[5] = sum5;
		probSumH1[6] = sum6;
		probSumH1[7] = sum7;
		probSumT1 = sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
	}
}

inline
void haplotype_segment::COLLAPSE(bool forward, bool paired) {
	if (paired) {
		float tmp_prob0 = M.nt[curr_abs_locus-forward] / probSumT1;
		float tmp_prob1 = M.t[curr_abs_locus-forward] / n_cond_haps;
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			float factor = probSumK1[k] * tmp_prob0 + tmp_prob1;
			prob2[i + 0] *= factor;
			prob2[i + 1] *= factor;
			prob2[i + 2] *= factor;
			prob2[i + 3] *= factor;
			prob2[i + 4] *= factor;
			prob2[i + 5] *= factor;
			prob2[i + 6] *= factor;
			prob2[i + 7] *= factor;
		}
	} else {
		float tmp_prob0 = M.nt[curr_abs_locus-forward] / probSumT2;
		float tmp_prob1 = M.t[curr_abs_locus-forward] / n_cond_haps;
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			float factor = probSumK2[k] * tmp_prob0 + tmp_prob1;
			prob1[i + 0] *= factor;
			prob1[i + 1] *= factor;
			prob1[i + 2] *= factor;
			prob1[i + 3] *= factor;
			prob1[i + 4] *= factor;
			prob1[i + 5] *= factor;
			prob1[i + 6] *= factor;
			prob1[i + 7] *= factor;
		}
	}
}

inline
void haplotype_segment::RUN(bool forward, bool paired) {
	if (paired) {
		float nt = M.nt[curr_abs_locus-forward] / probSumT1;
		float tfreq = M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT1);
		float tFreq0 = probSumH1[0] * tfreq;
		float tFreq1 = probSumH1[1] * tfreq;
		float tFreq2 = probSumH1[2] * tfreq;
		float tFreq3 = probSumH1[3] * tfreq;
		float tFreq4 = probSumH1[4] * tfreq;
		float tFreq5 = probSumH1[5] * tfreq;
		float tFreq6 = probSumH1[6] * tfreq;
		float tFreq7 = probSumH1[7] * tfreq;
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			prob2[i + 0] *= prob1[i + 0] * nt + tFreq0;
			prob2[i + 1] *= prob1[i + 1] * nt + tFreq1;
			prob2[i + 2] *= prob1[i + 2] * nt + tFreq2;
			prob2[i + 3] *= prob1[i + 3] * nt + tFreq3;
			prob2[i + 4] *= prob1[i + 4] * nt + tFreq4;
			prob2[i + 5] *= prob1[i + 5] * nt + tFreq5;
			prob2[i + 6] *= prob1[i + 6] * nt + tFreq6;
			prob2[i + 7] *= prob1[i + 7] * nt + tFreq7;
		}
	} else {
		float nt = M.nt[curr_abs_locus-forward] / probSumT2;
		float tfreq = M.t[curr_abs_locus-forward] / (n_cond_haps * probSumT2);
		float tFreq0 = probSumH2[0] * tfreq;
		float tFreq1 = probSumH2[1] * tfreq;
		float tFreq2 = probSumH2[2] * tfreq;
		float tFreq3 = probSumH2[3] * tfreq;
		float tFreq4 = probSumH2[4] * tfreq;
		float tFreq5 = probSumH2[5] * tfreq;
		float tFreq6 = probSumH2[6] * tfreq;
		float tFreq7 = probSumH2[7] * tfreq;
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			prob1[i + 0] *= prob2[i + 0] * nt + tFreq0;
			prob1[i + 1] *= prob2[i + 1] * nt + tFreq1;
			prob1[i + 2] *= prob2[i + 2] * nt + tFreq2;
			prob1[i + 3] *= prob2[i + 3] * nt + tFreq3;
			prob1[i + 4] *= prob2[i + 4] * nt + tFreq4;
			prob1[i + 5] *= prob2[i + 5] * nt + tFreq5;
			prob1[i + 6] *= prob2[i + 6] * nt + tFreq6;
			prob1[i + 7] *= prob2[i + 7] * nt + tFreq7;
		}
	}
}

inline
bool haplotype_segment::TRANSH() {
	sumHProbs = 0.0;
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		float sum0 = 0.0, sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0;

		float fact1 = M.nt[curr_abs_locus-1] / AlphaSumSum[curr_rel_segment_index - 1];
		float fact2 = (AlphaSum[curr_rel_segment_index - 1][h1]/AlphaSumSum[curr_rel_segment_index - 1]) * M.t[curr_abs_locus - 1] / n_cond_haps;

		for (int k = 0 ; k < n_cond_haps ; k ++) {
			//float alpha = Alpha[curr_rel_segment_index - 1][k*HAP_NUMBER + h1] * M.nt[curr_abs_locus-1] + AlphaSum[curr_rel_segment_index - 1][h1] * M.tfreq[curr_abs_locus - 1];
			//float alpha = Alpha[curr_rel_segment_index - 1][k*HAP_NUMBER + h1] * M.nt[curr_abs_locus-1] + AlphaSum[curr_rel_segment_index - 1][h1] * M.t[curr_abs_locus - 1] / n_cond_haps;
			float alpha = Alpha[curr_rel_segment_index - 1][k*HAP_NUMBER + h1] * fact1 + fact2;
			sum0 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 0];
			sum1 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 1];
			sum2 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 2];
			sum3 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 3];
			sum4 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 4];
			sum5 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 5];
			sum6 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 6];
			sum7 += alpha * Beta[curr_rel_segment_index][k*HAP_NUMBER + 7];
		}
		HProbs[h1*HAP_NUMBER+0] = sum0;
		HProbs[h1*HAP_NUMBER+1] = sum1;
		HProbs[h1*HAP_NUMBER+2] = sum2;
		HProbs[h1*HAP_NUMBER+3] = sum3;
		HProbs[h1*HAP_NUMBER+4] = sum4;
		HProbs[h1*HAP_NUMBER+5] = sum5;
		HProbs[h1*HAP_NUMBER+6] = sum6;
		HProbs[h1*HAP_NUMBER+7] = sum7;
		sumHProbs += sum0 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6 + sum7;
	}
	return (isnan(sumHProbs) || sumHProbs < numeric_limits<float>::min());
}

#endif


inline
bool haplotype_segment::TRANSD(int & n_underflows_recovered) {
	sumDProbs= 0.0;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					unsigned int prev_hap0 = DIP_HAP0(pd);
					unsigned int prev_hap1 = DIP_HAP1(pd);
					unsigned int next_hap0 = DIP_HAP0(nd);
					unsigned int next_hap1 = DIP_HAP1(nd);
					DProbs[t] = (double)(HProbs[prev_hap0*HAP_NUMBER+next_hap0] * scaling) * (double)(HProbs[prev_hap1*HAP_NUMBER+next_hap1] * scaling);;
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	if (sumDProbs < numeric_limits<double>::min()) {
		sumDProbs = 0.0;
		n_underflows_recovered++;
		for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
			if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
				for (int nd = 0 ; nd < 64 ; ++nd) {
					if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
						unsigned int prev_hap0 = DIP_HAP0(pd);
						unsigned int prev_hap1 = DIP_HAP1(pd);
						unsigned int next_hap0 = DIP_HAP0(nd);
						unsigned int next_hap1 = DIP_HAP1(nd);
						DProbs[t] = (double)(HProbs[prev_hap0*HAP_NUMBER+next_hap0] * scaling) + (double)(HProbs[prev_hap1*HAP_NUMBER+next_hap1] * scaling);;
						sumDProbs += DProbs[t];
						t++;
					}
				}
			}
		}
	}
	return (isnan(sumDProbs) || sumDProbs < numeric_limits<double>::min());
}

#endif
