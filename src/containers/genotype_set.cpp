#include <containers/genotype_set.h>

genotype_set::genotype_set() {
	n_site = 0;
	n_ind = 0;
}

genotype_set::~genotype_set() {
	for (int i = 0 ; i< vecG.size() ; i ++) delete vecG[i];
	vecG.clear();
	n_site = 0;
	n_ind = 0;
}

void genotype_set::imputeMonomorphic(variant_map & V) {
	for (unsigned int v = 0 ; v < V.size() ; v ++) {
		if (V.vec_pos[v]->isMonomorphic()) {
			bool uallele = (V.vec_pos[v]->cref)?false:true;
			for (unsigned int i = 0 ; i < vecG.size() ; i ++) {
				VAR_SET_HOM(MOD2(v), vecG[i]->Variants[DIV2(v)]);
				uallele?VAR_SET_HAP0(MOD2(v), vecG[i]->Variants[DIV2(v)]):VAR_CLR_HAP0(MOD2(v), vecG[i]->Variants[DIV2(v)]);
				uallele?VAR_SET_HAP1(MOD2(v), vecG[i]->Variants[DIV2(v)]):VAR_CLR_HAP1(MOD2(v), vecG[i]->Variants[DIV2(v)]);
			}
			if (uallele) V.vec_pos[v]->cref = 0;
			else V.vec_pos[v]->calt = 0;
			V.vec_pos[v]->cmis = 0;
		}
	}
}

unsigned int genotype_set::largestNumberOfTransitions() {
	unsigned int maxT = 0;
	for (int i = 0 ; i < n_ind ; i ++) {
		unsigned int nTrans = vecG[i]->n_transitions;
		if (nTrans > maxT) maxT = nTrans;
	}
	return maxT;
}

unsigned long genotype_set::numberOfSegments() {
	unsigned long size = 0;
	for (int i = 0 ; i < n_ind ; i ++) size += vecG[i]->n_segments;
	return size;
}

void genotype_set::masking() {
	for (int i = 0 ; i < n_ind ; i ++) vecG[i]->mask();
}

void genotype_set::solve() {
	tac.clock();
	for (int i = 0 ; i < vecG.size() ; i ++) vecG[i]->solve();
	vrb.bullet("HAP solving (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
