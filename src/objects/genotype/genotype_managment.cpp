#include <objects/genotype/genotype_header.h>

genotype::genotype(unsigned int _index) {
	index = _index;
	n_segments = 0;
	n_variants = 0;
	n_ambiguous = 0;
	std::fill(curr_dipcodes, curr_dipcodes + 64, 0);
	this->name = "";
}

genotype::~genotype() {
	free();
}

void genotype::free() {
	std::fill(curr_dipcodes, curr_dipcodes + 64, 0);
	name = "";
	vector < unsigned char > ().swap(Variants);
	vector < unsigned char > ().swap(Ambiguous);
	vector < unsigned long > ().swap(Diplotypes);
	vector < unsigned short > ().swap(Lengths);
}

void genotype::make(vector < unsigned char > & DipSampled) {
	for (unsigned int s = 0, vabs = 0, a = 0 ; s < n_segments ; s ++) {
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			if (VAR_GET_AMB(MOD2(vabs), Variants[DIV2(vabs)])) {
				HAP_GET(Ambiguous[a], hap0)?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
				HAP_GET(Ambiguous[a], hap1)?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
				a++;
			}
		}
	}
}
