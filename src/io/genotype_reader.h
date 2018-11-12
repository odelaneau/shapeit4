#ifndef _GENOTYPE_READER_H
#define _GENOTYPE_READER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>

class genotype_reader {
public:
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;
	string region;
	//COUNTS
	bool use_PS_field;
	unsigned long n_variants;
	unsigned long n_main_samples;
	unsigned long n_ref_samples;
	unsigned long n_geno_tot;
	unsigned long n_geno_het;
	unsigned long n_geno_hom;
	unsigned long n_geno_ips;
	unsigned long n_geno_sca;
	unsigned long n_geno_mis;
	//PHASESETS
	unordered_map < int, int > PSmap;
	vector < int > PScodes;

	//CONSTRUCTORS/DESCTRUCTORS
	genotype_reader(haplotype_set &, genotype_set &, variant_map &, string regions, bool use_PS_field);
	~genotype_reader();

	//IO
	void scanGenotypes(string funphased);
	void scanGenotypes(string funphased, string fphased);
	void allocateGenotypes();
	void readGenotypes0(string);
	void readGenotypes1(string, string);
	void readGenotypes2(string, string);
	void readGenotypes3(string, string, string);
	void setPScodes(int * ps_arr, int nps);
};

#endif
