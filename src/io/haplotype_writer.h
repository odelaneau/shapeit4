#ifndef _HAPLOTYPE_WRITER_H
#define _HAPLOTYPE_WRITER_H

#include <utils/otools.h>

#include <containers/variant_map.h>
#include <containers/haplotype_set.h>
#include <containers/genotype_set.h>


class haplotype_writer {
public:
	//DATA
	haplotype_set & H;
	genotype_set & G;
	variant_map & V;

	//CONSTRUCTORS/DESCTRUCTORS
	haplotype_writer(haplotype_set &, genotype_set &, variant_map &);
	~haplotype_writer();

	//IO
	void writeHaplotypes(string foutput);
};

#endif
