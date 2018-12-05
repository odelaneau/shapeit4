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
#include <io/genotype_reader.h>

genotype_reader::genotype_reader(haplotype_set & _H, genotype_set & _G, variant_map & _V, string _region, bool _use_PS_field) : H(_H), G(_G), V(_V) {
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	region = _region;
	n_geno_het = 0;
	n_geno_hom = 0;
	n_geno_ips = 0;
	n_geno_sca = 0;
	n_geno_mis = 0;
	use_PS_field = _use_PS_field;
}

genotype_reader::~genotype_reader() {
	n_variants = 0;
	n_main_samples = 0;
	n_ref_samples = 0;
	region = "";
}

void genotype_reader::allocateGenotypes() {
	assert(n_variants != 0 && (n_main_samples+n_ref_samples) != 0);
	//Genotypes
	G.vecG = vector < genotype * > (n_main_samples);
	for (unsigned int i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i] = new genotype (i);
		G.vecG[i]->n_variants = n_variants;
		G.vecG[i]->Variants = vector < unsigned char > (DIV2(n_variants) + MOD2(n_variants), 0);
	}
	G.n_ind = n_main_samples;
	G.n_site = n_variants;
	//Haplotypes
	H.n_ind = n_main_samples;
	H.n_hap = 2 * (n_main_samples + n_ref_samples);
	H.n_site = n_variants;
	H.H_opt_var.allocate(H.n_site, H.n_hap);
	H.H_opt_hap.allocate(H.n_hap, H.n_site);
}

void genotype_reader::setPScodes(int * ps_arr, int nps) {
	if (nps != n_main_samples) PScodes.clear();
	else {
		int ps_idx = 0;
		unordered_map < int , int > :: iterator it;
		PScodes = vector < int > (n_main_samples, 0);
		for (int i = 0 ; i < n_main_samples ; i ++) {
			if (ps_arr[i] != bcf_int32_missing) {
				it = PSmap.find(ps_arr[i]);
				if (it == PSmap.end()) {
					ps_idx = PSmap.size()+1;
					PSmap.insert(pair < int, int > (ps_arr[i], ps_idx));
				} else ps_idx = it->second;
				PScodes[i] = ps_idx;
			}
		}
	}
}

void genotype_reader::scanGenotypes(string fmain) {
	vrb.wait("  * VCF/BCF scanning");
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");
	n_variants = 0;
	n_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
	bcf1_t * line;
	while(bcf_sr_next_line (sr)) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) n_variants++;
	}
	bcf_sr_destroy(sr);
	if (n_variants == 0) vrb.error("No variants to be phased in [" + fmain + "]");
	vrb.bullet("VCF/BCF scanning [N=" + stb.str(n_main_samples) + " / L=" + stb.str(n_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void genotype_reader::scanGenotypes(string fmain, string fref) {
	vrb.wait("  * VCF/BCF scanning");
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fmain.c_str()))) vrb.error("Problem opening index file for [" + fmain + "]");
	if(!(bcf_sr_add_reader (sr, fref.c_str()))) vrb.error("Problem opening index file for [" + fref + "]");
	n_variants = 0;
	n_main_samples = bcf_hdr_nsamples(sr->readers[0].header);
	n_ref_samples = bcf_hdr_nsamples(sr->readers[1].header);
	int nset;
	bcf1_t * line_main, * line_ref;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_main =  bcf_sr_get_line(sr, 0);
			line_ref =  bcf_sr_get_line(sr, 1);
			if (line_main->n_allele == 2 && line_ref->n_allele == 2) n_variants ++;
		}
	}
	bcf_sr_destroy(sr);
	if (n_variants == 0) vrb.error("No variants to be phased in files");
	vrb.bullet("VCF/BCF scanning [Nm=" + stb.str(n_main_samples) + " / Nr=" + stb.str(n_ref_samples) + " / L=" + stb.str(n_variants) + " / Reg=" + region + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
