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

//**********************************************************************************//
//								ONE VCF/BCF PROCESSED								//
//								1. main genotype data								//
//**********************************************************************************//
void genotype_reader::readGenotypes0(string funphased) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader(sr, funphased.c_str());
	for (int i = 0 ; i < n_main_samples ; i ++) G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
	bcf1_t * line;
	int ngt_main, *gt_arr_main = NULL, ngt_arr_main = 0;
	int nps_main, *ps_arr_main = NULL, nps_arr_main = 0;
	unsigned int i_variant = 0;
	while(bcf_sr_next_line (sr)) {
		line =  bcf_sr_get_line(sr, 0);
		if (line->n_allele == 2) {
			bcf_unpack(line, BCF_UN_STR);
			string chr = bcf_hdr_id2name(sr->readers[0].header, line->rid);
			unsigned int pos = line->pos + 1;
			string id = string(line->d.id);
			string ref = string(line->d.allele[0]);
			string alt = string(line->d.allele[1]);
			variant * newV = new variant (chr, pos, id, ref, alt, V.size());
			unsigned int cref = 0, calt = 0, cmis = 0;
			ngt_main = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr_main, &ngt_arr_main);
			assert(ngt_main == 2 * n_main_samples);
			if (use_PS_field) {
				nps_main = bcf_get_format_int32(sr->readers[0].header, line, "PS", &ps_arr_main, &nps_arr_main);
				setPScodes(ps_arr_main, nps_main);
			}
			for(int i = 0 ; i < 2 * n_main_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_main[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_main[i+1])==1);
				bool mi = (gt_arr_main[i+0] == bcf_gt_missing || gt_arr_main[i+1] == bcf_gt_missing);
				bool he = !mi && a0 != a1;
				bool ho = !mi && a0 == a1;
				//bool ps = (mi || he) && use_PS_field;
				bool ps = he && use_PS_field;
				bool ph = (bcf_gt_is_phased(gt_arr_main[i+0]) || bcf_gt_is_phased(gt_arr_main[i+1])) && he && PScodes.size() > 0;
				if (a0) VAR_SET_HAP0(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (a1) VAR_SET_HAP1(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (mi) VAR_SET_MIS(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (he) VAR_SET_HET(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (ps) G.vecG[DIV2(i)]->pushPS(a0, a1, ph?PScodes[i/2]:0);
				if (!mi) { a0?calt++:cref++; a1?calt++:cref++; }
				else cmis ++;
				n_geno_het += he;
				n_geno_hom += ho;
				n_geno_mis += mi;
				n_geno_ips += ph;
			}
			newV->cref = cref;newV->calt = calt;newV->cmis = cmis;
			V.push(newV);
			i_variant ++;
			vrb.progress("  * VCF/BCF parsing", i_variant*1.0/n_variants);
		}
	}
	free(gt_arr_main);
	if (ps_arr_main) free(ps_arr_main);
	bcf_sr_destroy(sr);
	// Report
	n_geno_tot = n_main_samples*n_variants;
	string str0 = "Hom=" + stb.str(n_geno_hom*100.0/n_geno_tot, 1) + "%";
	string str1 = "Het=" + stb.str(n_geno_het*100.0/n_geno_tot, 1) + "%" + (use_PS_field?(" / Pha=" + stb.str(n_geno_ips*100.0/n_geno_tot, 3) + "%"):(""));
	string str2 = "Mis=" + stb.str(n_geno_mis*100.0/n_geno_tot, 1) + "%";
	string str3 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+"] ("+str3+")");
}

//**********************************************************************************//
//								TWO VCF/BCF PROCESSED								//
//								1. main genotype data								//
//								2. reference haplotype data							//
//**********************************************************************************//
void genotype_reader::readGenotypes1(string funphased, string freference) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, funphased.c_str());
	bcf_sr_add_reader (sr, freference.c_str());
	for (int i = 0 ; i < n_main_samples ; i ++) G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
	unsigned int i_variant = 0, nset = 0, n_ref_missing = 0, n_ref_unphased = 0;
	int ngt_main, *gt_arr_main = NULL, ngt_arr_main = 0;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	int nps_main, *ps_arr_main = NULL, nps_arr_main = 0;
	bcf1_t * line_main, * line_ref;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_main =  bcf_sr_get_line(sr, 0);
			line_ref =  bcf_sr_get_line(sr, 1);
			if (line_main->n_allele == 2 && line_ref->n_allele == 2) {
				bcf_unpack(line_main, BCF_UN_STR);
				string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
				int pos = line_main->pos + 1;
				string id = string(line_main->d.id);
				string ref = string(line_main->d.allele[0]);
				string alt = string(line_main->d.allele[1]);
				variant * newV = new variant (chr, pos, id, ref, alt, V.size());
				unsigned int cref = 0, calt = 0, cmis = 0;
				ngt_main = bcf_get_genotypes(sr->readers[0].header, line_main, &gt_arr_main, &ngt_arr_main); assert(ngt_main == 2 * n_main_samples);
				ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref); assert(ngt_ref == 2 * n_ref_samples);
				if (use_PS_field) {
					nps_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PS", &ps_arr_main, &nps_arr_main);
					setPScodes(ps_arr_main, nps_main);
				}
				for(int i = 0 ; i < 2 * n_main_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(gt_arr_main[i+0])==1);
					bool a1 = (bcf_gt_allele(gt_arr_main[i+1])==1);
					bool mi = (gt_arr_main[i+0] == bcf_gt_missing || gt_arr_main[i+1] == bcf_gt_missing);
					bool he = !mi && a0 != a1;
					bool ho = !mi && a0 == a1;
					//bool ps = (mi || he) && use_PS_field;
					bool ps = he && use_PS_field;
					bool ph = (bcf_gt_is_phased(gt_arr_main[i+0]) || bcf_gt_is_phased(gt_arr_main[i+1])) && he && PScodes.size() > 0;
					if (a0) VAR_SET_HAP0(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
					if (a1) VAR_SET_HAP1(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
					if (mi) VAR_SET_MIS(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
					if (he) VAR_SET_HET(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
					if (ps) G.vecG[DIV2(i)]->pushPS(a0, a1, ph?PScodes[i/2]:0);
					if (!mi) { a0?calt++:cref++; a1?calt++:cref++; }
					else cmis ++;
					n_geno_het += he;
					n_geno_hom += ho;
					n_geno_mis += mi;
					n_geno_ips += ph;
				}
				for(int i = 0 ; i < 2 * n_ref_samples ; i += 2) {
					bool a0 = (bcf_gt_allele(gt_arr_ref[i+0])==1);
					bool a1 = (bcf_gt_allele(gt_arr_ref[i+1])==1);
					n_ref_missing += (gt_arr_ref[i+0] == bcf_gt_missing || gt_arr_ref[i+1] == bcf_gt_missing);
					n_ref_unphased += !bcf_gt_is_phased(gt_arr_ref[i+1]);
					H.H_opt_hap.set(i+2*n_main_samples+0, i_variant, a0);
					H.H_opt_hap.set(i+2*n_main_samples+1, i_variant, a1);
					a0?calt++:cref++;
					a1?calt++:cref++;
				}
				newV->cref = cref;newV->calt = calt;newV->cmis = cmis;
				i_variant ++;
				V.push(newV);
				vrb.progress("  * VCF/BCF parsing", i_variant*1.0/n_variants);
			}
		}
	}
	free(gt_arr_main);
	free(gt_arr_ref);
	if (ps_arr_main) free(ps_arr_main);
	bcf_sr_destroy(sr);
	// Report
	n_geno_tot = n_main_samples*n_variants;
	string str0 = "Hom=" + stb.str(n_geno_hom*100.0/n_geno_tot, 1) + "%";
	string str1 = "Het=" + stb.str(n_geno_het*100.0/n_geno_tot, 1) + "%" + (use_PS_field?(" / Pha=" + stb.str(n_geno_ips*100.0/n_geno_tot, 3) + "%"):(""));
	string str2 = "Mis=" + stb.str(n_geno_mis*100.0/n_geno_tot, 1) + "%";
	string str3 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+"] ("+str3+")");
	if (n_ref_missing > 0) vrb.warning(stb.str(n_ref_missing) + " missing genotypes in the reference panel (randomly imputed)");
	if (n_ref_unphased > 0) vrb.warning(stb.str(n_ref_unphased) + " unphased genotypes in the reference panel (randomly phased)");
}

//**********************************************************************************//
//								TWO VCF/BCF PROCESSED								//
//								1. main genotype data								//
//								2. scaffold haplotype data							//
//**********************************************************************************//
void genotype_reader::readGenotypes2(string funphased, string fphased) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, funphased.c_str());
	if (!bcf_sr_add_reader (sr, fphased.c_str())) vrb.error("Problem opening index file for [" + fphased + "]");

	// Mapping scaffolded samples
	map < string, int > map_names;
	for (int i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
		map_names.insert(pair < string, int > (G.vecG[i]->name, i));
	}
	int n_scaf_samples = bcf_hdr_nsamples(sr->readers[1].header);
	vector < int > mappingS2G = vector < int > (n_scaf_samples, -1);
	for (int i = 0 ; i < n_scaf_samples ; i ++) {
		string scaf_name = string(sr->readers[1].header->samples[i]);
		map < string, int > :: iterator it = map_names.find(scaf_name);
		if (it != map_names.end()) mappingS2G[i] = it->second;
	}

	unsigned int i_variant = 0, nset = 0;
	int ngt_main, *gt_arr_main = NULL, ngt_arr_main = 0;
	int ngt_scaf, *gt_arr_scaf = NULL, ngt_arr_scaf = 0;
	int nps_main, *ps_arr_main = NULL, nps_arr_main = 0;
	bcf1_t * line_main, * line_scaf;
	while ((nset = bcf_sr_next_line (sr))) {
		if ((line_main=bcf_sr_get_line(sr, 0))&&(line_main->n_allele == 2)) {
			bcf_unpack(line_main, BCF_UN_STR);
			string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
			int pos = line_main->pos + 1;
			string id = string(line_main->d.id);
			string ref = string(line_main->d.allele[0]);
			string alt = string(line_main->d.allele[1]);
			variant * newV = new variant (chr, pos, id, ref, alt, V.size());
			unsigned int cref = 0, calt = 0, cmis = 0;
			ngt_main = bcf_get_genotypes(sr->readers[0].header, line_main, &gt_arr_main, &ngt_arr_main); assert(ngt_main == 2 * n_main_samples);
			if (use_PS_field) {
				nps_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PS", &ps_arr_main, &nps_arr_main);
				setPScodes(ps_arr_main, nps_main);
			}
			for(int i = 0 ; i < 2 * n_main_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_main[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_main[i+1])==1);
				bool mi = (gt_arr_main[i+0] == bcf_gt_missing || gt_arr_main[i+1] == bcf_gt_missing);
				bool he = !mi && a0 != a1;
				bool ho = !mi && a0 == a1;
				//bool ps = (mi || he) && use_PS_field;
				bool ps = he && use_PS_field;
				bool ph = (bcf_gt_is_phased(gt_arr_main[i+0]) || bcf_gt_is_phased(gt_arr_main[i+1])) && he && PScodes.size() > 0;
				if (a0) VAR_SET_HAP0(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (a1) VAR_SET_HAP1(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (mi) VAR_SET_MIS(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (he) VAR_SET_HET(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (ps) G.vecG[DIV2(i)]->pushPS(a0, a1, ph?PScodes[i/2]:0);
				if (!mi) { a0?calt++:cref++; a1?calt++:cref++; }
				else cmis ++;
				n_geno_het += he;
				n_geno_hom += ho;
				n_geno_mis += mi;
				n_geno_ips += ph;
			}
			if (line_scaf=bcf_sr_get_line(sr, 1)) {
				ngt_scaf = bcf_get_genotypes(sr->readers[1].header, line_scaf, &gt_arr_scaf, &ngt_arr_scaf); assert(ngt_scaf == 2 * n_scaf_samples);
				for(int i = 0 ; i < 2 * n_scaf_samples ; i += 2) {
					int ind = mappingS2G[DIV2(i)];
					if (ind>=0) {
						bool s0 = (bcf_gt_allele(gt_arr_scaf[i+0])==1);
						bool s1 = (bcf_gt_allele(gt_arr_scaf[i+1])==1);
						bool he = (s0 != s1);
						bool ph = (bcf_gt_is_phased(gt_arr_scaf[i+0]) || bcf_gt_is_phased(gt_arr_scaf[i+1]));
						bool mi = (gt_arr_scaf[i+0] == bcf_gt_missing || gt_arr_scaf[i+1] == bcf_gt_missing);
						if (he && !mi && ph) {
							bool a0 = VAR_GET_HAP0(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
							bool a1 = VAR_GET_HAP1(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
							if (a0!=a1) {
								VAR_SET_SCA(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
								s0?VAR_SET_HAP0(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]):VAR_CLR_HAP0(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
								s1?VAR_SET_HAP1(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]):VAR_CLR_HAP1(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
								n_geno_sca ++;
							}
						}
					}
				}
			}
			newV->cref = cref;newV->calt = calt;newV->cmis = cmis;
			i_variant ++;
			V.push(newV);
			vrb.progress("  * VCF/BCF parsing", i_variant*1.0/n_variants);
		}
	}
	free(gt_arr_main);
	free(gt_arr_scaf);
	if (ps_arr_main) free(ps_arr_main);
	bcf_sr_destroy(sr);
	// Report
	n_geno_tot = n_main_samples*n_variants;
	string str0 = "Hom=" + stb.str(n_geno_hom*100.0/n_geno_tot, 1) + "%";
	string str1 = "Het=" + stb.str(n_geno_het*100.0/n_geno_tot, 1) + "%" + (use_PS_field?(" / Pha=" + stb.str(n_geno_ips*100.0/n_geno_tot, 3) + "%"):(""));
	string str2 = "Sca=" + stb.str(n_geno_sca*100.0/n_geno_tot, 3) + "%";
	string str3 = "Mis=" + stb.str(n_geno_mis*100.0/(n_main_samples*n_variants), 1) + "%";
	string str4 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+" / "+str3+"] ("+str4+")");
}

//**********************************************************************************//
//								THREE VCF/BCF PROCESSED								//
//								1. main genotype data								//
//								2. reference haplotype data							//
//								3. scaffold haplotype data							//
//**********************************************************************************//
void genotype_reader::readGenotypes3(string funphased, string freference, string fphased) {
	tac.clock();
	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	bcf_sr_set_regions(sr, region.c_str(), 0);
	bcf_sr_add_reader (sr, funphased.c_str());
	bcf_sr_add_reader (sr, freference.c_str());
	if (!bcf_sr_add_reader (sr, fphased.c_str())) vrb.error("Problem opening index file for [" + fphased + "]");

	// Mapping scaffolded samples
	map < string, int > map_names;
	for (int i = 0 ; i < n_main_samples ; i ++) {
		G.vecG[i]->name = string(sr->readers[0].header->samples[i]);
		map_names.insert(pair < string, int > (G.vecG[i]->name, i));
	}
	int n_scaf_samples = bcf_hdr_nsamples(sr->readers[2].header);
	vector < int > mappingS2G = vector < int > (n_scaf_samples, -1);
	for (int i = 0 ; i < n_scaf_samples ; i ++) {
		string scaf_name = string(sr->readers[2].header->samples[i]);
		map < string, int > :: iterator it = map_names.find(scaf_name);
		if (it != map_names.end()) mappingS2G[i] = it->second;
	}

	unsigned int i_variant = 0, nset = 0, n_ref_missing = 0, n_ref_unphased = 0;
	int ngt_main, *gt_arr_main = NULL, ngt_arr_main = 0;
	int ngt_scaf, *gt_arr_scaf = NULL, ngt_arr_scaf = 0;
	int nps_main, *ps_arr_main = NULL, nps_arr_main = 0;
	int ngt_ref, *gt_arr_ref = NULL, ngt_arr_ref = 0;
	bcf1_t * line_main, * line_scaf, * line_ref;
	while ((nset = bcf_sr_next_line (sr))) {
		if ((line_main=bcf_sr_get_line(sr, 0))&&(line_ref=bcf_sr_get_line(sr, 1))&&(line_main->n_allele == 2)) {
			bcf_unpack(line_main, BCF_UN_STR);
			string chr = bcf_hdr_id2name(sr->readers[0].header, line_main->rid);
			int pos = line_main->pos + 1;
			string id = string(line_main->d.id);
			string ref = string(line_main->d.allele[0]);
			string alt = string(line_main->d.allele[1]);
			variant * newV = new variant (chr, pos, id, ref, alt, V.size());
			unsigned int cref = 0, calt = 0, cmis = 0;
			ngt_main = bcf_get_genotypes(sr->readers[0].header, line_main, &gt_arr_main, &ngt_arr_main); assert(ngt_main == 2 * n_main_samples);
			ngt_ref = bcf_get_genotypes(sr->readers[1].header, line_ref, &gt_arr_ref, &ngt_arr_ref); assert(ngt_ref == 2 * n_ref_samples);
			if (use_PS_field) {
				nps_main = bcf_get_format_int32(sr->readers[0].header, line_main, "PS", &ps_arr_main, &nps_arr_main);
				setPScodes(ps_arr_main, nps_main);
			}
			for(int i = 0 ; i < 2 * n_main_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_main[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_main[i+1])==1);
				bool mi = (gt_arr_main[i+0] == bcf_gt_missing || gt_arr_main[i+1] == bcf_gt_missing);
				bool he = !mi && a0 != a1;
				bool ho = !mi && a0 == a1;
				//bool ps = (mi || he) && use_PS_field;
				bool ps = he && use_PS_field;
				bool ph = (bcf_gt_is_phased(gt_arr_main[i+0]) || bcf_gt_is_phased(gt_arr_main[i+1])) && he && PScodes.size() > 0;
				if (a0) VAR_SET_HAP0(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (a1) VAR_SET_HAP1(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (mi) VAR_SET_MIS(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (he) VAR_SET_HET(MOD2(i_variant), G.vecG[DIV2(i)]->Variants[DIV2(i_variant)]);
				if (ps) G.vecG[DIV2(i)]->pushPS(a0, a1, ph?PScodes[i/2]:0);
				if (!mi) { a0?calt++:cref++; a1?calt++:cref++; }
				else cmis ++;
				n_geno_het += he;
				n_geno_hom += ho;
				n_geno_mis += mi;
				n_geno_ips += ph;
			}
			for(int i = 0 ; i < 2 * n_ref_samples ; i += 2) {
				bool a0 = (bcf_gt_allele(gt_arr_ref[i+0])==1);
				bool a1 = (bcf_gt_allele(gt_arr_ref[i+1])==1);
				n_ref_missing += (gt_arr_ref[i+0] == bcf_gt_missing || gt_arr_ref[i+1] == bcf_gt_missing);
				n_ref_unphased += !bcf_gt_is_phased(gt_arr_ref[i+1]);
				H.H_opt_hap.set(i+2*n_main_samples+0, i_variant, a0);
				H.H_opt_hap.set(i+2*n_main_samples+1, i_variant, a1);
				a0?calt++:cref++;
				a1?calt++:cref++;
			}
			if (line_scaf=bcf_sr_get_line(sr, 2)) {
				ngt_scaf = bcf_get_genotypes(sr->readers[2].header, line_scaf, &gt_arr_scaf, &ngt_arr_scaf); assert(ngt_scaf == 2 * n_scaf_samples);
				for(int i = 0 ; i < 2 * n_scaf_samples ; i += 2) {
					int ind = mappingS2G[DIV2(i)];
					if (ind>=0) {
						bool s0 = (bcf_gt_allele(gt_arr_scaf[i+0])==1);
						bool s1 = (bcf_gt_allele(gt_arr_scaf[i+1])==1);
						bool he = (s0 != s1);
						bool ph = (bcf_gt_is_phased(gt_arr_scaf[i+0]) || bcf_gt_is_phased(gt_arr_scaf[i+1]));
						bool mi = (gt_arr_scaf[i+0] == bcf_gt_missing || gt_arr_scaf[i+1] == bcf_gt_missing);
						if (he && !mi && ph) {
							bool a0 = VAR_GET_HAP0(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
							bool a1 = VAR_GET_HAP1(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
							if (a0!=a1) {
								VAR_SET_SCA(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
								s0?VAR_SET_HAP0(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]):VAR_CLR_HAP0(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
								s1?VAR_SET_HAP1(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]):VAR_CLR_HAP1(MOD2(i_variant), G.vecG[ind]->Variants[DIV2(i_variant)]);
								n_geno_sca ++;
							}
						}
					}
				}
			}
			newV->cref = cref;newV->calt = calt;newV->cmis = cmis;
			i_variant ++;
			V.push(newV);
			vrb.progress("  * VCF/BCF parsing", i_variant*1.0/n_variants);
		}
	}
	free(gt_arr_main);
	free(gt_arr_scaf);
	free(gt_arr_ref);
	if (ps_arr_main) free(ps_arr_main);
	bcf_sr_destroy(sr);
	// Report
	n_geno_tot = n_main_samples*n_variants;
	string str0 = "Hom=" + stb.str(n_geno_hom*100.0/(n_main_samples*n_variants), 1) + "%";
	string str1 = "Het=" + stb.str(n_geno_het*100.0/n_geno_tot, 1) + "%" + (use_PS_field?(" / Pha=" + stb.str(n_geno_ips*100.0/n_geno_tot, 3) + "%"):(""));
	string str2 = "Sca=" + stb.str(n_geno_sca*100.0/n_geno_tot, 3) + "%";
	string str3 = "Mis=" + stb.str(n_geno_mis*100.0/(n_main_samples*n_variants), 1) + "%";
	string str4 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing ["+str0+" / "+str1+" / "+str2+" / "+str3+"] ("+str4+")");
	if (n_ref_missing > 0) vrb.warning(stb.str(n_ref_missing) + " missing genotypes in the reference panel (randomly imputed)");
	if (n_ref_unphased > 0) vrb.warning(stb.str(n_ref_unphased) + " unphased genotypes in the reference panel (randomly phased)");
}
