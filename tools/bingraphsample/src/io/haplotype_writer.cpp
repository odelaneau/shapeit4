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
#include <io/haplotype_writer.h>

#define OFILE_VCFU	0
#define OFILE_VCFC	1
#define OFILE_BCFC	2

#define MODE_SPL	0
#define MODE_SOL	1
#define MODE_COL	2



haplotype_writer::haplotype_writer(genotype_set & _G, variant_map & _V): G(_G), V(_V) {
}

haplotype_writer::~haplotype_writer() {
}

void haplotype_writer::writeHaplotypes(string fname, int mode, int seed) {
	// Init
	tac.clock();
	string file_format = "w";
	unsigned int file_type = OFILE_VCFU;
	if (fname.size() > 6 && fname.substr(fname.size()-6) == "vcf.gz") { file_format = "wz"; file_type = OFILE_VCFC; }
	if (fname.size() > 3 && fname.substr(fname.size()-3) == "bcf") { file_format = "wb"; file_type = OFILE_BCFC; }
	htsFile * fp = hts_open(fname.c_str(),file_format.c_str());
	bcf_hdr_t * hdr = bcf_hdr_init("w");
	bcf1_t *rec = bcf_init1();

	// Create VCF header
	bcf_hdr_append(hdr, string("##fileDate="+tac.date()).c_str());
	string str_mode = "";
	switch (mode) {
	case MODE_SPL:	str_mode = "##source=shapeit4-sample-" + stb.str(seed); break;
	case MODE_SOL:	str_mode = "##source=shapeit4-best_guess"; break;
	case MODE_COL:	str_mode = "##source=shapeit4-collapse"; break;
	}
	bcf_hdr_append(hdr, str_mode.c_str());
	bcf_hdr_append(hdr, string("##contig=<ID="+ V.vec_pos[0]->chr + ">").c_str());
	bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">");
	bcf_hdr_append(hdr, "##INFO=<ID=AC,Number=1,Type=Integer,Description=\"Allele count\">");
	bcf_hdr_append(hdr, "##INFO=<ID=CM,Number=A,Type=Float,Description=\"Interpolated cM position\">");
	bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Phased genotypes\">");

	//Add samples
	for (int i = 0 ; i < G.n_ind ; i ++) bcf_hdr_add_sample(hdr, G.vecG[i]->name.c_str());
	bcf_hdr_add_sample(hdr, NULL);      // to update internal structures
	bcf_hdr_write(fp, hdr);

	//Add records
	if (mode != MODE_COL) {
		int * genotypes = (int*)malloc(bcf_hdr_nsamples(hdr)*2*sizeof(int));
		for (int l = 0 ; l < V.size() ; l ++) {
			bcf_clear1(rec);
			rec->rid = bcf_hdr_name2id(hdr, V.vec_pos[l]->chr.c_str());
			rec->pos = V.vec_pos[l]->bp - 1;
			bcf_update_id(hdr, rec, V.vec_pos[l]->id.c_str());
			string alleles = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
			bcf_update_alleles_str(hdr, rec, alleles.c_str());
			int count_alt = 0;
			for (int i = 0 ; i < G.n_ind ; i++) {
				bool a0 = G.vecG[i]->H0[l];
				bool a1 = G.vecG[i]->H1[l];
				count_alt += a0+a1;
				genotypes[2*i+0] = bcf_gt_phased(a0);
				genotypes[2*i+1] = bcf_gt_phased(a1);
			}
			bcf_update_info_int32(hdr, rec, "AC", &count_alt, 1);
			float freq_alt = count_alt * 1.0 / (2 * G.n_ind);
			bcf_update_info_float(hdr, rec, "AF", &freq_alt, 1);
			if (V.vec_pos[l]->cm >= 0) {
				float val = (float)V.vec_pos[l]->cm;
				bcf_update_info_float(hdr, rec, "CM", &val, 1);
			}
			bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2);
			bcf_write1(fp, hdr, rec);
			vrb.progress("  * VCF writing", (l+1)*1.0/V.size());
		}
		free(genotypes);
	} else {
		float * posteriors = (float*)malloc(bcf_hdr_nsamples(hdr)*4*sizeof(float));
		for (int l = 0 ; l < V.size() ; l ++) {
			bcf_clear1(rec);
			rec->rid = bcf_hdr_name2id(hdr, V.vec_pos[l]->chr.c_str());
			rec->pos = V.vec_pos[l]->bp - 1;
			bcf_update_id(hdr, rec, V.vec_pos[l]->id.c_str());
			string alleles = V.vec_pos[l]->ref + "," + V.vec_pos[l]->alt;
			bcf_update_alleles_str(hdr, rec, alleles.c_str());
			float count_alt = 0;
			for (int i = 0 ; i < G.n_ind ; i++) {
				float rr = G.vecG[i]->Collapsed[4*l+0];
				float ra = G.vecG[i]->Collapsed[4*l+1];
				float ar = G.vecG[i]->Collapsed[4*l+2];
				float aa = G.vecG[i]->Collapsed[4*l+3];
				float scale = 1.0/(rr+ra+ar+aa);
				rr *= scale;
				ra *= scale;
				ar *= scale;
				aa *= scale;
				count_alt += ra+ar+2.0*aa;
				posteriors[4*i+0] = rr;
				posteriors[4*i+1] = ra;
				posteriors[4*i+2] = ar;
				posteriors[4*i+3] = aa;
			}
			int count_alt_int = (int)count_alt;
			bcf_update_info_int32(hdr, rec, "AC", &count_alt_int, 1);
			float freq_alt = count_alt * 1.0 / (2 * G.n_ind);
			bcf_update_info_float(hdr, rec, "AF", &freq_alt, 1);
			if (V.vec_pos[l]->cm >= 0) {
				float val = (float)V.vec_pos[l]->cm;
				bcf_update_info_float(hdr, rec, "CM", &val, 1);
			}
			bcf_update_format_float(hdr, rec, "PO", posteriors, bcf_hdr_nsamples(hdr)*4);
			bcf_write1(fp, hdr, rec);
			vrb.progress("  * VCF writing", (l+1)*1.0/V.size());
		}
		free(posteriors);
	}
	bcf_destroy1(rec);
	bcf_hdr_destroy(hdr);
	if (hts_close(fp)) vrb.error("Non zero status when closing VCF/BCF file descriptor");
	switch (file_type) {
	case OFILE_VCFU: vrb.bullet("VCF writing [Uncompressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_VCFC: vrb.bullet("VCF writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	case OFILE_BCFC: vrb.bullet("BCF writing [Compressed / N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)"); break;
	}
}
