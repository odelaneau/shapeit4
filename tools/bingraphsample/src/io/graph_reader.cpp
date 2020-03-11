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
#include <io/graph_reader.h>

graph_reader::graph_reader(genotype_set & _G, variant_map & _V): G(_G), V(_V) {
}

graph_reader::~graph_reader() {
}

void graph_reader::binary_read(input_file & fin, vector < bool > & x) {
	vector < bool >::size_type n;
	fin.read((char*)&n, sizeof(vector<bool>::size_type));
	x.resize(n);
    for(vector<bool>::size_type i = 0; i < n;) {
    	unsigned char aggr;
        fin.read((char*)&aggr, sizeof(unsigned char));
        for(unsigned char mask = 1; mask > 0 && i < n; ++i, mask <<= 1) x.at(i) = aggr & mask;
	}
}

void graph_reader::string_read(input_file & fin, string & x) {
	size_t str_len;
	fin.read((char*)&str_len, sizeof(size_t));
	char * temp = new char[str_len+1];
	fin.read(temp, str_len);
	temp[str_len] = '\0';
	x = temp;
	delete [] temp;
}

void graph_reader::readGraphs(string fname) {
	// Init
	tac.clock();
	input_file fd (fname);

	//Read variant map
	int n_variants;
	fd.read((char*)&n_variants, sizeof(int));
	for (int l = 0 ; l < n_variants ; l ++) {
		variant * v = new variant();
		string_read(fd, v->chr);
		fd.read(reinterpret_cast<char*>(&v->bp), sizeof(v->bp));
		string_read(fd, v->id);
		string_read(fd, v->ref);
		string_read(fd, v->alt);
		fd.read(reinterpret_cast<char*>(&v->idx), sizeof(v->idx));
		V.push(v);
	}
	G.n_site = V.size();

	//Read genotype graphs
	int n_samples;
	fd.read((char*)&n_samples, sizeof(int));
	for (int i  = 0 ; i < n_samples ; i++) {
		genotype * g = new genotype(i);

		// name
		string_read(fd, g->name);

		// integers
		fd.read(reinterpret_cast<char*>(&g->index), sizeof(g->index));
		fd.read(reinterpret_cast<char*>(&g->n_segments), sizeof(g->n_segments));
		fd.read(reinterpret_cast<char*>(&g->n_variants), sizeof(g->n_variants));
		fd.read(reinterpret_cast<char*>(&g->n_ambiguous), sizeof(g->n_ambiguous));
		fd.read(reinterpret_cast<char*>(&g->n_missing), sizeof(g->n_missing));
		fd.read(reinterpret_cast<char*>(&g->n_transitions), sizeof(g->n_transitions));
		fd.read(reinterpret_cast<char*>(&g->n_stored_transitionProbs), sizeof(g->n_stored_transitionProbs));
		fd.read(reinterpret_cast<char*>(&g->n_storage_events), sizeof(g->n_storage_events));

		// allocation
		g->Variants = vector < unsigned char > (DIV2(n_variants) + MOD2(n_variants), 0);
		g->Ambiguous = vector < unsigned char > (g->n_ambiguous,0);
		g->Diplotypes = vector < unsigned long > (g->n_segments, 0);
		g->Lengths = vector < unsigned short > (g->n_segments, 0);
		g->ProbStored = vector < float > (g->n_stored_transitionProbs, 0.0f);
		g->ProbMissing = vector < float > (g->n_missing * HAP_NUMBER, 0.0f);

		// vectors
		fd.read(reinterpret_cast<char*>(&g->Variants[0]), g->Variants.size());
		fd.read(reinterpret_cast<char*>(&g->Ambiguous[0]), g->Ambiguous.size());
		fd.read(reinterpret_cast<char*>(&g->Diplotypes[0]), g->Diplotypes.size() * sizeof(unsigned long));
		fd.read(reinterpret_cast<char*>(&g->Lengths[0]), g->Lengths.size() * sizeof(unsigned short));
		binary_read(fd, g->ProbMask);
		fd.read(reinterpret_cast<char*>(&g->ProbStored[0]), g->ProbStored.size() * sizeof(float));
		fd.read(reinterpret_cast<char*>(&g->ProbMissing[0]), g->ProbMissing.size() * sizeof(float));

		//
		G.vecG.push_back(g);
		G.n_ind ++;
	}
	vrb.bullet("BIN reading [N=" + stb.str(G.n_ind) + " / L=" + stb.str(V.size()) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");
}
