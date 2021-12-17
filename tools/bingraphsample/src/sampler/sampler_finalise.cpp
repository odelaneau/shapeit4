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
#include <sampler/sampler_header.h>

#include <io/haplotype_writer.h>

#define MODE_SPL	0
#define MODE_SOL	1
#define MODE_COL	2

void sampler::write_files_and_finalise() {
	vrb.title("Finalization:");

	//
	int mode;
	G.init();
	if (options.count("sample")) {
		G.sample();
		mode = MODE_SPL;
	}
	if (options.count("solve")) {
		G.solve();
		mode = MODE_SOL;
	}
	if (options.count("collapse")) {
		G.collapse(options["collapse"].as < int > (), options["thread"].as < int > ());
		mode = MODE_COL;
	}



	//step1: writing best guess haplotypes in VCF/BCF file
	haplotype_writer(G, V).writeHaplotypes(options["output"].as < string > (), mode, options["seed"].as < int > ());

	if (options["thread"].as < int > () > 1) pthread_mutex_destroy(&G.mutex_workers);

	//step2: Measure overall running time
	vrb.bullet("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
