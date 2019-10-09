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
#include <phaser/phaser_header.h>

#include <io/genotype_reader.h>
#include <io/haplotype_writer.h>
#include <io/gmap_reader.h>

#include <modules/builder.h>
#include <modules/pbwt_solver.h>

void phaser::read_files_and_initialise() {
	vrb.title("Initialization:");

	//step0: Initialize seed and multi-threading
	rng.setSeed(options["seed"].as < int > ());
	if (options["thread"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step2: Read input files
	genotype_reader readerG(H, G, V, options["region"].as < string > (), options.count("use-PS"));
	if (!options.count("reference")) readerG.scanGenotypes(options["input"].as < string > ());
	else readerG.scanGenotypes(options["input"].as < string > (), options["reference"].as < string > ());
	readerG.allocateGenotypes();
	if (!options.count("reference") && !options.count("scaffold")) readerG.readGenotypes0(options["input"].as < string > ());
	if ( options.count("reference") && !options.count("scaffold")) readerG.readGenotypes1(options["input"].as < string > (), options["reference"].as < string > ());
	if (!options.count("reference") &&  options.count("scaffold")) readerG.readGenotypes2(options["input"].as < string > (), options["scaffold"].as < string > ());
	if ( options.count("reference") &&  options.count("scaffold")) readerG.readGenotypes3(options["input"].as < string > (), options["reference"].as < string > (), options["scaffold"].as < string > ());
	G.imputeMonomorphic(V);

	//step3: Read and initialise genetic map
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	M.initialise(V, options["effective-size"].as < int > (), (readerG.n_main_samples+readerG.n_ref_samples)*2);

	//step4: Initialize haplotypes
	H.parametrizePBWT(options["pbwt-depth"].as < int > (), options["pbwt-modulo"].as < double > (), options["pbwt-mac"].as < int > (), options["pbwt-mdr"].as < double > (), options["thread"].as < int > ());
	H.initializePBWTmapping(V);
	H.allocatePBWTarrays();
	H.updateHaplotypes(G, true);
	H.transposeHaplotypes_H2V(true);
	H.searchIBD2matching(V, options["ibd2-length"].as < double > (), options["window"].as < double > (), options["ibd2-maf"].as < double > (), options["ibd2-mdr"].as < double > (), options["ibd2-count"].as < int > ());
	pbwt_solver solver = pbwt_solver(H);
	solver.sweep(G);
	solver.free();

	//step5: Initialize genotype structures
	builder(G, options["thread"].as < int > ()).build();
	if (options.count("use-PS")) G.masking();

	//step6: Allocate data structures for computations
	unsigned int max_number_transitions = G.largestNumberOfTransitions();
	threadData = vector < compute_job >(options["thread"].as < int > (), compute_job(V, G, H, max_number_transitions));
}
