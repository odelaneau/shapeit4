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
	gmap_reader readerGM;
	readerGM.readGeneticMapFile(options["map"].as < string > ());
	V.setGeneticMap(readerGM);
	M.initialise(V, options["effective-size"].as < int > (), (readerG.n_main_samples+readerG.n_ref_samples*2), 100UL, options.count("map"));

	//step4: Initialize haplotypes
	H.allocate(V, options["pbwt-modulo"].as < int > (), options["pbwt-depth"].as < int > ());
	H.update(G, true);
	H.transposeH2V(true);
	H.searchIBD2((int)round((options["window"].as < double > () * V.size()) / V.length()));
	if (!options.count("pbwt-disable-init")) {
		pbwt_solver solver = pbwt_solver(H);
		solver.sweep(G);
		solver.free();
	}

	//step5: Initialize genotype structures
	builder(G, options["thread"].as < int > ()).build();
	if (options.count("use-PS")) G.masking();

	//step6: Allocate data structures for computations
	unsigned int max_number_transitions = G.largestNumberOfTransitions();
	threadData = vector < compute_job >(options["thread"].as < int > (), compute_job(V, G, H, max_number_transitions));
}
