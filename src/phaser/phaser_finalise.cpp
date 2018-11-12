#include <phaser/phaser_header.h>

#include <io/haplotype_writer.h>

void phaser::write_files_and_finalise() {
	vrb.title("Finalization:");

	//step0: multi-threading
	if (options["thread"].as < int > () > 1) pthread_mutex_destroy(&mutex_workers);

	//
	G.solve();
	H.update(G);
	H.transposeH2V(false);

	//step1: writing best guess haplotypes in VCF/BCF file
	haplotype_writer(H, G, V).writeHaplotypes(options["output"].as < string > ());

	//step2: Measure overall running time
	vrb.bullet("Total running time = " + stb.str(tac.abs_time()) + " seconds");
}
