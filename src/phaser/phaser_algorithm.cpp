#include <phaser/phaser_header.h>

#include <io/haplotype_writer.h>

void * phaseWindow_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_worker, id_job;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_ind) vrb.progress("  * HMM computations", id_job*1.0/S->G.n_ind);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_ind) S->phaseWindow(id_worker, id_job);
		else pthread_exit(NULL);
	}
}

void phaser::phaseWindow(int id_worker, int id_job) {
	threadData[id_worker].make(id_job, options["window"].as < double > ());
	for (int w = 0 ; w < threadData[id_worker].size() ; w ++) {
		if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
		statH.push(threadData[id_worker].Kvec[w].size()*1.0);
		statS.push((V.vec_pos[threadData[id_worker].C[w].stop_locus]->bp - V.vec_pos[threadData[id_worker].C[w].start_locus]->bp + 1) * 1.0 / 1e6);
		if (options.count("mcmc-store-K")) storedKsizes.push_back(threadData[id_worker].Kvec[w].size()*1.0);
		if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);
		assert(threadData[id_worker].Kvec[w].size()>0);

		haplotype_segment HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kvec[w], threadData[id_worker].C[w], M);
		int outcome = HS.expectation(threadData[id_worker].T);
		if (outcome < 0) vrb.error("Underflow impossible to recover");
		else n_underflow_recovered += outcome;
	}

	if (options.count("use-PS") && G.vecG[id_job]->ProbabilityMask.size() > 0) threadData[id_worker].maskingTransitions(id_job, options["use-PS"].as < double > ());

	vector < bool > flagMerges;
	switch (iteration_types[iteration_stage]) {
	case STAGE_BURN:	G.vecG[id_job]->sample(threadData[id_worker].T);
						break;
	case STAGE_PRUN:	G.vecG[id_job]->sample(threadData[id_worker].T);
						G.vecG[id_job]->mapMerges(threadData[id_worker].T, options["mcmc-prune"].as < double > (), flagMerges);
						G.vecG[id_job]->performMerges(threadData[id_worker].T, flagMerges);
						break;
	case STAGE_MAIN:	G.vecG[id_job]->sample(threadData[id_worker].T);
						G.vecG[id_job]->store(threadData[id_worker].T);
						break;
	}
}

void phaser::phaseWindow() {
	tac.clock();
	int n_thread = options["thread"].as < int > ();
	n_underflow_recovered = 0;
	i_workers = 0; i_jobs = 0;
	statH.clear(); statS.clear();
	storedKsizes.clear();
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, phaseWindow_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_ind ; i ++) {
		phaseWindow(0, i);
		vrb.progress("  * HMM computations", (i+1)*1.0/G.n_ind);
	}
	if (n_underflow_recovered) vrb.bullet("HMM computations [K=" + stb.str(statH.mean(), 1) + "+/-" + stb.str(statH.sd(), 1) + " / W=" + stb.str(statS.mean(), 2) + "Mb / U=" + stb.str(n_underflow_recovered) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	else vrb.bullet("HMM computations [K=" + stb.str(statH.mean(), 1) + "+/-" + stb.str(statH.sd(), 1) + " / W=" + stb.str(statS.mean(), 2) + "Mb] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void phaser::phase() {
	unsigned long n_old_segments = G.numberOfSegments(), n_new_segments = 0, current_iteration = 0;
	for (iteration_stage = 0 ; iteration_stage < iteration_counts.size() ; iteration_stage ++) {
		for (int iter = 0 ; iter < iteration_counts[iteration_stage] ; iter ++) {
			switch (iteration_types[iteration_stage]) {
			case STAGE_BURN:	vrb.title("Burn-in iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			case STAGE_PRUN:	vrb.title("Pruning iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			case STAGE_MAIN:	vrb.title("Main iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			}
			H.transposeV2H(false);
			H.select();
			H.transposeC2H();
			phaseWindow();
			if (options.count("mcmc-store-K")) {
				string filename = options["mcmc-store-K"].as < string > () + stb.str(current_iteration) + ".txt.gz";
				output_file outputK(filename);
				for (int e = 0 ; e < storedKsizes.size() ; e++) outputK << storedKsizes[e] << endl;
				outputK.close();
				current_iteration++;
			}

			H.update(G);
			H.transposeH2V(false);
			if (iteration_types[iteration_stage] == STAGE_PRUN) {
				n_new_segments = G.numberOfSegments();
				//vrb.bullet("Pruning info [old=" + stb.str(n_old_segments) + " / new=" + stb.str(n_new_segments) + " / compression=" + stb.str((1-n_new_segments*1.0/n_old_segments)*100, 2) + "%]");
				vrb.bullet("Pruning outcome [compression=" + stb.str((1-n_new_segments*1.0/n_old_segments)*100, 2) + "%]");
				if (options.count("use-PS")) G.masking();
			}
		}
	}
}
