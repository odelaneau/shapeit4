//$Id: builder.cpp 640 2012-11-28 17:13:55Z koskos $

#include <modules/builder.h>


builder::builder(genotype_set & _G, int n_thread): G(_G) {
	this->n_thread = n_thread;
	if (n_thread > 1) {
		i_workers = 0;
		id_workers = vector < pthread_t > (n_thread);
		pthread_mutex_init(&mutex_workers, NULL);
	}
}

builder::~builder() {
	if (n_thread > 1) {
		i_workers = 0;
		pthread_mutex_destroy(&mutex_workers);
		id_workers.clear();
	}
}

void * builder_callback(void * ptr) {
	builder * B = static_cast< builder * >( ptr );
	for(;;) {
		pthread_mutex_lock( &B->mutex_workers );
		int curr_ind_to_process = B->i_workers++;
		pthread_mutex_unlock( &B->mutex_workers);
		if (curr_ind_to_process < B->G.n_ind) B->build(curr_ind_to_process);
		else pthread_exit(NULL);
	}
	return NULL;
}

void builder::build(int ind) {
	G.vecG[ind]->build();
}

void builder::build() {
	tac.clock();
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, builder_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i  <  G.n_ind ; i ++) build(i);
	long int n_segments = G.numberOfSegments();
	vrb.bullet("Build genotype graphs [seg=" + stb.str(n_segments) + "] (" + stb.str(tac.rel_time()*0.001, 2) + "s)");
}

