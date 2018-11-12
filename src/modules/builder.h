#ifndef _BUILDER_H
#define _BUILDER_H

#include <utils/otools.h>
#include <containers/genotype_set.h>

class builder {
public:
	//DATA
	genotype_set & G;

	//MULTI-THREADING
	int i_workers;
	int n_thread;
	pthread_mutex_t mutex_workers;
	vector < pthread_t > id_workers;

	//CONSTRUCTOR/DESTRUCTOR
	builder(genotype_set &, int n_thread = 1);
	~builder();

	//METHODS
	void build();
	void build(int);
};

#endif
