//$Id: phaser.h 617 2012-10-10 12:57:33Z koskos $

#ifndef _PHASER_H
#define _PHASER_H

#include <utils/otools.h>
#include <objects/hmm_parameters.h>
#include <models/haplotype_segment.h>

#include <containers/genotype_set.h>
#include <containers/haplotype_set.h>
#include <containers/variant_map.h>


#define STAGE_BURN	0
#define STAGE_PRUN	1
#define STAGE_MAIN	2

class phaser {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//INTERNAL DATA
	haplotype_set H;
	genotype_set G;
	hmm_parameters M;
	variant_map V;

	//MULTI-THREADING
	int i_workers, i_jobs;
	vector < pthread_t > id_workers;
	pthread_mutex_t mutex_workers;
	vector < compute_job > threadData;

	//MCMC
	vector < unsigned int > iteration_types;
	vector < unsigned int > iteration_counts;
	unsigned int iteration_stage;
	int n_underflow_recovered;

	//
	basic_stats statH,statS;
	vector < double > storedKsizes;

	//CONSTRUCTOR
	phaser();
	~phaser();

	//METHODS
	void phase();
	void phaseWindow(int, int);
	void phaseWindow();

	//PARAMETERS
	void declare_options();
	void parse_command_line(vector < string > &);
	void parse_iteration_scheme(string);
	string get_iteration_scheme();
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void read_files_and_initialise();
	void phase(vector < string > &);
	void write_files_and_finalise();
};


#endif


