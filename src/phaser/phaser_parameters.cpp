#include <phaser/phaser_header.h>

void phaser::declare_options() {
	bpo::options_description opt_base ("Basic options");
	opt_base.add_options()
			("help", "Produce help message")
			("seed", bpo::value<int>()->default_value(15052011), "Seed of the random number generator")
			("thread,T", bpo::value<int>()->default_value(1), "Number of thread used");

	bpo::options_description opt_input ("Input files");
	opt_input.add_options()
			("input,I", bpo::value< string >(), "Genotypes to be phased in VCF/BCF format")
			("reference,H", bpo::value< string >(), "Reference panel of haplotypes in VCF/BCF format")
			("scaffold,S", bpo::value< string >(), "Scaffold of haplotypes in VCF/BCF format")
			("map,M", bpo::value< string >(), "Genetic map")
			("region,R", bpo::value< string >(), "Target region")
			("use-PS", bpo::value<double>(), "Informs phasing using PS field from read based phasing");

	bpo::options_description opt_mcmc ("MCMC parameters");
	opt_mcmc.add_options()
			("mcmc-iterations", bpo::value<string>()->default_value("5b,1p,1b,1p,1b,1p,5m"), "Iteration scheme of the MCMC")
			("mcmc-prune", bpo::value<double>()->default_value(0.999), "Pruning threshold")
			("mcmc-store-K", bpo::value<string>(), "Store K sizes in last iterations");

	bpo::options_description opt_pbwt ("PBWT parameters");
	opt_pbwt.add_options()
			("pbwt-disable-init", "Do not initialise haplotypes by PBWT (rephase input haplotype data)")
			("pbwt-modulo", bpo::value<int>()->default_value(8), "Storage frequency of PBWT indexes in variant numbers (i.e. 16 means storage every 16 variants)")
			("pbwt-depth", bpo::value<int>()->default_value(4), "Depth of PBWT indexes to condition on");
	
	bpo::options_description opt_hmm ("HMM parameters");
	opt_hmm.add_options()
			("window,W", bpo::value<double>()->default_value(2e6), "Minimal size of the phasing window")
			("effective-size", bpo::value<int>()->default_value(15000), "Effective size of the population");

	bpo::options_description opt_output ("Output files");
	opt_output.add_options()
			("output,O", bpo::value< string >(), "Phased haplotypes in VCF/BCF format")
			("log", bpo::value< string >(), "Log file");

	descriptions.add(opt_base).add(opt_input).add(opt_mcmc).add(opt_pbwt).add(opt_hmm).add(opt_output);
}

void phaser::parse_command_line(vector < string > & args) {
	try {
		bpo::store(bpo::command_line_parser(args).options(descriptions).run(), options);
		bpo::notify(options);
	} catch ( const boost::program_options::error& e ) { cerr << "Error parsing command line arguments: " << string(e.what()) << endl; exit(0); }

	if (options.count("help")) { cout << descriptions << endl; exit(0); }

	if (options.count("log") && !vrb.open_log(options["log"].as < string > ()))
		vrb.error("Impossible to create log file [" + options["log"].as < string > () +"]");

	vrb.title("SHAPEIT");
	vrb.bullet("Authors       : Olivier DELANEAU");
	vrb.bullet("Copyright (C) : University of Lausanne");
	vrb.bullet("Contact       : olivier.delaneau@gmail.com");
	vrb.bullet("Version       : 4.0.0");
	vrb.bullet("Run date      : " + tac.date());
}

void phaser::check_options() {
	if (!options.count("input"))
		vrb.error("You must specify one input file using --input");

	if (!options.count("region"))
		vrb.error("You must specify a region or chromosome to phase using --region");

	if (!options.count("output"))
		vrb.error("You must specify a phased output file with --output");

	if (!options.count("map"))
		vrb.error("You must specify a genetic map file with --map");

	if (options.count("seed") && options["seed"].as < int > () < 0)
		vrb.error("Random number generator needs a positive seed value");

	if (options.count("thread") && options["thread"].as < int > () < 1)
		vrb.error("You must use at least 1 thread");

	if (!options["thread"].defaulted() && !options["seed"].defaulted())
		vrb.warning("Using multi-threading prevents reproducing a run by specifying --seed");

	if (!options["effective-size"].defaulted() && options["effective-size"].as < int > () < 1)
		vrb.error("You must specify a positive effective size");

	if (!options["window"].defaulted() && options["window"].as < double > () < 1e5)
		vrb.error("You must specify a window size of at least 0.1 Mb");

	parse_iteration_scheme(options["mcmc-iterations"].as < string > ());
}

void phaser::verbose_files() {
	vrb.title("Files:");
	vrb.bullet("Input VCF     : [" + options["input"].as < string > () + "]");
	if (options.count("reference")) vrb.bullet("Reference VCF : [" + options["reference"].as < string > () + "]");
	if (options.count("scaffold")) vrb.bullet("Scaffold VCF  : [" + options["scaffold"].as < string > () + "]");
	vrb.bullet("Genetic Map   : [" + options["map"].as < string > () + "]");
	vrb.bullet("Output VCF    : [" + options["output"].as < string > () + "]");
	if (options.count("log")) vrb.bullet("Output LOG    : [" + options["log"].as < string > () + "]");
}

void phaser::verbose_options() {
	vrb.title("Parameters:");
	vrb.bullet("Seed    : " + stb.str(options["seed"].as < int > ()));
	vrb.bullet("Threads : " + stb.str(options["thread"].as < int > ()) + " threads");
	vrb.bullet("MCMC    : " + get_iteration_scheme());
	if (options.count("pbwt-disable-init")) vrb.bullet("PBWT    : No PBWT initialization");
	vrb.bullet("PBWT    : Store indexes every " + stb.str(options["pbwt-modulo"].as < int > ()) + " variants");
	vrb.bullet("PBWT    : Depth of PBWT neighbours to condition on: " + stb.str(options["pbwt-depth"].as < int > ()));
	vrb.bullet("HMM     : K is variable / min W is " + stb.str(options["window"].as < double > ()/1e6, 2) + "Mb / Ne is "+ stb.str(options["effective-size"].as < int > ()));
	if (options.count("use-PS")) vrb.bullet("HMM     : Inform phasing using VCF/PS field / Error rate of PS field is " + stb.str(options["use-PS"].as < double > ()));
}
