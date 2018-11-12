#ifndef _OLIVIER_TOOLS_H
#define _OLIVIER_TOOLS_H

//INCLUDE STANDARD TEMPLATE LIBRARY USEFULL STUFFS (STL)
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <bitset>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <string>
#include <exception>
#include <cassert>
#include <limits>
#include <cstdint>

//INCLUDE BOOST USEFULL STUFFS (BOOST)
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>

//INCLUDE HTS LIBRARY
#include <htslib/hts.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
extern "C" {
	#include <htslib/vcf_sweep.h>
	#include <htslib/synced_bcf_reader.h>
	#include <htslib/vcf.h>
	#include <htslib/vcfutils.h>
}

//INCLUDES BASE STUFFS
#include <utils/compressed_io.h>
#include <utils/random_number.h>
#include <utils/basic_stats.h>
#include <utils/basic_algos.h>
#include <utils/string_utils.h>
#include <utils/timer.h>
#include <utils/verbose.h>

//MACROS
#define DIV2(v)	(v>>1)
#define MOD2(v)	(v&1)

//NAMESPACE
using namespace std;
namespace bio = boost::iostreams;
namespace bpo = boost::program_options;
namespace bid = boost::uuids;

//MAKE SOME TOOL FULLY ACCESSIBLE THROUGHOUT THE SOFTWARE
#ifdef _DECLARE_TOOLBOX_HERE
	random_number_generator rng;	//Random number generator
	string_utils stb;				//String manipulation
	basic_algos alg;				//Basic algorithms
	verbose vrb;					//Verbose
	timer tac;						//Timer
#else
	extern random_number_generator rng;
	extern string_utils stb;
	extern basic_algos alg;
	extern verbose vrb;
	extern timer tac;
#endif

#endif
