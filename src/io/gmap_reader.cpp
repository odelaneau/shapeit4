//$Id: gmap_reader.cpp 645 2012-12-03 15:22:36Z koskos $

#include <io/gmap_reader.h>

gmap_reader::gmap_reader() {
}

gmap_reader::~gmap_reader() {
	vector < int > ().swap(pos_bp);
	vector < double > ().swap(pos_cm);
}

void gmap_reader::readGeneticMapFile(string fmap) {
	tac.clock();
	string buffer;
	vector < string > tokens;
	int line = 0;
	input_file fd_gmap(fmap);
	getline(fd_gmap, buffer, '\n');
	int prev_bp = 0;
	double prev_cm = 0;
	while (getline(fd_gmap, buffer, '\n')) {
		if (stb.split(buffer, tokens) == 3) {
			int curr_bp = atoi(tokens[0].c_str());
			double curr_cm = atof(tokens[2].c_str());
			if (curr_bp < prev_bp || curr_cm < prev_cm)
				vrb.error("Wrong order in your genetic map file " + stb.str(prev_bp) + "bp / " + stb.str(prev_cm,5) + "cM > " + stb.str(curr_bp) + "bp / " + stb.str(curr_cm,5) + "cM");
			pos_bp.push_back(curr_bp);
			pos_cm.push_back(curr_cm);
			prev_bp = curr_bp;
			prev_cm = curr_cm;
		} else vrb.error("Parsing line " + stb.str(line) + " : incorrect number of columns, observed: " + stb.str(tokens.size()) + " expected: 3");
		line++;
	}
	fd_gmap.close();
	vrb.bullet("GMAP parsing [n=" + stb.str(line) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}
