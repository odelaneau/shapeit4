#ifndef _GMAP_DATA_H
#define _GMAP_DATA_H

#include <utils/otools.h>

class gmap_reader {
public:
	//DATA
	vector < int > pos_bp;
	vector < double > pos_cm;

	//CONSTRUCTOR/DESTRUCTOR
	gmap_reader();
	~gmap_reader();

	//IO
	void readGeneticMapFile(string);
};

#endif
