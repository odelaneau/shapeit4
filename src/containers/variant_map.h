#ifndef _SNP_SET_H
#define _SNP_SET_H

#include <utils/otools.h>
#include <objects/variant.h>
#include <io/gmap_reader.h>

class variant_map {
public :
	//DATA
	vector < variant * > vec_pos;			//vector of variants ordered by position in bp
	multimap < int, variant * > map_pos;	//associative container of variant with position in bp

	//CONSTRUCTOR/DESTRUCTOR
	variant_map();
	~variant_map();

	//METHODS
	int size();
	vector < variant * > getByPos(int);
	vector < variant * > getByRef(int, string &, string &);
	variant * getByIndex(int);
	void push(variant *);
	void setGeneticMap(gmap_reader&);
	int setCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM);
	int interpolateCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM);
	unsigned int length();
};

#endif
