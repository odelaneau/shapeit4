//$Id: chaplotype.h 757 2014-02-24 07:30:13Z koskos $

#ifndef _BITVECTOR_H
#define _BITVECTOR_H

#include <utils/otools.h>

class bitvector {
public:
	vector < unsigned int > data;
	unsigned int n;

	//CONSTRUCTOR/DESTRUCTOR
	bitvector(unsigned int);
	~bitvector();

	//METHODS
	bool operator [] (const unsigned int);
	void set(const unsigned int);
	void clr(const unsigned int);
};

inline
bool bitvector::operator [] (const unsigned int l) {
	return data[l>>5]&(1<<(l&31));
}

inline
void bitvector::set(const unsigned int l) {
	data[l>>5]|=(1<<(l&31));
}

inline
void bitvector::clr(const unsigned int l) {
	data[l>>5]&=~(1<<(l&31));
}

#endif
