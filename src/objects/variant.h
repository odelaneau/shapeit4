//$Id: snp.h 640 2012-11-28 17:13:55Z koskos $

#ifndef _VARIANT_H
#define _VARIANT_H

#include <utils/otools.h>

class variant {
public :
	//DATA
	string chr;
	int bp;
	string id;
	string ref;
	string alt;
	double cm;
	int idx;
	unsigned int cref;
	unsigned int calt;
	unsigned int cmis;
	unsigned int fbin;

	//CONSTRUCTOR/DESTRUCTOR
	variant(string & chr, int bp, string & id, string & ref, string & alt, int idx);
	~variant();

	bool isSingleton();
	bool isMonomorphic();
	unsigned int getMAC();
};

#endif
