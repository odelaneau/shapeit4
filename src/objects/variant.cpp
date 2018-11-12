//$Id: snp.cpp 640 2012-11-28 17:13:55Z koskos $

#define _GLOBAL
#include <objects/variant.h>

variant::variant(string & chr, int bp, string & id, string & ref, string & alt, int idx) {
	this->chr = chr;
	this->bp = bp;
	this->id = id;
	this->ref = ref;
	this->alt = alt;
	this->idx = idx;
	this->cref = 0;
	this->calt = 0;
	this->cmis = 0;
	this->fbin = 0;
	cm = -1;
}

variant::~variant() {
}

unsigned int variant::getMAC() {
	return min(cref, calt);
}

bool variant::isSingleton() {
	return (calt == 1 || cref == 1);
}

bool variant::isMonomorphic() {
	return (calt == 0 || cref == 0);
}
