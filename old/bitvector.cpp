#include <objects/bitvector.h>

bitvector::bitvector(unsigned int _n) {
	n = _n;
	data = vector < unsigned int > ((unsigned int)ceil(n * 1.0 / 32), 0);
}

bitvector::~bitvector() {
	n = 0;
	vector < unsigned int > ().swap(data);
}
