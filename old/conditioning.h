#ifndef _CONDITIONING_H
#define _CONDITIONING_H

#include <utils/otools.h>

#define ALLOC_FIRST	64
#define ALLOC_NEXT	16

struct cstate {
	unsigned int index;
	unsigned int parent : 1;
	unsigned int length : 31;
};

class conditioning {
public:
	vector < cstate > CH;
	int last_parent0, last_parent1;


	conditioning() {
		last_parent0 = -1;
		last_parent1 = -1;
		CH.reserve(ALLOC_FIRST);
	};

	~conditioning() { hardClear(); };

	void push0(unsigned int h, unsigned int gap);
	void push1(unsigned int h, unsigned int gap);
	void fill(vector < unsigned int > &, vector < unsigned int > &);
	bool check(unsigned);
	void softClear();
	void hardClear();
};

inline
void conditioning::push0(unsigned int h, unsigned int gap) {
	if (CH.size() == CH.capacity()) CH.reserve(CH.size() + ALLOC_NEXT);
	last_parent0 = CH.size();
	CH.push_back(cstate());
	CH[last_parent0].index = h;
	CH[last_parent0].parent = 0;
	CH[last_parent0].length = gap;
}

inline
void conditioning::push1(unsigned int h, unsigned int gap) {
	if (CH.size() == CH.capacity()) CH.reserve(CH.size() + ALLOC_NEXT);
	last_parent1 = CH.size();
	CH.push_back(cstate());
	CH[last_parent1].index = h;
	CH[last_parent1].parent = 1;
	CH[last_parent1].length = gap;
}

inline
void conditioning::fill(vector < unsigned int > & v0, vector < unsigned int > & v1) {
	unsigned int labs0 = 0, labs1 = 0;
	for (unsigned int c = 0 ; c < CH.size() ; c++) {
		if (CH[c].parent) {
			std::fill(v1.begin() + labs1, v1.begin() + labs1 + CH[c].length, CH[c].index);
			labs1 += CH[c].length;
		} else {
			std::fill(v0.begin() + labs0, v0.begin() + labs0 + CH[c].length, CH[c].index);
			labs0 += CH[c].length;
		}
	}
	for (; labs0 < v0.size() ; labs0 ++) v0[labs0] = CH[last_parent0].index;
	for (; labs1 < v1.size() ; labs1 ++) v1[labs1] = CH[last_parent1].index;
}

inline
bool conditioning::check(unsigned size) {
	unsigned int labs0 = 0, labs1 = 0;
	for (unsigned int c = 0 ; c < CH.size() ; c++) {
		if (CH[c].parent) labs1 += CH[c].length;
		else labs0 += CH[c].length;
	}
	if (labs0 != size || labs1 != size) {
		cout << labs0 << " " << labs1 << " / " << size << endl;
	}
	return (labs0 != size || labs1 != size);
}

inline
void conditioning::softClear() {
	CH.clear();
	last_parent0 = -1;
	last_parent1 = -1;
}

inline
void conditioning::hardClear() {
	vector < cstate > ().swap(CH);
	last_parent0 = -1;
	last_parent1 = -1;

}

#endif

