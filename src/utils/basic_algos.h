#ifndef _BASIC_ALGOS_H
#define _BASIC_ALGOS_H

#include <vector>

class basic_algos {
public:
	basic_algos () {};
	~basic_algos () {};

	template < class T >
	unsigned int imax(std::vector < T > & vec) {
		T maxValue = vec[0];
		int maxIndex = 0;
		for (unsigned int i = 1; i < vec.size() ; i ++) if (vec[i] > maxValue) {
			maxValue = vec[i];
			maxIndex = i;
		}
		return maxIndex;
	}
};

#endif

