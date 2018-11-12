#ifndef _BASIC_STATS_H
#define _BASIC_STATS_H

#include <vector>

class basic_stats {
protected:
	uint32_t m_n;
	double m_oldM;
	double m_newM;
	double m_oldS;
	double m_newS;

public:
	basic_stats() {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
	}

	template <class T>
	basic_stats(std::vector < T > & X) {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
		for (uint32_t e = 0 ; e < X.size() ; e ++) push(X[e]);
	}

	void clear() {
		m_n = 0;
		m_oldM = 0;
		m_newM = 0;
		m_oldS = 0;
		m_newS = 0;
	}

	template <class T>
	void push(T x) {
		m_n++;
		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		} else {
			m_newM = m_oldM + (x - m_oldM)/m_n;
            m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
            m_oldM = m_newM;
            m_oldS = m_newS;
		}
	}

	int size() const {
		return m_n;
	}

	double mean() const {
		return (m_n > 0) ? m_newM : 0.0;
	}

	double variance() const {
		return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
	}

	double sd() const {
		return sqrt( variance() );
	}
};

#endif
