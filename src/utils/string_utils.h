#ifndef _STRING_UTILS_H
#define _STRING_UTILS_H

#include <sstream>
#include <iomanip>
#include <string>
#include <vector>

//using namespace std;

class string_utils {
public:
	string_utils () {};
	~string_utils () {};

	int split(const std::string & str, std::vector < std::string > & tokens, std::string sep = " 	", unsigned int n_max_tokens = 1000000) {
		tokens.clear();
		if (str == ""){
			tokens.push_back("");
			return tokens.size();
		}
		std::string::size_type p_last = str.find_first_not_of(sep, 0);
		std::string::size_type p_curr = str.find_first_of(sep, p_last);
		while ((std::string::npos != p_curr || std::string::npos != p_last) && tokens.size() < n_max_tokens) {
			tokens.push_back(str.substr(p_last, p_curr - p_last));
			p_last = str.find_first_not_of(sep, p_curr);
			p_curr = str.find_first_of(sep, p_last);
		}
		if (tokens.back()[tokens.back().size()-1] == '\r') tokens.back() = tokens.back().substr(0, tokens.back().size()-1);
		return tokens.size();
	}

	bool numeric(std::string & str) {
		double n;
		std::istringstream in(str);
		if (!(in >> n)) return false;
		return true;
    }

	template < class T >
	std::string str(T n, int prec = -1) {
		std::ostringstream ss( std::stringstream::out );
		if (prec >= 0) { ss << setiosflags( std::ios::fixed ); ss.precision(prec); }
		ss << n;
		return ss.str();
	}

	template < class T >
	std::string str(std::vector < T > & v, int prec = -1) {
		std::ostringstream ss( std::stringstream::out );
		if (prec >= 0) { ss << setiosflags( std::ios::fixed ); ss.precision(prec); }
		for (int e = 0 ; e < v.size() ; e ++) ss << (e>0?" ":"") << v[e] ;
		return ss.str();
	}
};

#endif
