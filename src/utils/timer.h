#ifndef _TIMER_H
#define _TIMER_H

#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include <string>

class timer {
protected:
	std::chrono::time_point<std::chrono::high_resolution_clock> start_timing_clock, prev_timing_clock;

public:
	timer () {
		start_timing_clock = std::chrono::high_resolution_clock::now();
	}

	~timer() {
	}

	void clock() {
		prev_timing_clock = std::chrono::high_resolution_clock::now();
	}

	unsigned int rel_time() {
		return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - prev_timing_clock).count();
	}

	unsigned int abs_time() {
		return std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_timing_clock).count();
	}

	std::string date() {
		auto now = std::chrono::system_clock::now();
		auto in_time_t = std::chrono::system_clock::to_time_t(now);
		std::stringstream ss;
	    ss << std::put_time(std::localtime(&in_time_t), "%d/%m/%Y - %X");
	    return ss.str();
	}
};

#endif
