/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#ifndef _VERBOSE_H
#define _VERBOSE_H

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;

class verbose {
protected:
	ofstream log;
	bool verbose_on_screen;
	bool verbose_on_log;
	int prev_percent;

public:
	verbose() {
		verbose_on_screen = true;
		verbose_on_log = false;
		prev_percent = -1;
	}

	~verbose() {
		close_log();
	}

	bool open_log(string fname) {
		log.open(fname.c_str());
		if (log.fail()) return false;
		else return (verbose_on_log = true);
	}

	void close_log() {
		log.close();
	}

	void set_silent() {
		verbose_on_screen = false;
	}

	void print(string s) {
		if (verbose_on_screen) cout << s << endl;
		if (verbose_on_log) log << s << endl;
	}

	void ctitle(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[32m" << s <<  "\033[0m" << endl;
		if (verbose_on_log) log << endl << s << endl;
	}

	void title(string s) {
		if (verbose_on_screen) cout << endl << s << endl;
		if (verbose_on_log) log << endl << s << endl;
	}

	void bullet(string s) {
		if (verbose_on_screen) cout << "  * " << s << endl;
		if (verbose_on_log) log << "  * " << s << endl;
	}

	void warning(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[33m" << "WARNING: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "WARNING: " << s << endl;
	}

	void leave(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[33m" << "EXITED: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "EXITED: " << s << endl;
		exit(EXIT_SUCCESS);
	}

	void error(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[31m" << "ERROR: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "ERROR: " << s << endl;
		exit(EXIT_FAILURE);
	}

	void done(string s) {
		if (verbose_on_screen) cout << endl << "\x1B[32m" << "DONE: " <<  "\033[0m" << s << endl;
		if (verbose_on_log) log << endl << "DONE: " << s << endl;
		exit(EXIT_SUCCESS);
	}

	void wait(string s) {
		if (verbose_on_screen) {
			cout << s << " ...\r";
			cout.flush();
		}
	}

	void progress(string prefix, float percent) {
		if (verbose_on_screen) {
			int curr_percent = int(percent * 100.0);
			if (prev_percent > curr_percent) prev_percent = -1;
			if (curr_percent > prev_percent) {
				cout << prefix << " [" << curr_percent << "%]\r";
				cout.flush();
				prev_percent = curr_percent;
			}
		}
	}
};
#endif
