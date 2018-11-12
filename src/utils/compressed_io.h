#ifndef _COMPRESSED_IO_H
#define _COMPRESSED_IO_H

//STL INCLUDES
#include <iostream>
#include <sstream>
#include <fstream>

//BOOST INCLUDES
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

class input_file : public boost::iostreams::filtering_istream {
protected:
	std::ifstream file_descriptor;

public:
	input_file(std::string filename) {
		if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
			file_descriptor.open(filename.c_str(), std::ios::in | std::ios::binary);
			push(boost::iostreams::gzip_decompressor());
		} else if (filename.substr(filename.find_last_of(".") + 1) == "bz2") {
			file_descriptor.open(filename.c_str(), std::ios::in | std::ios::binary);
			push(boost::iostreams::bzip2_decompressor());
		} else file_descriptor.open(filename.c_str());
		if (!file_descriptor.fail()) push(file_descriptor);
	}

	~input_file() {
		close();
	}

	bool fail() {
		return file_descriptor.fail();
	}

	void close() {
		if (!file_descriptor.fail()) {
			if (!empty()) reset();
			file_descriptor.close();
		}
	}
};

class output_file : public boost::iostreams::filtering_ostream {
protected:
	std::ofstream file_descriptor;

public:
	output_file(std::string filename) {
		if (filename.substr(filename.find_last_of(".") + 1) == "gz") {
			file_descriptor.open(filename.c_str(), std::ios::out | std::ios::binary);
			push(boost::iostreams::gzip_compressor());
		} else if (filename.substr(filename.find_last_of(".") + 1) == "bz2") {
			file_descriptor.open(filename.c_str(), std::ios::out | std::ios::binary);
			push(boost::iostreams::bzip2_compressor());
		} else file_descriptor.open(filename.c_str());
		if (!file_descriptor.fail()) push(file_descriptor);
	}

	~output_file() {
		close();
	}

	bool fail() {
		return file_descriptor.fail();
	}

	void close() {
		if (!file_descriptor.fail()) {
			if (!empty()) reset();
			file_descriptor.close();
		}
	}
};

#endif
