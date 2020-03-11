////////////////////////////////////////////////////////////////////////////////
// Copyright (C) 2018 Olivier Delaneau, University of Lausanne
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
////////////////////////////////////////////////////////////////////////////////

#include <containers/variant_map.h>

variant_map::variant_map() {
}

variant_map::~variant_map() {
	for (int s = 0 ; s < vec_pos.size() ; s++) delete vec_pos[s];
	vec_pos.clear();
	map_pos.clear();
}

int variant_map::size() {
	return vec_pos.size();
}

variant * variant_map::getByIndex (int i) {
	return vec_pos[i];
}

vector < variant * > variant_map::getByPos (int pos) {
	vector < variant * > vecS;
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}

vector < variant * > variant_map::getByRef(int pos, string & ref, string & alt) {
	vector < variant * > vecS = vector < variant * >();
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) {
		if (it->second->ref == ref && it->second->alt == alt) vecS.push_back(it->second);
	}
	return vecS;
}

void variant_map::push(variant * v) {
	vec_pos.push_back(v);
	map_pos.insert(pair < int , variant * > (v->bp, v));
}


unsigned int variant_map::length() {
	return vec_pos.back()->bp - vec_pos[0]->bp + 1;
}
