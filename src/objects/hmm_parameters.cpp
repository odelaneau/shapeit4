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

#include <objects/hmm_parameters.h>

hmm_parameters::hmm_parameters() {
	ed = 0.0001;
	ee = 0.9999;
}

hmm_parameters::~hmm_parameters() {
	t.clear();
	nt.clear();
}

void hmm_parameters::initialise(variant_map & V, int Neff, int Nhap) {
	t = vector < double > (V.size() - 1, 0.0);
	nt = vector < double > (V.size() - 1, 0.0);
	for (int l = 1 ; l < V.size() ; l ++) {
		double dist_cm = V.vec_pos[l]->cm - V.vec_pos[l-1]->cm;
		if (dist_cm <= 0) dist_cm = 0.00001;
		t[l-1] = -1.0 * expm1(-0.04 * Neff * dist_cm / Nhap);
		nt[l-1] = 1-t[l-1];
	}
}

