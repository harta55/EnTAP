/****
DIAMOND protein aligner
Copyright (C) 2013-2017 Benjamin Buchfink <buchfink@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
****/

#include "dp.h"

double background_scores[20];
const double background_freq[] = { 0.0844581,0.0581912,0.0421072,0.0546748,0.0146359,0.040118,0.0621211,0.0669379,0.0225159,0.0547866,0.0957934,0.0523275,0.0218629,0.038769,0.0505311,
0.0760908,0.0573267,0.0127314,0.0295317,0.0644889 };

void init_cbs()
{	
	for (unsigned i = 0; i < 20; ++i) {
		background_scores[i] = 0;
		for (unsigned j = 0; j < 20; ++j)
			background_scores[i] += background_freq[j] * score_matrix(i, j);
	}
}

struct Vector_scores
{
	Vector_scores()
	{
		memset(scores, 0, sizeof(scores));
	}
	Vector_scores& operator+=(Letter l)
	{
		for (unsigned i = 0; i < 20; ++i)
			scores[i] += score_matrix(l, i);
		return *this;
	}
	Vector_scores& operator-=(Letter l)
	{
		for (unsigned i = 0; i < 20; ++i)
			scores[i] -= score_matrix(l, i);
		return *this;
	}
	int scores[20];
};

Bias_correction::Bias_correction(const sequence &seq):
	vector<float>(seq.length())
{
	Vector_scores scores;
	const unsigned window_half = std::min(20u, (unsigned)seq.length() - 1);
	unsigned n = 0;
	unsigned h = 0, m = 0, t = 0, l = (unsigned)seq.length();
	while (n < window_half && h < l) {
		++n;
		scores += seq[h];
		++h;
	}
	while (n < 41 && h < l) {
		++n;
		scores += seq[h];
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++h;
		++m;
	}
	while (h < l) {
		scores += seq[h];
		scores -= seq[t];
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++h;
		++t;
		++m;
	}
	while (m < l && n>21) {
		--n;
		scores -= seq[t];
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++t;
		++m;
	}
	while (m < l) {
		const Letter r = seq[m];
		if (r < 20)
			this->operator[](m) = (float)background_scores[(int)r] - float(scores.scores[(int)r] - score_matrix(r, r)) / (n - 1);
		++m;
	}
}

int Bias_correction::operator()(const Hsp_data &hsp) const
{
	float s = 0;
	for (Hsp_data::Iterator i = hsp.begin(); i.good(); ++i) {
		switch (i.op()) {
		case op_match:
		case op_substitution:
			s += (*this)[i.query_pos];
		default:
			;
		}
	}
	return (int)s;
}

int Bias_correction::operator()(const Diagonal_segment &d) const
{
	float s = 0;
	const int end = d.query_end();
	for (int i = d.i; i < end; ++i)
		s += (*this)[i];
	return (int)s;
}