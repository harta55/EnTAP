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

#include <limits>
#include <iostream>
#include <sstream>
#include <set>
#include "../basic/config.h"
#include "reference.h"
#include "../basic/statistics.h"
#include "load_seqs.h"
#include "../util/seq_file_format.h"
#include "../util/log_stream.h"
#include "../basic/masking.h"

String_set<0>* ref_ids::data_ = 0;
Ref_map ref_map;
Partitioned_histogram ref_hst;
unsigned current_ref_block;
Reference_header ref_header;
Sequence_set* ref_seqs::data_ = 0;
bool blocked_processing;

using std::cout;
using std::endl;

string* get_allseqids(const char *s)
{
	string *r = new string;
	const vector<string> t(tokenize(s, "\1"));
	for (vector<string>::const_iterator i = t.begin(); i != t.end(); ++i) {
		if (i != t.begin())
			r->append("\1");
		r->append(i->substr(0, find_first_of(i->c_str(), Const::id_delimiters)));
	}
	return r;
}

struct Pos_record
{
	Pos_record()
	{}
	Pos_record(uint64_t pos, size_t len):
		pos(pos),
		seq_len(uint32_t(len))
	{}
	uint64_t pos;
	uint32_t seq_len;
};

void push_seq(const sequence &seq, const sequence &id, uint64_t &offset, vector<Pos_record> &pos_array, Output_stream &out, size_t &letters, size_t &n_seqs)
{	
	pos_array.push_back(Pos_record(offset, seq.length()));
	out.write("\xff", 1);
	out.write(seq.data(), seq.length());
	out.write("\xff", 1);
	out.write(id.data(), id.length() + 1);
	letters += seq.length();
	++n_seqs;
	offset += seq.length() + id.length() + 3;
}

void make_db()
{
	message_stream << "Database file: " << config.input_ref_file << endl;
	
	Timer total;
	total.start();
	if (config.input_ref_file == "")
		std::cerr << "Input file parameter (--in) is missing. Input will be read from stdin." << endl;
	task_timer timer("Opening the database file", true);
	auto_ptr<Input_stream> db_file (Compressed_istream::auto_detect(config.input_ref_file));
	
	Output_stream out(config.database);
	out.write(&ref_header, 1);

	size_t letters = 0, n = 0, n_seqs = 0;
	uint64_t offset = sizeof(ref_header);
	Sequence_set *seqs;
	String_set<0> *ids;
	const FASTA_format format;
	vector<Pos_record> pos_array;

	try {
		while ((timer.go("Loading sequences"), n = load_seqs(*db_file, format, &seqs, ids, 0, (size_t)(1e9), string())) > 0) {
			if (config.masking == 1) {
				timer.go("Masking sequences");
				mask_seqs(*seqs, Masking::get(), false);
			}
			timer.go("Writing sequences");
			for (size_t i = 0; i < n; ++i) {
				sequence seq = (*seqs)[i];
				if (seq.length() == 0)
					throw std::runtime_error("File format error: sequence of length 0 at line " + to_string(db_file->line_count));
				push_seq(seq, (*ids)[i], offset, pos_array, out, letters, n_seqs);
			}
			delete seqs;
			delete ids;
		}
	}
	catch (std::exception&) {
		out.close();
		out.remove();
		throw;
	}
	
	timer.go("Writing trailer");
	ref_header.pos_array_offset = offset;
	pos_array.push_back(Pos_record(offset, 0));
	out.write(pos_array);

	timer.go("Closing the input file");
	db_file->close();
	
	timer.go("Closing the database file");
	ref_header.letters = letters;
	ref_header.sequences = n_seqs;
	out.seekp(0);
	out.write(&ref_header, 1);
	out.close();

	timer.finish();
	message_stream << "Processed " << n_seqs << " sequences, " << letters << " letters." << endl;
	message_stream << "Total time = " << total.getElapsedTimeInSec() << "s" << endl;
}

bool Database_file::load_seqs()
{
	task_timer timer("Loading reference sequences");
	const size_t max_letters = (size_t)(config.chunk_size*1e9);
	seek(pos_array_offset);
	size_t letters = 0, seqs = 0, id_letters = 0;

	ref_seqs::data_ = new Sequence_set;
	ref_ids::data_ = new String_set<0>;

	Pos_record r;
	read(&r, 1);
	size_t start_offset = r.pos;

	while (r.seq_len > 0 && letters < max_letters) {
		Pos_record r_next;
		read(&r_next, 1);
		letters += r.seq_len;
		ref_seqs::data_->reserve(r.seq_len);
		const size_t id = r_next.pos - r.pos - r.seq_len - 3;
		id_letters += id;
		ref_ids::data_->reserve(id);
		pos_array_offset += sizeof(Pos_record);
		++seqs;
		r = r_next;
	}

	if (seqs == 0) {
		delete ref_seqs::data_;
		delete ref_ids::data_;
		return false;
	}

	ref_seqs::data_->finish_reserve();
	ref_ids::data_->finish_reserve();
	seek(start_offset);
	size_t masked = 0;

	for (size_t n = 0; n < seqs; ++n) {
		read(ref_seqs::data_->ptr(n) - 1, ref_seqs::data_->length(n) + 2);
		read(ref_ids::data_->ptr(n), ref_ids::data_->length(n) + 1);
		if (config.masking == 1)
			Masking::get().bit_to_hard_mask(ref_seqs::data_->ptr(n), ref_seqs::data_->length(n), masked);
		else
			Masking::get().remove_bit_mask(ref_seqs::data_->ptr(n), ref_seqs::data_->length(n));
		if (!config.sfilt.empty() && strstr(ref_ids::get()[n].c_str(), config.sfilt.c_str()) == 0)
			memset(ref_seqs::data_->ptr(n), value_traits.mask_char, ref_seqs::data_->length(n));
	}
	timer.finish();
	ref_seqs::get().print_stats();
	log_stream << "Masked letters = " << masked << endl;	

	blocked_processing = seqs < ref_header.sequences;
	return true;
}

void Database_file::get_seq()
{
	vector<Letter> seq;
	string id;
	char c;
	bool all = config.seq_no.size() == 0;
	std::set<size_t> seqs;
	if (!all)
		for (vector<string>::const_iterator i = config.seq_no.begin(); i != config.seq_no.end(); ++i)
			seqs.insert(atoi(i->c_str()) - 1);
	for (size_t n = 0; n < ref_header.sequences; ++n) {
		read(&c, 1);
		while (read(&c, 1), c != '\xff')
			seq.push_back(c);
		while (read(&c, 1), c != '\0')
			id.push_back(c);
		if (all || seqs.find(n) != seqs.end()) {
			cout << '>' << id << endl;
			if (config.reverse) {
				sequence(seq).print(cout, value_traits, sequence::Reversed());
				cout << endl;
			}
			else
				cout << sequence(seq) << endl;
		}
		seq.clear();
		id.clear();
	}
}