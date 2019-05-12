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

#include <iostream>
#include <limits>
#include "../data/reference.h"
#include "../data/queries.h"
#include "../basic/statistics.h"
#include "../basic/shape_config.h"
#include "../search/align_range.h"
#include "../util/seq_file_format.h"
#include "../data/load_seqs.h"
#include "../output/output_format.h"
#include "../data/frequent_seeds.h"
#include "../output/daa_write.h"
#include "../data/taxonomy.h"
#include "../basic/masking.h"

using std::endl;
using std::cout;

struct Search_context
{
	Search_context(unsigned sid, const sorted_list &ref_idx, const sorted_list &query_idx) :
		sid(sid),
		ref_idx(ref_idx),
		query_idx(query_idx)
	{ }
	void operator()(unsigned thread_id, unsigned seedp) const
	{
		Statistics stat;
		align_partition(seedp,
			stat,
			sid,
			ref_idx.get_partition_cbegin(seedp),
			query_idx.get_partition_cbegin(seedp),
			thread_id);
		statistics += stat;
	}
	const unsigned sid;
	const sorted_list &ref_idx;
	const sorted_list &query_idx;
};

void process_shape(unsigned sid,
	unsigned query_chunk,
	char *query_buffer,
	char *ref_buffer)
{
	using std::vector;

	::partition<unsigned> p(Const::seedp, config.lowmem);
	for (unsigned chunk = 0; chunk < p.parts; ++chunk) {

		message_stream << "Processing query chunk " << query_chunk << ", reference chunk " << current_ref_block << ", shape " << sid << ", index chunk " << chunk << '.' << endl;
		const seedp_range range(p.getMin(chunk), p.getMax(chunk));
		current_range = range;

		task_timer timer("Building reference index", true);
		sorted_list *ref_idx;
		if (config.algo == Config::query_indexed)
			ref_idx = new sorted_list(ref_buffer,
				*ref_seqs::data_,
				sid,
				ref_hst.get(sid),
				range,
				ref_hst.partition(),
				query_seeds);
		else if (query_seeds_hashed != 0)
			ref_idx = new sorted_list(ref_buffer,
				*ref_seqs::data_,
				sid,
				ref_hst.get(sid),
				range,
				ref_hst.partition(),
				query_seeds_hashed);
		else
			ref_idx = new sorted_list(ref_buffer,
				*ref_seqs::data_,
				sid,
				ref_hst.get(sid),
				range,
				ref_hst.partition(),
				&no_filter);

		timer.go("Building query index");
		sorted_list query_idx(query_buffer,
			*query_seqs::data_,
			sid,
			query_hst.get(sid),
			range,
			query_hst.partition(),
			&no_filter);

		timer.go("Building seed filter");
		frequent_seeds.build(sid, range, *ref_idx, query_idx);

		timer.go("Searching alignments");
		Search_context context(sid, *ref_idx, query_idx);
		launch_scheduled_thread_pool(context, Const::seedp, config.threads_);
		delete ref_idx;
	}
}

void run_ref_chunk(Database_file &db_file,
	Timer &total_timer,
	unsigned query_chunk,
	pair<size_t, size_t> query_len_bounds,
	char *query_buffer,
	Output_stream &master_out,
	vector<Temp_file> &tmp_file)
{
	task_timer timer("Building reference histograms");
	if(config.algo==Config::query_indexed)
		ref_hst = Partitioned_histogram(*ref_seqs::data_, false, query_seeds);
	else if(query_seeds_hashed != 0)
		ref_hst = Partitioned_histogram(*ref_seqs::data_, true, query_seeds_hashed);
	else
		ref_hst = Partitioned_histogram(*ref_seqs::data_, false, &no_filter);

	ref_map.init(safe_cast<unsigned>(ref_seqs::get().get_length()));

	timer.go("Allocating buffers");
	char *ref_buffer = sorted_list::alloc_buffer(ref_hst);

	timer.go("Initializing temporary storage");
	Trace_pt_buffer::instance = new Trace_pt_buffer(query_seqs::data_->get_length() / align_mode.query_contexts,
		config.tmpdir,
		config.query_bins);
	timer.finish();
	
	for (unsigned i = 0; i < shapes.count(); ++i)
		process_shape(i, query_chunk, query_buffer, ref_buffer);

	timer.go("Deallocating buffers");
	delete[] ref_buffer;

	Output_stream* out;
	if (blocked_processing) {
		timer.go("Opening temporary output file");
		tmp_file.push_back(Temp_file());
		out = new Output_stream(tmp_file.back());
	}
	else
		out = &master_out;

	timer.go("Computing alignments");
	align_queries(*Trace_pt_buffer::instance, out);
	delete Trace_pt_buffer::instance;

	if (blocked_processing) {
		Intermediate_record::finish_file(*out);
		delete out;
	}

	timer.go("Deallocating reference");
	delete ref_seqs::data_;
	delete ref_ids::data_;
	timer.finish();
}

void run_query_chunk(Database_file &db_file,
	Timer &total_timer,
	unsigned query_chunk,
	Output_stream &master_out,
	Output_stream *unaligned_file)
{
	static const double max_coverage = 0.15;

	task_timer timer("Building query seed set");
	if (query_chunk == 0)
		setup_search_cont();
	if (config.algo == -1) {
		query_seeds = new Seed_set(query_seqs::get(), max_coverage);
		timer.finish();
		log_stream << "Seed space coverage = " << query_seeds->coverage() << endl;
		if (query_seeds->coverage() >= max_coverage) {
			config.algo = Config::double_indexed;
			delete query_seeds;
			query_seeds = 0;
		}
		else
			config.algo = Config::query_indexed;
	}
	else if (config.algo == Config::query_indexed) {
		query_seeds = new Seed_set(query_seqs::get(), 2);
		timer.finish();
		log_stream << "Seed space coverage = " << query_seeds->coverage() << endl;
	}
	else
		timer.finish();
	if (query_chunk == 0)
		setup_search();
	if (config.algo == Config::double_indexed && config.small_query) {
		timer.go("Building query seed hash set");
		query_seeds_hashed = new Hashed_seed_set(query_seqs::get());
	}

	timer.go("Building query histograms");
	const pair<size_t, size_t> query_len_bounds = query_seqs::data_->len_bounds(shapes[0].length_);
	setup_search_params(query_len_bounds, 0);
	query_hst = Partitioned_histogram(*query_seqs::data_, false, &no_filter);
	timer.finish();
	//const bool long_addressing_query = query_seqs::data_->raw_len() > (size_t)std::numeric_limits<uint32_t>::max();

	timer.go("Allocating buffers");
	char *query_buffer = sorted_list::alloc_buffer(query_hst);
	vector<Temp_file> tmp_file;
	query_aligned.clear();
	query_aligned.insert(query_aligned.end(), query_ids::get().get_length(), false);
	db_file.rewind();
	timer.finish();
	
	for (current_ref_block = 0; db_file.load_seqs(); ++current_ref_block)
		run_ref_chunk(db_file, total_timer, query_chunk, query_len_bounds, query_buffer, master_out, tmp_file);

	timer.go("Deallocating buffers");
	delete[] query_buffer;
	delete query_seeds;
	query_seeds = 0;

	if (blocked_processing) {
		timer.go("Joining output blocks");
		join_blocks(current_ref_block, master_out, tmp_file);
	}

	if (unaligned_file) {
		timer.go("Writing unaligned queries");
		write_unaligned(unaligned_file);
	}

	timer.go("Deallocating queries");
	delete query_seqs::data_;
	delete query_ids::data_;
	delete query_source_seqs::data_;
}

void master_thread(Database_file &db_file, Timer &total_timer)
{
	if(config.query_file.empty())
		std::cerr << "Query file parameter (--query/-q) is missing. Input will be read from stdin." << endl;
	task_timer timer("Opening the input file", true);
	auto_ptr<Input_stream> query_file(Compressed_istream::auto_detect(config.query_file));
	const Sequence_file_format *format_n(guess_format(*query_file));

	current_query_chunk = 0;

	timer.go("Opening the output file");
	auto_ptr<Output_stream> master_out(config.compression == 1
		? new Compressed_ostream(config.output_file)
		: new Output_stream(config.output_file));
	if (*output_format == Output_format::daa)
		init_daa(*master_out);
	auto_ptr<Output_stream> unaligned_file;
	if (!config.unaligned.empty())
		unaligned_file = auto_ptr<Output_stream>(new Output_stream(config.unaligned));
	timer.finish();

	for (;; ++current_query_chunk) {
		task_timer timer("Loading query sequences", true);
		size_t n_query_seqs;
		n_query_seqs = load_seqs(*query_file, *format_n, &query_seqs::data_, query_ids::data_, &query_source_seqs::data_, (size_t)(config.chunk_size * 1e9), config.qfilt);
		if (n_query_seqs == 0)
			break;
		timer.finish();
		query_seqs::data_->print_stats();

		if (current_query_chunk == 0 && *output_format != Output_format::daa)
			output_format->print_header(*master_out, align_mode.mode, config.matrix.c_str(), score_matrix.gap_open(), score_matrix.gap_extend(), config.max_evalue, query_ids::get()[0].c_str(),
				unsigned(align_mode.query_translated ? query_source_seqs::get()[0].length() : query_seqs::get()[0].length()));

		if (config.masking == 1) {
			timer.go("Masking queries");
			mask_seqs(*query_seqs::data_, Masking::get());
			timer.finish();
		}		

		run_query_chunk(db_file, total_timer, current_query_chunk, *master_out, unaligned_file.get());
	}

	timer.go("Closing the input file");
	query_file->close();

	timer.go("Closing the output file");
	if (*output_format == Output_format::daa)
		finish_daa(*master_out);
	else
		output_format->print_footer(*master_out);
	master_out->close();
	if (unaligned_file.get())
		unaligned_file->close();
	
	timer.go("Closing the database file");
	db_file.close();

	timer.finish();
	message_stream << "Total time = " << total_timer.getElapsedTimeInSec() << "s" << endl;
	statistics.print();
}

void master_thread_di()
{
	Timer timer2;
	timer2.start();

	align_mode = Align_mode(Align_mode::from_command(config.command));
	init_output();

	message_stream << "Temporary directory: " << Temp_file::get_temp_dir() << endl;

	if (config.mode_very_sensitive) {
		Config::set_option(config.chunk_size, 0.4);
		Config::set_option(config.lowmem, 1u);
	}
	else {
		Config::set_option(config.chunk_size, 2.0);
		Config::set_option(config.lowmem, 4u);
	}

	task_timer timer("Opening the database", 1);
	Database_file db_file;
	timer.finish();
	verbose_stream << "Reference = " << config.database << endl;
	verbose_stream << "Sequences = " << ref_header.sequences << endl;
	verbose_stream << "Letters = " << ref_header.letters << endl;
	verbose_stream << "Block size = " << (size_t)(config.chunk_size * 1e9) << endl;
	Config::set_option(config.db_size, (uint64_t)ref_header.letters);

	set_max_open_files(config.query_bins * config.threads_ + unsigned(ref_header.letters / (size_t)(config.chunk_size * 1e9)) + 16);
	taxonomy.init();

	master_thread(db_file, timer2);
}